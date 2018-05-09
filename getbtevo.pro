pro getbtevo, obslist, evolist, output, $
              NITER = niter, CIRCULARIZE = circularize

  if NOT keyword_set(CIRCULARIZE) then circularize = 0 else circularize = 1
  
  if NOT keyword_set(NITER) then niter = 20
  if n_elements(output) eq 0 then suffix = '_lgn_structEvo.fits' $
  else suffix = output
  
  summaryData = mrdfits('objectSummaryData.fits', 1)
  finalMasses = mrdfits('finalMasses.fits', 1)

  readcol, obslist, obsfiles, f = 'A'
  readcol, evolist, evofiles, f = 'A'
  ngals = n_elements(obsfiles)

  files = strarr(ngals)

  test     = mrdfits(obsfiles[0], 1)
  regions  = test.REGION
;  nregions = n_elements(regions)
  
  extract = where(regions eq 'INNFL' OR $
                  regions eq 'INTER' OR $
                  regions eq 'OUTER')
  master  = where(regions eq 'OPTFL')
  regions = regions[[extract,master]]
  nregions = n_elements(extract) + 1 ;; INN INT OUT OPT

  test = mrdfits(evofiles[0], 1)
  eextract = where(test.REGION eq 'INNFL' OR $
                   test.REGION eq 'INTER' OR $
                   test.REGION eq 'OUTER')
  emaster  = where(test.REGION eq 'OPTFL')

  dr = 0.02
  rad = findgen(2./dr+1) * dr
  for ii = 0, ngals - 1 do begin
     
     files[ii] = strmid(obsfiles[ii], 0, 7)
     q = where(summaryData.ID eq files[ii])
     mag = summaryData[q].MU
     emag = summaryData[q].EMU / mag / alog(10)
     
     odat   = mrdfits(obsfiles[ii], 1)
     evodat = mrdfits(evofiles[ii], 1)

     odat   = odat[[extract,master]]
     evodat = evodat[[eextract,emaster]]
     
     time = evodat[emaster].TIME

     odex = value_locate(time, evodat[-1].TOBS)
     print, odat[-1].LMASS[0] - alog10(mag)
     print, alog10(evodat[-1].MGH[odex,1]) - alog10(mag)
     print, alog10(total(evodat[0:2].MGH[odex,1] * [1,2,2])) - alog10(mag), $
            sqrt(total((0.5/alog(10) * (evodat[0:2].MGH[odex,2] - evodat[0:2].MGH[odex,0]) * sqrt([1,2,2]) / evodat[0:2].MGH[odex,1])^2) + emag^2)
     print, n_elements(time)
     
     savedata = {REGION     : string(0), $
                 TIME       : time, $
                 REDSHIFT   : evodat[-1].REDSHIFT, $
                 TOBS       : evodat[-1].TOBS, $
                 OMASS      : fltarr(n_elements(TIME), 4), $
                 OSFR       : fltarr(n_elements(TIME), 4), $
                 TMASS      : fltarr(n_elements(TIME), 4), $
                 TSFR       : fltarr(n_elements(TIME), 4), $
                 SIGMA_SFR  : fltarr(n_elements(TIME), 2), $
                 SIGMA_MSTEL: fltarr(n_elements(TIME), 2), $
                 BT         : fltarr(n_elements(TIME), 2), $
                 HM_RAD     : fltarr(n_elements(TIME), 2), $
                 BT_OBS     : [0.,0.], $
                 HM_RAD_OBS : [0.,0.], $
                 MPROF_INT  : fltarr(n_elements(TIME), 2)}
     savedata = replicate(savedata, nregions)
     savedata.REGION = regions
     
     ;; Do optimal + total histories first
     savedata[-1].OMASS[*,0] = alog10(evodat[-1].MGH[*,1] / mag) ;; This is the median...
     savedata[-1].OMASS[*,2] = alog10(evodat[-1].MGH[*,0] / mag)
     savedata[-1].OMASS[*,3] = alog10(evodat[-1].MGH[*,2] / mag)
     
     err = 0.5/alog(10) * (evodat[-1].MGH[*,2] - evodat[-1].MGH[*,0]) $
           / evodat[-1].MGH[*,1]
;     err = 0.5 * alog10(evodat[-1].MGH[*,2]/evodat[-1].MGH[*,0])
     err = sqrt(err^2 + emag^2)
     savedata[-1].OMASS[*,1] = err

     savedata[-1].OSFR[*,0]  = alog10(evodat[-1].SFH[*,1] / mag) ;; This is the median...
     savedata[-1].OSFR[*,2]  = alog10(evodat[-1].SFH[*,0] / mag)
     savedata[-1].OSFR[*,3]  = alog10(evodat[-1].SFH[*,2] / mag)
     
     err = 0.5/alog(10) * (evodat[-1].SFH[*,2] - evodat[-1].SFH[*,0]) / evodat[-1].SFH[*,1]
;     err = 0.5 * alog10(evodat[-1].SFH[*,2]/evodat[-1].SFH[*,0])
     err = sqrt(err^2 + emag^2)
     savedata[-1].OSFR[*,1]  = err

     for jj = 0, n_elements(time) - 1 do begin
        savedata[*].TMASS[jj,0]  = alog10(total(evodat[0:2].MGH[jj,1] * [1,2,2]) / mag)
        savedata[*].TMASS[jj,2]  = alog10(total(evodat[0:2].MGH[jj,0] * [1,2,2]) / mag)
        savedata[*].TMASS[jj,3]  = alog10(total(evodat[0:2].MGH[jj,2] * [1,2,2]) / mag)
        
        err = 0.5/alog(10) * (evodat[0:2].MGH[jj,2] $
                              - evodat[0:2].MGH[jj,0]) * sqrt([1,2,2]) / $
              evodat[0:2].MGH[jj,1]  ;; LEA 20170823. 1) Caught bug in
;        this line; 2) to avoid log-error blow-ups at low SFRs,
;        switching to slightly different error definition. (see
;        plotmasscomptime.pro on LEA's personal laptop). No
;        qualitative changes.
;        err = 0.5 * alog10(evodat[0:2].MGH[jj,2]/evodat[0:2].MGH[jj,0]) * sqrt([1,2,2])
        savedata[*].TMASS[jj,1]  = sqrt(total(err^2) + emag^2)

        savedata[*].TSFR[jj,0]   = alog10(total(evodat[0:2].SFH[jj,1] * [1,2,2] / mag))
        savedata[*].TSFR[jj,2]   = alog10(total(evodat[0:2].SFH[jj,0] * [1,2,2] / mag))
        savedata[*].TSFR[jj,3]   = alog10(total(evodat[0:2].SFH[jj,2] * [1,2,2] / mag))
        
        err = 0.5/alog(10) * (evodat[0:2].SFH[jj,2] - $
                              evodat[0:2].SFH[jj,0]) * sqrt([1,2,2]) / $
              evodat[0:2].SFH[jj,1]
;        err = 0.5 * alog10(evodat[0:2].SFH[jj,2]/evodat[0:2].SFH[jj,0]) * sqrt([1,2,2])
        savedata[*].TSFR[jj,1]   = sqrt(total(err^2) + emag^2)
     endfor
        
     ;; Do regions
     for jj = 0, n_elements(regions) - 2 do begin
        savedata[jj].SIGMA_SFR[*,0]   = alog10(reform(evodat[jj].SFH[*,1]) / odat[jj].AREA_PHYS)
        err = 0.5/alog(10) * (evodat[jj].SFH[*,2] - evodat[jj].SFH[*,0]) / evodat[jj].SFH[*,1]
;        err = sqrt(err^2 + emag^2) ;; area subject to uncertainty through mag factor
        savedata[jj].SIGMA_SFR[*,1]   = err
        savedata[jj].SIGMA_MSTEL[*,0] = alog10(reform(evodat[jj].MGH[*,1]) / odat[jj].AREA_PHYS)
        err = 0.5/alog(10) * (evodat[jj].MGH[*,2] - evodat[jj].MGH[*,0]) / evodat[jj].MGH[*,1]
;        err = sqrt(err^2 + emag^2) ;; area subject to uncertainty through mag factor
        savedata[jj].SIGMA_MSTEL[*,1] = err
     endfor
     
     ;; Do B/T, etc
     tdat   = odat[0:-2]
     mmaster  = evodat[-1].MGH[*,1] / mag
     emmaster = 0.5 * (evodat[-1].MGH[*,2] - evodat[-1].MGH[*,0]) / mag
     femmaster = emmaster / mmaster
     
     if NOT circularize then begin
;        m    = alog10(evodat[eextract].MGH[*,1]) ;; FOLDED -- must multiply inter and outer by 2
;        em   = 0.5/alog(10) * (evodat[eextract].MGH[*,2] - evodat[eextract].MGH[*,0]) / 10.^m
     endif else begin
        m      = savedata[0:-2].SIGMA_MSTEL[*,0]
        em     = savedata[0:-2].SIGMA_MSTEL[*,1]
     endelse
        
     tr     = tdat.RNORM
     rp     = tdat[0].RE_PHYS     
     for kk = 0, n_elements(evodat[-1].TIME) - 1 do begin
        mprofs = fltarr(n_elements(rad), niter)
        tbts   = fltarr(niter,2)
        thmrs  = fltarr(niter,2)
        for ll = 0, niter - 1 do begin          
           if NOT circularize then begin
              tm = (m[kk,*] + alog10([1,2,2])) + randomn(seed, 3) * (em[kk,*] + alog10([1,sqrt(2), sqrt(2)]))
              t = total(10.^tm, /cum) ;
              mprofs[*,ll] = interpol(t, tr, rad)                     ;
              tmp1 = mprofs[0,ll] / mprofs[value_locate(rad, tr[-1]),ll] ;10.^odat[-1].LMASS[0];mprofs[-1,ll]
              tmp2 = mprofs[0,ll] / (mmaster[kk] * (1 + randomn(seed, 1) * emmaster))[0]              
              tbts[ll,*]  = [tmp1,tmp2] ;[tmp1[value_locate(rad, 1.0)], tmp2[value_locate(rad, 1.0)]]
              thmrs[ll,*] = [rad[value_locate(mprofs[*,ll], 0.5 * mprofs[-1,ll])], $
                             rad[value_locate(mprofs[*,ll], 0.5 * mmaster[kk])]] * rp / sqrt(mag)
;           plot, rad, tmp & wait, 0.001
           endif else begin
              tm = m[kk,*] + randomn(seed, 3) * em[kk,*]
              t = interpol(tm, tr, rad)
              mprofs[*,ll] = total(2 * !pi * rad * dr * rp^2 * 10.^t, /cum)
              mprofs[*,ll] /= mprofs[-1,ll] / (mmaster[kk] * (1 + randomn(seed, 1) * emmaster))[0]              
              tmp1 = mprofs[value_locate(rad*rp,2.5),ll] / mprofs[-1,ll]
              tmp2 = mprofs[value_locate(rad*rp,2.5),ll] / 10.^(mmaster + randomn(seed, 1) * emmaster)
              tbts[ll,*]  = [tmp1,tmp2] ;[tmp1[value_locate(rad, 1.0)], tmp2[value_locate(rad, 1.0)]]
              thmrs[ll,*] = [rad[value_locate(mprofs[*,ll], 0.5 * mprofs[-1,ll])], $
                             rad[value_locate(mprofs[*,ll], 0.5 * mmaster[kk])]] * rp / sqrt(mag)
           endelse
        endfor
        savedata[*].BT[kk,0] = median(tbts[*,0]); + median(tbts[*,1]))
        savedata[*].BT[kk,1] = stddev(tbts[*,0],/nan)
        savedata[*].HM_RAD[kk,0] = median(thmrs[*,0]); + median(thmrs[*,1]))
        savedata[*].HM_RAD[kk,1] = stddev(thmrs[*,0],/nan)
        savedata[*].MPROF_INT[kk,0] = median(mprofs[value_locate(rad, tr[-1]),*])
        savedata[*].MPROF_INT[kk,1] = stddev(mprofs[value_locate(rad, tr[-1]),*], /nan)
     endfor

     if NOT circularize then begin
        m   = tdat.LMASS[0] + [0, alog10(2), alog10(2)] ;; Folded; must double inter & outer stuff
        em  = 0.5 * (tdat.LMASS[2] - tdat.LMASS[1]) + [0, alog10(sqrt(2)), alog10(sqrt(2))]
     endif else begin
        m    = tdat.LMASS[0] - alog10(tdat.AREA_PHYS)
        em   = 0.5 * (tdat.LMASS[2] - tdat.LMASS[1])
     endelse
     mmaster = odat[-1].LMASS[0]
     mprofs = fltarr(n_elements(rad), niter)
     tbts   = fltarr(niter, 2)
     thmrs  = fltarr(niter, 2)
     for ll = 0, niter - 1 do begin
        tm = m + randomn(seed, n_elements(em)) * em
        if NOT circularize then begin
           t = total(10.^tm, /cum) ; interpol(tm, tr, rad)     
           mprofs[*,ll] = interpol(t, tr, rad)
           tmp1 = mprofs[value_locate(rad,0),ll] / mprofs[value_locate(rad,tr[-1]),ll]
           tmp2 = mprofs[value_locate(rad,0),ll] / 10.^mmaster
           tbts[ll,*] = [tmp1, tmp2]
           thmrs[ll,*] = [rad[value_locate(mprofs[*,ll] / mprofs[-1,ll], 0.5)], $
                          rad[value_locate(mprofs[*,ll] / 10.^mmaster, 0.5)]] * rp / sqrt(mag)
        endif else begin
           t = interpol(tm, tr, rad)
           mprofs[*,ll] = total(2 * !pi * rad * dr * rp^2 * 10.^t, /cum) ;; MAG INDEPENDENT
           mprofs[*,ll] /= mprofs[-1,ll]/10.^(mmaster + randomn(seed, 1) * emmaster)[0]
           tmp1 = mprofs[value_locate(rad*rp,2.5),ll] / mprofs[-1,ll]
           tmp2 = mprofs[value_locate(rad*rp,2.5),ll] / 10.^(mmaster + randomn(seed, 1) * emmaster)[0]
           tbts[ll,*] = [tmp1, tmp2]
           thmrs[ll,*] = [rad[value_locate(mprofs[*,ll] / mprofs[-1,ll], 0.5)], $
                          rad[value_locate(mprofs[*,ll] / 10.^mmaster, 0.5)]] * rp / sqrt(mag)
        endelse
     endfor
     savedata[*].BT_OBS[0] = median(tbts[*,0]); + median(tbts[*,1]))
     savedata[*].BT_OBS[1] = stddev(tbts[*,0],/nan)
     savedata[*].HM_RAD_OBS[0] = median(thmrs[*,0]); + median(thmrs[*,1]))
     savedata[*].HM_RAD_OBS[1] = stddev(thmrs[*,0],/nan)
     
     mwrfits, savedata, files[ii]+suffix, /create
     
  endfor
  
end
;getbtevo, 'study5_lgn_Bestfits.list', 'study5_lgn_Recons.list', niter
;= 100

pro doCirc
  getbtevo, 'study5_lgn_Bestfits.list', 'study5_lgn_Recons.list', $
            niter = 100, /circularize
end
