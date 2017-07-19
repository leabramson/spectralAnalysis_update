;; Plot the observed profiles of the galaxies in any quantity
pro plotprofiles, lgnlist, explist, output, $
                  QUANTITY = quantity, $
                  NITER = niter

  if NOT keyword_set(NITER) then niter = 100
  
  summaryData = mrdfits('objectSummaryData.fits', 1)
  finalMasses = mrdfits('finalMasses.fits', 1)

  readcol, explist, expfiles, f = 'A'
  readcol, lgnlist, lgnfiles, f = 'A'
  ngals = n_elements(lgnfiles)
  files = strarr(ngals)

  test     = mrdfits(lgnfiles[0], 1)
  regions  = test.REGION
  nregions = n_elements(regions)
  
  extract = where(regions eq 'OPTFL' OR $
                  regions eq 'INNFL' OR $
                  regions eq 'INTER' OR $
                  regions eq 'OUTER')
  master  = where(regions eq 'OPTFL')
  
  sizes   = fltarr(ngals, 2)
  esizes  = fltarr(ngals, 2)
  ssizes  = fltarr(ngals, 2)
  essizes = fltarr(ngals, 2)
  areas   = fltarr(ngals, nregions, 2)
  eareas  = fltarr(ngals, nregions, 2)
  masses  = fltarr(ngals, nregions, 2)
  emasses = fltarr(ngals, nregions, 2)
  avs     = fltarr(ngals, nregions, 2)
  eavs    = fltarr(ngals, nregions, 2)
  sfrs    = fltarr(ngals, nregions, 2)
  esfrs   = fltarr(ngals, nregions, 2)
  lsfrs   = fltarr(ngals, nregions, 2)
  t0s      = fltarr(ngals, nregions)
  et0s     = fltarr(ngals, nregions)
  taus     = fltarr(ngals, nregions)
  etaus    = fltarr(ngals, nregions)
  pktimes  = fltarr(ngals, nregions)
  epktimes = fltarr(ngals, nregions)
  ages     = fltarr(ngals, nregions)
  eages    = fltarr(ngals, nregions)
  efolds   = fltarr(ngals, nregions)
  
  extrads = fltarr(ngals, nregions)
    
  ;; Masses are consistent w/in the errors including magnification
  
  mags = fltarr(ngals)
  emags = fltarr(ngals)
  bts  = fltarr(ngals, 2)
  ebts = fltarr(ngals,2)
  files = strarr(ngals)

  dr = 0.05
  rad = findgen(2/dr+1) * dr
  
  msigmaprofs = fltarr(ngals, n_elements(rad), niter, 2)
  mprofs = fltarr(ngals, n_elements(rad), niter, 2)
  dustprofs = fltarr(ngals, n_elements(rad), niter, 2)
  hmrs = fltarr(ngals, niter, 2)
  prs  = fltarr(ngals, n_elements(rad))

  nsig = 2.
  for ii = 0, ngals - 1 do begin
     for jj = 0, 1 do begin
        case jj of
           0: begin
              file = lgnfiles[ii]
              files[ii] = strmid(file, 0, 7)
              qui = where(summaryData.ID eq files[ii])
              mag  = summaryData[qui].MU
              emag = 1./alog(10) * summaryData[qui].EMU / mag
              mags[ii] = mag
              emags[ii] = emag
           end
           1: file = expfiles[ii]
        endcase
                
        data = mrdfits(file, 1)

        foo = calcbt(data, mags[ii], emags[ii])
        bts[ii,jj] = foo[0]
        ebts[ii,jj] = foo[1]
        
        ro = data[0].RE_OBS
        rp = data[0].RE_PHYS
        r  = (mrdfits(data[0].FITSFILE, 1, /silent)).RE
        rs = (mrdfits(data[0].FITSFILE, 1, /silent)).RE_SEX
        if jj eq 0 then begin
           extrads[ii,*] = data.RNORM
           t0s[ii,*] = alog(data.T0[0]/1d9)
           et0s[ii,*] = 0.5 * (data.T0[2] - data.T0[1]) / data.T0[0]
           taus[ii,*] = data.TAU[0]
           etaus[ii,*] = 0.5 * (data.Tau[2] - data.Tau[1])
           pktimes[ii,*] = t0s[ii,*] - taus[ii,*]^2
           epktimes[ii,*] = sqrt(et0s[ii,*]^2 + 2 * etaus[ii,*]^2)
        endif else if jj eq 1 then begin
           ages[ii,*] = data.AGE[0]/1d9
           eages[ii,*] = 0.5 * (data.AGE[2] - data.AGE[1]) / 1d9
           efolds[ii,*] = data.AGE[0] / data.TAU[0]
        endif
           
        prs[ii,*] = rp * rad
        
        pscale = ro / r
        kscale = rp / r
        
        ;; Correct for magnification
        data.LMASS -= alog10(mag)
        data.SFR   /= mag
        data.AREA_PHYS /= mag
;        data[where(data.REGION eq 'INNFL')].AREA_PHYS /= 2 ;; Need to account for FOLDING!
        
        masses[ii,*,jj]  = data.LMASS[0]
        emasses[ii,*,jj] = sqrt((0.5 * (data.LMASS[2] - data.LMASS[1]))^2 + emag^2)
        
        areas[ii,*,jj]  = alog10(data.AREA_PHYS)
        eareas[ii,*,jj] = emag

        avs[ii,*,jj]  = data.A_V[0]
        eavs[ii,*,jj] = 0.5 * (data.A_V[2] - data.A_V[1])
        
        sizes[ii,jj]  = alog10(rp / sqrt(mag))
        esizes[ii,jj] = 0.5 * emag

        ssizes[ii,jj]  = alog10(rs * kscale / sqrt(mag))
        essizes[ii,jj] = 0.5 * emag

        sfrs[ii,*,jj]  = data.SFR[0]
        esfrs[ii,*,jj] = 0.5 * (data.SFR[2] - data.SFR[1]) ;; mag errs added below

        uls = where(sfrs[ii,*,jj] lt nsig * (data.SFR[0] - data.SFR[1]), nuls) ;; lower error bar
        if nuls gt 0 then begin
           sfrs[ii,uls,jj] = 0
           lsfrs[ii,uls,jj] = (data.SFR[0] - data.SFR[1])[uls]
           esfrs[ii,uls,jj] = (data.SFR[0] - data.SFR[1])[uls]
        endif
           
        ;; Estimate B/T
        ;; Mass w/in r_e over mass w/in 2 * re
        tdat = data[extract[1:*]]
        tdat = tdat[sort(tdat.RNORM)]
        m    = tdat.LMASS[0] - areas[ii,*,jj]
        em   = sqrt((0.5 * (tdat.LMASS[2] - tdat.LMASS[1]))^2 + eareas[ii,*,jj]^2 + emag^2)
        s    = tdat.SFR[0] - areas[ii,*,jj]
        es   = sqrt((0.5*(tdat.SFR[2] - tdat.SFR[1]))^2 + eareas[ii,*,jj]^2 + emag^2)
        ssfr = alog10(tdat.SFR[0]) - tdat.LMASS[0]
        essfr = sqrt((0.5/alog(10)*(tdat.SFR[2] - tdat.SFR[1])/tdat.SFR[0])^2 + (0.5 * (tdat.LMASS[2] - tdat.LMASS[1]))^2)
        av   = tdat.A_V[0]
        eav  = 0.5 * (tdat.A_V[2] - tdat.A_V[1])
        
        tr   = tdat.RNORM
        for ll = 0, niter - 1 do begin
           tm = m + randomn(seed, 3) * em
           ts = s + randomn(seed, 3) * es
           tssfr = ssfr + randomn(seed, 3) * essfr
           tav = av + randomn(seed, 3) * eav
           
           msigmaprofs[ii,*,ll,jj] = interpol(tm, tr, rad)
           dustprofs[ii,*,ll,jj] = interpol(tav, tr, rad)     

           mprofs[ii,*,ll,jj] = total(2 * !pi * rad * dr * rp^2 * 10.^msigmaprofs[ii,*,ll,jj], /cum) ;; mag independent!
           hmrs[ii,ll,jj] = rad[value_locate(mprofs[ii,*,ll,jj], 0.5 * mprofs[ii,-1,ll,jj])]           
        endfor
        
     endfor
     
  endfor
  tsfrs = mean(sfrs[*,master,*], dim = 3, /nan)
;  tsfrs = sfrs[*,master,0];, dim = 3, /nan)
  tsfrs[where(tsfrs eq 0)] = nsig*esfrs[where(tsfrs eq 0),master,0]
  resreg = extract[value_locate(regions[extract], ['INNFL', 'INTER', 'OUTER'])]
  
  cgloadct, 33, /rev, ncolors = 9, $
            clip = [10,240]

  ssfrs = alog10(tsfrs/10.^masses[*,master,0])
  minssfr = -11
  maxssfr = -8
  dssfr   = (maxssfr - minssfr) / 8.
  tssfrs  = findgen(9)*dssfr + minssfr
  cols  = []
  for ii = 0, ngals - 1 do begin
     foo = min(abs(tssfrs - ssfrs[ii]), hit)
     cols = [cols, hit[0]]
  endfor

  ;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;
  ;; Mass Density ;;
  ;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;
  
  set_plot, 'PS'
  xsize = 8.5
  ysize = 4
  ar = ysize / xsize
  device, filename = 'massProfiles.eps', $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = xsize, ysize = ysize, /in
  plot, rad, median(reform(msigmaprofs[0,*,*,0]), dim = 2), /nodat, $
        xtitle = '!18r/r!D!Xe!N', xran = [-0.1,1.8], $
        yran = [7.5,10], $
        ytitle = 'log '+greek('Sigma')+'!D!18M!X!L*!N [M!D'+sunsymbol()+'!N kpc!E-2!N]', $
        xthick = 3, ythick = 3, charsize = 1.25, charthick = 4, $
        pos = [0.1,0.15,0.5,0.9], /xsty
  xxx = [rad, reverse(rad)]
  plotsym, 0, /fill
  for ii = 0, ngals - 1 do begin
     tre = (median(hmrs[ii,*,0], dim = 2))[0]
     sre = (stddev(hmrs[ii,*,0], dim = 2, /nan))[0]
     oplot, replicate(tre,2), !Y.CRANGE[0] + [0,0.3], thick = 10
     oplot, tre + sre * [-1,1], $
            !Y.CRANGE[0] + [0.3,0.3], thick = 10
     oplot, replicate(tre,2), !Y.CRANGE[0] + [0,0.3], col = cgcolor(string(cols[ii])), thick = 6
     oplot, tre + sre * [-1,1], $
            !Y.CRANGE[0] + [0.3,0.3], thick = 6, col = cgcolor(string(cols[ii]))
     ;; do detections and upper limits...
     tmasses = reform(masses[ii,tr,0])
     etmasses = reform(emasses[ii,tr,0])
     tr = resreg[sort(extrads[ii,resreg])]
     ttr = reform(extrads[ii,tr])
     nuse = n_elements(tr)
     oplot, ttr, tmasses - areas[ii,tr,0], $
            thick = 10
     oplot, ttr, tmasses - areas[ii,tr,0], $
            thick = 4, col = cgcolor(string(cols[ii]))
     oploterror, ttr, tmasses - areas[ii,tr,0], $
                 replicate(emags[ii],nuse), reform(sqrt(etmasses + emags[ii]^2)), $
                 psym = 8, symsize = 1.3, /nohat
     oplot, ttr, tmasses - areas[ii,tr,0], $
            col = cgcolor(string(cols[ii])), psym = 8

;     oplot, extrads[ii,tr], masses[ii,tr,1] - areas[ii,tr,1], $ ;;
;                                              Not meaningfully different
;            thick = 10, linesty = 2, col = '777777'x
;     oplot, extrads[ii,tr], masses[ii,tr,1] - areas[ii,tr,1], $
;            thick = 4, col = cgcolor(string(cols[ii])), linesty = 2
;     oploterror, extrads[ii,tr], masses[ii,tr,1] - areas[ii,tr,1], $
;                 [replicate(emags[ii],nuse)], reform(sqrt(emasses[ii,tr,1] + emags[ii]^2)), $
;                 psym = 8, symsize = 1.3, /nohat
;     oplot, extrads[ii,tr], masses[ii,tr,1] - areas[ii,tr,1], $
;            col = cgcolor(string(cols[ii])), psym = 8
;     mp = reform(msigmaprofs[ii,*,*,0])     
;     mmm = median(mp, dim = 2)
;     eee = stddev(mp, dim = 2)
;     oplot, rad, mmm+eee, thick = 4
;     oplot, rad, mmm-eee, thick = 4     
;     polyfill, xxx, [mmm-eee, reverse(mmm+eee)] > !Y.CRANGE[0], $
;               col = cgcolor(string(cols[ii])), /line_fill, $
;               thick = 1, orien = ii * 45, spacing = 0.025
  endfor
;  cgtext, 1.1, !Y.CRANGE[1]-0.1, 'half-mass radii', $
;          orien = 270, align = 0, charsize = 1, charthick = 3
  
  top = 0.3
  plot, rad, median(reform(msigmaprofs[0,*,*,0]), dim = 2), /nodat, $
        xtitle = '!18r/r!D!Xe!N', $
;        ytitle = greek('Sigma')+'!D!18M!X!L*!N [M!D'+sunsymbol()+'!N kpc!E-2!N]', $
        ytickname = replicate(' ', 60), $
        xran = !X.CRANGE+[1d-6,0], /xsty, $
        yran = [top-2.5,top], ysty = 8+1, $
        pos = [0.5,0.15,0.9,0.9], /noer, $
        xthick = 3, ythick = 3, charsize = 1.25, charthick = 4
  axis, yaxis = 1, ytitle = 'log '+greek('Sigma')+'!D!18M!X!L*!N!18/<!X'+greek('Sigma')+'!D!18M!X!L*!N(!18r!X<1 kpc)!18>!X', $;, $
        ythick = 3, charsize = 1.25, charthick = 4, $
        yran = !Y.CRANGE, /ysty
  for ii = 0, ngals - 1 do begin
;     tre = median(hmrs[ii,*,0], dim = 2)
;     norm = value_locate(rad, tre)
;     norm = norm[0]
;     norm = 0
     mp = reform(msigmaprofs[ii,*,*,0])     
     qui = where(prs[ii,*] le 1.) ;; Inner Kpc...
     norm = mean(mp[qui,*])
     tr = resreg[sort(extrads[ii,resreg])]
     nuse = n_elements(tr)
     ttr = reform(extrads[ii,tr])
     tmasses = reform(masses[ii,tr,0])
     etmasses = reform(emasses[ii,tr,0])
     oplot, ttr, tmasses - areas[ii,tr,0] - norm, $
            thick = 10
     oplot, ttr, tmasses - areas[ii,tr,0] - norm, $
            thick = 4, col = cgcolor(string(cols[ii]))
     oploterror, ttr, tmasses - areas[ii,tr,0] - norm, $
                 replicate(emags[ii],nuse), reform(sqrt(etmasses + emags[ii]^2)), $
                 psym = 8, symsize = 1.3, /nohat
     oplot, ttr, tmasses - areas[ii,tr,0] - norm, $
            col = cgcolor(string(cols[ii])), psym = 8
;     mmm = median(mp, dim = 2)
;     eee = stddev(mp, dim = 2)
;     oplot, rad, mmm - eee - norm, thick = 4
;     oplot, rad, mmm + eee - norm, thick = 4
;     polyfill, xxx, [mmm-eee, reverse(mmm+eee)] - norm > !Y.CRANGE[0], $
;               col = cgcolor(string(cols[ii])), /line_fill, $
;               thick = 1, orien = ii * 45, spacing = 0.025
  endfor
  rep = prs[*,value_locate(rad, 1.0)]
  kcols = []
  for ii = 0, ngals- 1 do $
     kcols = [kcols, cgcolor(string(cols[ii]))]
  plotsym, 0, /fill
  key = []
  for ii = 0, ngals - 1 do $
     key = [key, string(rep[ii], f = '(F3.1)')+' kpc,'+string(ssfrs[ii], f = '(F5.1)')+' yr!E-1!N']
  cgtext, !X.crange[0] + 0.1, !Y.CRANGE[0]+0.8, /data, $
          '(!18r!X!De!N, log !18sSFR!X)=', charsize = 1, charthick = 3
  legend, /bottom, /left, box = 0, $
          key, $
          psym = replicate(8,ngals), $
          col = kcols, pspacing = 0.5, $
          charsize = 1, charthick = 3
  
  device, /close
;  spawn, 'gv massProfiles.eps &'

  
  ;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;
  ;; SFR Density ;;
  ;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;
  
  set_plot, 'PS'
  xsize = 8.5
  ysize = 4
  ar = ysize / xsize
  device, filename = 'sfrProfiles.eps', $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = xsize, ysize = ysize, /in
  plot, rad, reform(sfrs[*,resreg,0]) - areas[*,resreg], /nodat, $
        xtitle = '!18r/r!D!Xe!N', xran = [-0.1, 1.8], /xsty, $
        yran = [-2.5,0], $
        ytitle = 'log '+greek('Sigma')+'!D!18SFR!X!N [M!D'+sunsymbol()+'!N yr !E-1!N kpc!E-2!N]', $
        xthick = 3, ythick = 3, charsize = 1.25, charthick = 4, $
        pos = [0.1,0.15,0.5,0.9]
  plotsym, 0, /fill
  for ii = 0, ngals - 1 do begin
     tre = (median(hmrs[ii,*,0], dim = 2))[0]
     sre = (stddev(hmrs[ii,*,0], dim = 2, /nan))[0]
     oplot, replicate(tre,2), !Y.CRANGE[0] + [0,0.2], thick = 10
     oplot, tre + sre * [-1,1], $
            !Y.CRANGE[0] + [0.2,0.2], thick = 10
     oplot, replicate(tre,2), !Y.CRANGE[0] + [0,0.2], col = cgcolor(string(cols[ii])), thick = 6
     oplot, tre + sre * [-1,1], $
            !Y.CRANGE[0] + [0.2,0.2], thick = 6, col = cgcolor(string(cols[ii]))

     ;; do detections and upper limits...
     tr = resreg[sort(extrads[ii,resreg])]
     use = where(sfrs[ii,tr,0] ne 0, nuse)
     uls = where(sfrs[ii,tr,0] eq 0, nuls)
     oplot, extrads[ii,tr], alog10(sfrs[ii,tr,0] > nsig * lsfrs[ii,tr,0])- areas[ii,tr,0], $
            thick = 10, linesty = 2
     oplot, extrads[ii,tr], alog10(sfrs[ii,tr,0] > nsig * lsfrs[ii,tr,0]) - areas[ii,tr,0], $
            thick = 4, col = cgcolor(string(cols[ii])), linesty = 2
     if nuls gt 0 then begin
        plotsym, 1, thick = 10
        oplot, [extrads[ii,tr[uls],0]], [alog10(nsig * lsfrs[ii,tr[uls],0]) - areas[ii,tr[uls],0]], $
               psym = 8, symsize = 2
        plotsym, 1, thick = 6
        oplot, [extrads[ii,tr[uls],0]], [alog10(nsig * lsfrs[ii,tr[uls],0]) - areas[ii,tr[uls],0]], $
               psym = 8, col = cgcolor(string(cols[ii])), symsize = 2
     endif
     if nuse gt 0 then begin
        plotsym, 0, /fill
        luse = [use[0]-1, use] > 0       
        tesfrs = 1./alog(10) * esfrs[ii,tr[use],0]/sfrs[ii,tr[use],0]
        oplot, [extrads[ii,tr[luse]]], reform([alog10(sfrs[ii,tr[luse],0]) - areas[ii,tr[luse],0]]), $
               thick = 10
        oplot, [extrads[ii,tr[luse]]], reform([alog10(sfrs[ii,tr[luse],0]) - areas[ii,tr[luse],0]]), $
               thick = 4, col = cgcolor(string(cols[ii]))
        oploterror, reform(extrads[ii,tr[use]]), reform(alog10(sfrs[ii,tr[use],0]) - areas[ii,tr[use],0]), $
                    [replicate(emags[ii],nuse)], reform(sqrt(tesfrs^2 + emags[ii]^2)), $
                    psym = 8, symsize = 1.3, /nohat
        oplot, [extrads[ii,tr[use]]], [alog10(sfrs[ii,tr[use],0]) - areas[ii,tr[use],0]], $
               col = cgcolor(string(cols[ii])), psym = 8
     endif
     
  endfor
;  cgtext, 1.1, !Y.CRANGE[1]-0.1, 'half-mass radii', $
;          orien = 270, align = 0, charsize = 1, charthick = 3
  
  plot, rad, sfrs[0,*,0], /nodat, $
        xtitle = '!18r/r!D!Xe!N', $
;        ytitle = greek('Sigma')+'!D!18M!X!L*!N [M!D'+sunsymbol()+'!N kpc!E-2!N]', $
        ytickname = replicate(' ', 60), $
        xran = !X.CRANGE+[1d-6,0], /xsty, $
        yran = [-11,-8.3], ysty = 8+1, $
        pos = [0.5,0.15,0.9,0.9], /noer, $
        xthick = 3, ythick = 3, charsize = 1.25, charthick = 4
  axis, yaxis = 1, ytitle = 'log !18SFR!N/M!X!D*!N [yr!E-1!N]', $;, $
        ythick = 3, charsize = 1.25, charthick = 4, $
        yran = !Y.CRANGE, /ysty
  for ii = 0, ngals - 1 do begin
     ;; do detections and upper limits...
     tr = resreg[sort(extrads[ii,resreg])]
     use = where(sfrs[ii,tr,0] ne 0, nuse)
     uls = where(sfrs[ii,tr,0] eq 0, nuls)
     oplot, extrads[ii,tr], alog10(sfrs[ii,tr,0] > nsig * lsfrs[ii,tr,0])- masses[ii,tr,0], $
            thick = 10, linesty = 2
     oplot, extrads[ii,tr], alog10(sfrs[ii,tr,0] > nsig * lsfrs[ii,tr,0]) - masses[ii,tr,0], $
            thick = 4, col = cgcolor(string(cols[ii])), linesty = 2
     if nuls gt 0 then begin
        plotsym, 1, thick = 10
        oplot, [extrads[ii,tr[uls],0]], [alog10(nsig * lsfrs[ii,tr[uls],0]) - masses[ii,tr[uls],0]], $
               psym = 8, symsize = 2
        plotsym, 1, thick = 6
        oplot, [extrads[ii,tr[uls],0]], [alog10(nsig * lsfrs[ii,tr[uls],0]) - masses[ii,tr[uls],0]], $
               psym = 8, col = cgcolor(string(cols[ii])), symsize = 2
     endif
     if nuse gt 0 then begin
        plotsym, 0, /fill
        luse = [use[0]-1, use] > 0       
        tesfrs = 1./alog(10) * esfrs[ii,tr[use],0]/sfrs[ii,tr[use],0]
        oplot, [extrads[ii,tr[luse]]], reform([alog10(sfrs[ii,tr[luse],0]) - masses[ii,tr[luse],0]]), $
               thick = 10
        oplot, [extrads[ii,tr[luse]]], reform([alog10(sfrs[ii,tr[luse],0]) - masses[ii,tr[luse],0]]), $
               thick = 4, col = cgcolor(string(cols[ii]))
        oploterror, reform(extrads[ii,tr[use]]), reform(alog10(sfrs[ii,tr[use],0]) - masses[ii,tr[use],0]), $
                    [emasses[ii,tr[use],0]], reform(sqrt(tesfrs^2 + emasses[ii,tr[use],0]^2)), $
                    psym = 8, symsize = 1.3, /nohat
        oplot, [extrads[ii,tr[use]]], [alog10(sfrs[ii,tr[use],0]) - masses[ii,tr[use],0]], $
               col = cgcolor(string(cols[ii])), psym = 8
     endif
  endfor
  rep = prs[*,value_locate(rad, 1.0)]
  kcols = []
  for ii = 0, ngals- 1 do $
     kcols = [kcols, cgcolor(string(cols[ii]))]
  plotsym, 0, /fill
  key = []
  for ii = 0, ngals - 1 do $
     key = [key, string(rep[ii], f = '(F3.1)')+' kpc,'+string(ssfrs[ii], f = '(F5.1)')+' yr!E-1!N']
;  cgtext, !X.crange[-1] - 0.1, !Y.CRANGE[0]+0.9, /data, $
;          '(!18r!X!De!N, log !18sSFR!X)=', charsize = 1, charthick = 3, align = 1
  legend, /bottom, /right, box = 0, $
          key, $
          psym = replicate(8,ngals), $
          col = kcols, pspacing = 0.5, $
          charsize = 1, charthick = 3
  
  device, /close
;  spawn, 'gv sfrProfiles.eps &'
;  stop

  ;;;;;;;;;;
  ;;;;;;;;;;
  ;; DUST ;;
  ;;;;;;;;;;
  ;;;;;;;;;;
  
  set_plot, 'PS'
  xsize = 8.5
  ysize = 4.
  ar = ysize / xsize
  device, filename = 'dustProfiles.eps', $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = xsize, ysize = ysize, /in
  plot, rad, median(reform(dustprofs[0,*,*,0]), dim = 2), /nodat, $
        xran = [-0.1, 1.8], /xsty, $
        xtitle = '!18r/r!D!Xe!N', yran = [0,2], $
        ytitle = 'A!D!18V!X,LogNormal!N [Mag]', $
        xthick = 3, ythick = 3, charsize = 1.25, charthick = 4, $
        pos = [0.1,0.15,0.5,0.9]
  plotsym, 0, /fill
  for ii = 0, ngals - 1 do begin

     tavs = reform(avs[ii,tr,0])
     etavs = reform(eavs[ii,tr,0])     
     ttr = reform(extrads[ii,tr])
     nuse = n_elements(tr)
     
     oplot, ttr, tavs, $
            thick = 10
     oplot, ttr, tavs, $
            thick = 4, col = cgcolor(string(cols[ii]))
     oploterror, ttr, tavs, $
                 replicate(emags[ii],nuse), etavs, $
                 psym = 8, symsize = 1.3, /nohat
     oplot, ttr, tavs, $
            col = cgcolor(string(cols[ii])), psym = 8

;     mp = reform(dustprofs[ii,*,*,0])     
;     mmm = median(mp, dim = 2)
;     eee = stddev(mp, dim = 2)
;     oplot, rad, mmm+eee, thick = 4
;     oplot, rad, mmm-eee, thick = 4     
;     polyfill, xxx < !X.CRANGE[1], [mmm-eee, reverse(mmm+eee)] < !Y.CRANGE[1], $
;               col = cgcolor(string(cols[ii])), /line_fill, $
;               thick = 1, orien = ii * 45, spacing = 0.025
  endfor
;  for ii = 0, ngals - 1 do begin
;     oplot, extrads[ii,resreg], avs[ii,resreg,1], col = '777777'x, thick = 4
;     oplot, extrads[ii,resreg], avs[ii,resreg,1], col = cgcolor(string(cols[ii])), thick = 2
;     oploterror, extrads[ii,resreg], avs[ii,resreg,1], eavs[ii,resreg,1], $
;                 psym = 8, symsize = 1.25, col = '777777'x, errcol = '777777'x, $
;                 errthick = 4
;     oplot, extrads[ii,resreg], avs[ii,resreg,1], psym = 8, col = cgcolor(string(cols[ii]))
;;     oplot, rad, mmm+eee, thick = 6, linesty = 2
;;     oplot, rad, mmm-eee, thick = 6, linesty = 2
;;     oplot, rad, mmm+eee, thick = 2, col = cgcolor(string(cols[ii])), linesty = 2
;;     oplot, rad, mmm-eee, thick = 2, col = cgcolor(string(cols[ii])), linesty = 2
;  endfor
     
;  cgtext, 1.1, !Y.CRANGE[1]-0.1, 'half-mass radii', $
;          orien = 270, align = 0, charsize = 1, charthick = 3
  
  plot, rad, sfrs[0,*,0], /nodat, $
        xtitle = '!18r/r!D!Xe!N', $
;        ytitle = greek('Sigma')+'!D!18M!X!L*!N [M!D'+sunsymbol()+'!N kpc!E-2!N]', $
        ytickname = replicate(' ', 60), $
        xran = !X.CRANGE+[1d-6,0], /xsty, $
        yran = !Y.CRANGE, ysty = 8+1, $
        pos = [0.5,0.15,0.9,0.9], /noer, $
        xthick = 3, ythick = 3, charsize = 1.25, charthick = 4
  axis, yaxis = 1, ytitle = 'A!D!18V!X,DelayedExponential!N [mag]', $;, $
        ythick = 3, charsize = 1.25, charthick = 4, $
        yran = !Y.CRANGE, /ysty
  for ii = 0, ngals - 1 do begin
       
     tavs = reform(avs[ii,tr,1])
     etavs = reform(eavs[ii,tr,1])
     ttr = reform(extrads[ii,tr])
     nuse = n_elements(tr)
     
     oplot, ttr, tavs, $
            thick = 10, col = '777777'x
     oplot, ttr, tavs, $
            thick = 4, col = cgcolor(string(cols[ii]))
     oploterror, ttr, tavs, $
                 replicate(emags[ii],nuse), etavs, $
                 psym = 8, symsize = 1.3, /nohat, col = '777777'x, errcol = '777777'x
     oplot, ttr, tavs, $
            col = cgcolor(string(cols[ii])), psym = 8

  endfor
  rep = prs[*,value_locate(rad, 1.0)]
  kcols = []
  for ii = 0, ngals- 1 do $
     kcols = [kcols, cgcolor(string(cols[ii]))]
  plotsym, 0, /fill
  key = []
  for ii = 0, ngals - 1 do $
     key = [key, string(rep[ii], f = '(F3.1)')+' kpc,'+string(ssfrs[ii], f = '(F5.1)')+' yr!E-1!N']
;  cgtext, !X.crange[0] + 0.1, !Y.CRANGE[1]-0.1, /data, $
;          '(!18r!X!De!N, log !18sSFR!X)=', charsize = 1, charthick = 3
  legend, /bottom, /right, box = 0, $
          key, $
          psym = replicate(8,ngals), $
          col = kcols, pspacing = 0.5, $
          charsize = 1, charthick = 3
  
  device, /close
;  spawn, 'gv dustProfiles.eps &'

  ;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;
  ;; SFH stuff ;;
  ;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;

  set_plot, 'PS'
  xsize = 8.5
  ysize = 4.
  ar = ysize / xsize
  device, filename = 'expProfiles.eps', $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = xsize, ysize = ysize, /in
  plot, rad, median(reform(dustprofs[0,*,*,0]), dim = 2), /nodat, $
        xran = [-0.1, 1.8], /xsty, $
        xtitle = '!18r/r!D!Xe!N', yran = [5,0], $
        ytitle = 'age [Gyr]', $
        xthick = 3, ythick = 3, charsize = 1.25, charthick = 4, $
        pos = [0.1,0.15,0.5,0.9], yminor = 2
  plotsym, 0, /fill
  for ii = 0, ngals - 1 do begin

     tavs = reform(ages[ii,tr])
     etavs = reform(eages[ii,tr])     
     ttr = reform(extrads[ii,tr])
     nuse = n_elements(tr)
     
     oplot, ttr, tavs, $
            thick = 10, col = '777777'x
     oplot, ttr, tavs, $
            thick = 4, col = cgcolor(string(cols[ii]))
     oploterror, ttr, tavs, $
                 replicate(emags[ii],nuse), etavs, $
                 psym = 8, symsize = 1.3, /nohat, col = '777777'x, errcol = '777777'x
     oplot, ttr, tavs, $
            col = cgcolor(string(cols[ii])), psym = 8

  endfor
  rep = prs[*,value_locate(rad, 1.0)]
  kcols = []
  for ii = 0, ngals- 1 do $
     kcols = [kcols, cgcolor(string(cols[ii]))]
  plotsym, 0, /fill
  key = []
  for ii = 0, ngals - 1 do $
     key = [key, string(rep[ii], f = '(F3.1)')+' kpc,'+string(ssfrs[ii], f = '(F5.1)')+' yr!E-1!N']
;  cgtext, !X.crange[0] + 0.1, !Y.CRANGE[1]-0.1, /data, $
;          '(!18r!X!De!N, log !18sSFR!X)=', charsize = 1, charthick = 3
  legend, /bottom, /right, box = 0, $
          key, $
          psym = replicate(8,ngals), $
          col = kcols, pspacing = 0.5, $
          charsize = 1, charthick = 3
  
  plot, rad, sfrs[0,*,0], /nodat, $
        xtitle = '!18r/r!D!Xe!N', $
;        ytitle = greek('Sigma')+'!D!18M!X!L*!N [M!D'+sunsymbol()+'!N kpc!E-2!N]', $
        ytickname = replicate(' ', 60), $
        xran = !X.CRANGE+[1d-6,0], /xsty, $
        yran = [0,10], ysty = 8+1, $
        pos = [0.5,0.15,0.9,0.9], /noer, $
        xthick = 3, ythick = 3, charsize = 1.25, charthick = 4, yminor = 2
  axis, yaxis = 1, ytitle = 'SFH !18e!X-foldings', $;, $
        ythick = 3, charsize = 1.25, charthick = 4, $
        yran = !Y.CRANGE, /ysty, yminor = 2
  for ii = 0, ngals - 1 do begin
       
     tavs = reform(efolds[ii,tr])
     if total(tavs lt !Y.CRANGE[1]) eq nuse then begin
        oplot, ttr, tavs, $
               thick = 10, col = '777777'x
        oplot, ttr, tavs, $
               thick = 4, col = cgcolor(string(cols[ii]))
        oplot, ttr, tavs, $
               psym = 8, symsize = 1.3, col = '777777'x
        oplot, ttr, tavs, $
               col = cgcolor(string(cols[ii])), psym = 8
     endif else begin
        oplot, ttr, replicate(!Y.CRANGE[1] - 1, nuse), $
               thick = 10, col = '777777'x
        oplot, ttr, replicate(!Y.CRANGE[1] - 1, nuse), $
               thick = 4, col = cgcolor(string(cols[ii]))
        plotsym, 2, thick = 10
        oplot, ttr, replicate(!Y.CRANGE[1] - 1, nuse), psym = 8, $
               thick = 10, col = '777777'x, symsize = 2
        plotsym, 2, thick = 6
        oplot, ttr, replicate(!Y.CRANGE[1] - 1, nuse), psym = 8, $
               thick = 6, col = cgcolor(string(cols[ii])), symsize = 2
        plotsym, 0, /fill
     endelse

  endfor
  
  device, /close
  set_plot, 'X'
  spawn, 'gv expProfiles.eps &'
  
  plotsym, 0, /fill
  xw = (0.9 - 0.1) / 3.
  set_plot, 'PS'
  xsize = 8.5
  ysize = 3.
  ar = ysize / xsize
  device, filename = 'lgnProfiles.eps', $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = xsize, ysize = ysize, /in
  !p.multi = [0,3,0]
  plot, ttr, t0s[0,*], /nodat, $
        xran = [-0.1, 1.8], /xsty, $
        xtitle = '!18r/r!D!Xe!N', yran = [0.5,2], $
        ytitle = '!18T!D!X0!N [ln(Gyr)]', $
        xthick = 3, ythick = 3, charsize = 2, charthick = 4;, $
;        pos = [0.1,0.15,0.1+xw,0.9]
  plotsym, 0, /fill
  for ii = 0, ngals - 1 do begin

     tavs = reform(t0s[ii,tr])
     etavs = reform(et0s[ii,tr])     
     ttr = reform(extrads[ii,tr])
     nuse = n_elements(tr)
     
     oplot, ttr, tavs, $
            thick = 10
     oplot, ttr, tavs, $
            thick = 4, col = cgcolor(string(cols[ii]))
     oploterror, ttr, tavs, $
                 replicate(emags[ii],nuse), etavs, $
                 psym = 8, symsize = 1.3, /nohat
     oplot, ttr, tavs, $
            col = cgcolor(string(cols[ii])), psym = 8

  endfor

  plot, ttr, taus[0,*], /nodat, $
        xtitle = '!18r/r!D!Xe!N', $
        ytitle =  greek('tau')+' [ln(Gyr)]', $
;        ytickname = replicate(' ', 60), $
        xran = !X.CRANGE+[1d-6,0], /xsty, $
        yran = [0,2.5], ysty = 1, $
;        pos = [!X.WINDOW[1],0.15,!X.WINDOW[1]+xw,0.9], $; /noer, $
        xthick = 3, ythick = 3, charsize = 2, charthick = 4
;  axis, yaxis = 1, ytitle = greek('tau')+' [ln(Gyr)]', $;, $
;        ythick = 3, charsize = 2, charthick = 4, $
;        yran = !Y.CRANGE, /ysty
  for ii = 0, ngals - 1 do begin
       
     tavs = reform(taus[ii,tr])
     etavs = reform(etaus[ii,tr])
     ttr = reform(extrads[ii,tr])
     nuse = n_elements(tr)
     
     oplot, ttr, tavs, $
            thick = 10
     oplot, ttr, tavs, $
            thick = 4, col = cgcolor(string(cols[ii]))
     oploterror, ttr, tavs, $
                 replicate(emags[ii],nuse), etavs, $
                 psym = 8, symsize = 1.3, /nohat
     oplot, ttr, tavs, $
            col = cgcolor(string(cols[ii])), psym = 8

  endfor

  plot, ttr, pktimes[0,*], /nodat, $
        xtitle = '!18r/r!D!Xe!N', $
        ytitle = '!18T!X!Dpeak!N [ln(Gyr)]', $
;        ytickname = replicate(' ', 60), $
        xran = !X.CRANGE+[1d-6,0], xsty = 1, $
        yran = [-1,alog(10)], ysty = 8+1, $
;        pos = [!X.WINDOW[1],0.15,!X.WINDOW[1]+xw,0.9], /noer, $
        xthick = 3, ythick = 3, charsize = 2, charthick = 4
  ty = findgen((exp(!Y.CRANGE[1]) - exp(!Y.CRANGE[0]))/0.1+1)*0.1 + exp(!Y.CRANGE[0])
  axis, yaxis = 1, ytitle = '!18z!X!Dpeak!N', $
        ythick = 3, charsize = 2, charthick = 4, $
;        yran = reverse(getredshift(minmax(ty))), /ysty, $
        ytickv = alog(ty[value_locate(getredshift(ty), [8,6,4,2,1])]), $
        ytickname = ['8','6','4','2','1'], yticks = 4, yminor = 2
  for ii = 0, ngals - 1 do begin
       
     tavs = reform(pktimes[ii,tr])
     etavs = reform(epktimes[ii,tr])
     ttr = reform(extrads[ii,tr])
     nuse = n_elements(tr)
     
     oplot, ttr, tavs, $
            thick = 10
     oplot, ttr, tavs, $
            thick = 4, col = cgcolor(string(cols[ii]))
     oploterror, ttr, tavs, $
                 replicate(emags[ii],nuse), etavs, $
                 psym = 8, symsize = 1.3, /nohat
     oplot, ttr, tavs, $
            col = cgcolor(string(cols[ii])), psym = 8

  endfor
  
;  rep = prs[*,value_locate(rad, 1.0)]
;  kcols = []
;  for ii = 0, ngals- 1 do $
;     kcols = [kcols, cgcolor(string(cols[ii]))]
;  plotsym, 8, /fill
;  key = []
;  for ii = 0, ngals - 1 do $
;     key = [key, string(rep[ii], f = '(F3.1)')+' kpc,'+string(ssfrs[ii], f = '(F5.1)')+' yr!E-1!N']
;;  cgtext, !X.crange[0] + 0.1, !Y.CRANGE[1]-0.1, /data, $
;;          '(!18r!X!De!N, log !18sSFR!X)=', charsize = 1, charthick = 3
;  legend, /bottom, /right, box = 0, $
;          key, $
;          psym = replicate(8,ngals), $
;          col = kcols, pspacing = 0.5, $
;          charsize = 1, charthick = 3
  
  device, /close
  spawn, 'gv lgnProfiles.eps &'
  set_plot, 'X'
  !P.multi = 0
  
  stop

;  plot, taus[*,master], bts[*,0], psym = 1, xran = [0,3], yran = [0,1]     
;  oplot, tdat.TAU, tdat.MBT, psym = 1, col = '777777'x                
;  oploterror, taus[*,master], bts[*,0], etaus[*,master], ebts[*,0]

;  case quantity of
;     'MASS'
;     'MDENSITY'
;     'SDENSITY'
;     'SSFR'
;     'A_V'
;     'AGE'
;     'T0'
;     'TAU'
;     'PEAKTIME'
;  endcase
  
  
  
end
;plotprofiles, 'study5_lgn_Bestfits.list', 'study5_exp_Bestfits.list'
