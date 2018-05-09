
;; Plot the observed profiles of the galaxies in any quantity
pro plotPaperProfileFigure, lgnlist, explist, output, $
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
  
  extract = where($;regions eq 'OPTFL' OR $
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
;  extareas = fltarr(ngals, nregions) ;; LEA 20170731
    
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
;  aprofs = fltarr(ngals, n_elements(rad), niter, 2) ;; LEA 20170731 -- need area as fn of r
  dustprofs = fltarr(ngals, n_elements(rad), niter, 2)
  hmrs = fltarr(ngals, niter, 2)
  prs  = fltarr(ngals, n_elements(rad))

  nsig = 2.
  names = []
  for ii = 0, ngals - 1 do begin
     for jj = 0, 1 do begin
        case jj of
           0: begin
              file = lgnfiles[ii]
              files[ii] = strmid(file, 0, 7)
              case files[ii] of
                 '00900_1': nn = 'SSF'
                 '00451_2': nn = 'CSF'
                 '01916_2': nn = 'PSB'
                 '00660_2': nn = 'PAS'
              endcase
              names = [names, nn]
              qui = where(summaryData.ID eq files[ii])
              mag  = summaryData[qui].MU
              emag = 1./alog(10) * summaryData[qui].EMU / mag
              mags[ii] = mag
              emags[ii] = emag
           end
           1: file = expfiles[ii]
        endcase

        ;; Read the summary fitting results
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
;           extareas[ii,*] = data.AREA_PHYS
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

        prs[ii,*] = rp * rad / sqrt(mag)
        
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
        m    = reform(masses[ii,extract,jj] - areas[ii,extract,jj])
        em   = reform(sqrt(emasses[ii,extract,jj]^2 + eareas[ii,extract,jj]^2))
        s    = reform(alog10(sfrs[ii,extract,jj]) - areas[ii,extract,jj])
        es   = reform(sqrt((1./alog(10)*esfrs[ii,extract,jj]/sfrs[ii,extract,jj])^2 + eareas[ii,extract,jj]^2))
        ssfr = s - m
        essfr = reform(sqrt((1./alog(10)*esfrs[ii,extract,jj]/sfrs[ii,extract,jj])^2 + emasses[ii,extract,jj]^2))
        av   = data[extract].A_V[0]
        eav  = 0.5 * (data[extract].A_V[2] - data[extract].A_V[1])
       
        tr   = data[extract].RNORM
;        a    = interpol(areas[ii,extract,jj], tr, rad)
        for ll = 0, niter - 1 do begin
           tm = m + randomn(seed, 3) * em
           ts = s + randomn(seed, 3) * es
           tssfr = ssfr + randomn(seed, 3) * essfr
           tav = av + randomn(seed, 3) * eav

           ttr = tr + randomn(seed, 3) * (10.^(0.5 * emags[ii]) - 1) * tr
           
           msigmaprofs[ii,*,ll,jj] = interpol(tm, ttr, rad)
           dustprofs[ii,*,ll,jj] = interpol(tav, ttr, rad)     

           mprofs[ii,*,ll,jj] = total(2 * !pi * rad * dr * rp^2 * 10.^msigmaprofs[ii,*,ll,jj], /cum) ;; mag independent!;total(10.^msigmaprofs[ii,*,ll,jj] * a, /cum);
           mprofs[ii,*,ll,jj] /= mprofs[ii,-1,ll,jj] / 10.^(masses[ii,master,jj] + randomn(seed, 1) * emasses[ii,master,jj])[0]
           
           hmrs[ii,ll,jj] = rad[value_locate(mprofs[ii,*,ll,jj], 0.5 * mprofs[ii,-1,ll,jj])]           
        endfor

        savedata = {R: rad, RP: rp / sqrt(mag), $
                    M_R: median(mprofs[ii,*,*,jj], dim = 3), $
                    EM_R: stddev(mprofs[ii,*,*,jj], dim = 3, /nan)}
        if jj eq 0 then $
           mwrfits, savedata, files[ii]+'_lgn_circprof.fits', /create $
        else $
           mwrfits, savedata, files[ii]+'_exp_circprof.fits', /create
                    
     endfor
     
  endfor

;  stop
  
;  tsfrs = mean(sfrs[*,master,*], dim = 3, /nan)
  tsfrs = sfrs[*,master,0];, dim = 3, /nan)
  tsfrs[where(tsfrs eq 0)] = nsig*esfrs[where(tsfrs eq 0),master,0]
  resreg = extract;[value_locate(regions[extract], ['INNFL', 'INTER', 'OUTER'])]
  
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

;  stop
  
  ;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;
  ;; Mass Density ;;
  ;;  SFR Density ;;
  ;;    Av(R)     ;;
  ;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;
  
  set_plot, 'PS'
  xsize = 5.5
  ysize = 8.5
  ar = ysize / xsize

  xm = 0.2
  xw = (1. - 2*xm) / 2.
  x0 = xm
  x1 = xm + xw

  ym = 0.1
  yw = (1. - (ym+0.05)) / 3. ;; top margin can be smaller
  y0 = ym + 2 * yw
  y1 = ym + yw
  y2 = ym
  
  device, filename = 'multiProfiles.eps', $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = xsize, ysize = ysize, /in
  plot, rad, median(reform(msigmaprofs[0,*,*,0]), dim = 2), /nodat, $
;        xtitle = '!18r/r!D!Xe!N', $
        xtickname = replicate(' ', 60), $
        xran = [-0.1,1.8], $
        yran = [7.5,10], $
        ytitle = 'log '+greek('Sigma')+'!D!18M!X!L*!N [M!D'+sunsymbol()+'!N kpc!E-2!N]', $
        xthick = 3, ythick = 3, charsize = 1.3, charthick = 4, $
        pos = [x0,y0,x1,y0+yw], /xsty
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
     tmasses = reform(masses[ii,extract,0])
     etmasses = reform(emasses[ii,extract,0])
     ttr = reform(extrads[ii,extract])
     nuse = n_elements(extract)     

     oplot, ttr, tmasses - reform(areas[ii,extract,0]), $
            thick = 10
     oplot, ttr, tmasses - reform(areas[ii,extract,0]), $
            thick = 4, col = cgcolor(string(cols[ii]))
     oploterror, ttr, tmasses - reform(areas[ii,extract,0]), $
                 (10.^(0.5*emags[ii])-1) * ttr, reform(sqrt(etmasses^2 + emags[ii]^2)), $
                 psym = 8, symsize = 1.3, /nohat
     oplot, ttr, tmasses - reform(areas[ii,extract,0]), $
            col = cgcolor(string(cols[ii])), psym = 8

  endfor
  
  top = 0.3
  plot, rad, median(reform(msigmaprofs[0,*,*,0]), dim = 2), /nodat, $
;        xtitle = '!18r/r!D!Xe!N', $
;        ytitle = greek('Sigma')+'!D!18M!X!L*!N [M!D'+sunsymbol()+'!N kpc!E-2!N]', $
        xtickname = replicate(' ', 60), $
        ytickname = replicate(' ', 60), $
        xran = !X.CRANGE+[1d-6,0], /xsty, $
        yran = [top-2,top], ysty = 8+1, $
        pos = [x1,y0,x1+xw,y0+yw], /noer, $
        xthick = 3, ythick = 3, charsize = 1.3, charthick = 4
  axis, yaxis = 1, ytitle = 'log '+greek('Sigma')+'!D!18M!X!L*!N!18/<!X'+greek('Sigma')+'!D!18M!X!L*!N(!18r!X<1 kpc)!18>!X', $
        ythick = 3, charsize = 1.3, charthick = 4, $
        yran = !Y.CRANGE, /ysty
  for ii = 0, ngals - 1 do begin
     mp = reform(msigmaprofs[ii,*,*,0])     
     qui = where(reform(prs[ii,*]) le 1.) ;; Inner Kpc...
     norm = mean(mp[qui,*])
     tr = extract
     nuse = n_elements(tr)
     ttr = reform(extrads[ii,tr])
     tmasses = reform(masses[ii,tr,0])
     etmasses = reform(emasses[ii,tr,0])
     oplot, ttr, tmasses - reform(areas[ii,tr,0]) - norm, $
            thick = 10
     oplot, ttr, tmasses - reform(areas[ii,tr,0]) - norm, $
            thick = 4, col = cgcolor(string(cols[ii]))
     oploterror, ttr, tmasses - areas[ii,tr,0] - norm, $
                 (10.^(0.5*emags[ii])-1) * ttr, reform(sqrt(etmasses^2 + emags[ii]^2)), $
                 psym = 8, symsize = 1.3, /nohat
     oplot, ttr, tmasses - areas[ii,tr,0] - norm, $
            col = cgcolor(string(cols[ii])), psym = 8
  endfor
  rep = prs[*,value_locate(rad, 1.0)]
  kcols = []
  for ii = 0, ngals- 1 do $
     kcols = [kcols, cgcolor(string(cols[ii]))]
  plotsym, 0, /fill
  key = []
  for ii = 0, ngals - 1 do $
     key = [key, string(rep[ii], f = '(F3.1)')+','+string(ssfrs[ii], f = '(F5.1)')+'; '+names[ii]]
  cgtext, mean(!X.crange), -0.925, /data, align = 0.5, $
          '[!18r!X!De!N, log!18sSFR!X; obj.]', charsize = 1, charthick = 3
  legend, pos = [!X.CRANGE[0],-1], /data, /left, box = 0, $
          key, $
          psym = replicate(8,ngals), $
          col = kcols, pspacing = 0.5, $
          charsize = 1, charthick = 3

  ;;
  ;;
  ;;
  
  plot, rad, reform(sfrs[*,resreg,0]) - areas[*,resreg], /nodat, $
;        xtitle = '!18r/r!D!Xe!N', $
        xtickname = replicate(' ', 60), $
        xran = [-0.1, 1.8], /xsty, $
        yran = [-2.5,0.1], /ysty, $
        ytitle = 'log '+greek('Sigma')+'!D!18SFR!X!N [M!D'+sunsymbol()+'!N yr !E-1!N kpc!E-2!N]', $
        xthick = 3, ythick = 3, charsize = 1.3, charthick = 4, $
        pos = [x0,y1,x0+xw,y0], /noer
  plotsym, 0, /fill
  for ii = 0, ngals - 1 do begin

     ;; do detections and upper limits...
     tr = resreg[sort(extrads[ii,resreg])]
     use = where(sfrs[ii,tr,0] ne 0, nuse)
     uls = where(sfrs[ii,tr,0] eq 0, nuls)
     oplot, extrads[ii,tr], alog10(sfrs[ii,tr,0] > nsig * lsfrs[ii,tr,0])- areas[ii,tr,0], $
            thick = 10, linesty = 1
     oplot, extrads[ii,tr], alog10(sfrs[ii,tr,0] > nsig * lsfrs[ii,tr,0]) - areas[ii,tr,0], $
            thick = 4, col = cgcolor(string(cols[ii])), linesty = 1
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
                    (10.^(0.5*emags[ii])-1) * reform(extrads[ii,tr[use]]), reform(sqrt(tesfrs^2 + emags[ii]^2)), $
                    psym = 8, symsize = 1.3, /nohat
        oplot, [extrads[ii,tr[use]]], [alog10(sfrs[ii,tr[use],0]) - areas[ii,tr[use],0]], $
               col = cgcolor(string(cols[ii])), psym = 8
     endif
     
  endfor
  
  plot, rad, sfrs[0,*,0], /nodat, $
;        xtitle = '!18r/r!D!Xe!N', $
;        ytitle = greek('Sigma')+'!D!18M!X!L*!N [M!D'+sunsymbol()+'!N kpc!E-2!N]', $
        xtickname = replicate(' ', 60), $
        ytickname = replicate(' ', 60), $
        xran = !X.CRANGE+[1d-6,0], /xsty, $
        yran = [-11,-8.3], ysty = 8+1, $
        pos = [x1,y1,x1+xw,y0], /noer, $
        xthick = 3, ythick = 3, charsize = 1.3, charthick = 4
  axis, yaxis = 1, ytitle = 'log !18SFR!N/M!X!D*!N [yr!E-1!N]', $;, $
        ythick = 3, charsize = 1.3, charthick = 4, $
        yran = !Y.CRANGE, /ysty
  for ii = 0, ngals - 1 do begin
     ;; do detections and upper limits...
     tr = resreg[sort(extrads[ii,resreg])]
     use = where(sfrs[ii,tr,0] ne 0, nuse)
     uls = where(sfrs[ii,tr,0] eq 0, nuls)
     oplot, extrads[ii,tr], alog10(sfrs[ii,tr,0] > nsig * lsfrs[ii,tr,0])- masses[ii,tr,0], $
            thick = 10, linesty = 1
     oplot, extrads[ii,tr], alog10(sfrs[ii,tr,0] > nsig * lsfrs[ii,tr,0]) - masses[ii,tr,0], $
            thick = 4, col = cgcolor(string(cols[ii])), linesty = 1
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
                    (10.^(0.5*emags[ii])-1) * reform(extrads[ii,tr[use]]), $
                    reform(sqrt(tesfrs^2 + emasses[ii,tr[use],0]^2)), $
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
     key = [key, string(rep[ii], f = '(F3.1)')+' kpc,'+string(ssfrs[ii], f = '(F5.1)')+' yr!E-1!N; '+names[ii]]

  ;;
  ;;
  ;;

  plot, rad, median(reform(dustprofs[0,*,*,0]), dim = 2), /nodat, $
        xran = [-0.1, 1.8], /xsty, $
        xtitle = '!18r/r!D!Xe!N', yran = [0,2.1], ytickint = 0.5, /ysty, $
        ytitle = 'A!D!18V!X,LgNml!N [mag]', $
        xthick = 3, ythick = 3, charsize = 1.3, charthick = 4, $
        pos = [x0,y2,x0+xw,y1], /noer
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
                 (10.^(0.5 * emags[ii]) - 1) * ttr, etavs, $
                 psym = 8, symsize = 1.3, /nohat
     oplot, ttr, tavs, $
            col = cgcolor(string(cols[ii])), psym = 8
  endfor

  plot, rad, sfrs[0,*,0], /nodat, $
        xtitle = '!18r/r!D!Xe!N', $
;        ytitle = greek('Sigma')+'!D!18M!X!L*!N [M!D'+sunsymbol()+'!N kpc!E-2!N]', $
        ytickname = replicate(' ', 60), $
        xran = !X.CRANGE+[1d-6,0], /xsty, $
        yran = !Y.CRANGE, ysty = 8+1, $
        pos = [x1,y2,x1+xw,y1], /noer, $
        xthick = 3, ythick = 3, charsize = 1.3, charthick = 4
  axis, yaxis = 1, ytitle = 'A!D!18V!X,DelExp!N [mag]', $;, $
        ythick = 3, charsize = 1.3, charthick = 4, $
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
                 (10.^(0.5*emags[ii])-1) * ttr, etavs, $
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
     key = [key, string(rep[ii], f = '(F3.1)')+' kpc,'+string(ssfrs[ii], f = '(F5.1)')+' yr!E-1!N; '+names[ii]]
 
  device, /close
  spawn, 'gv multiProfiles.eps &'
  
end
;plotPaperProfileFigure, 'study5_lgn_Bestfits.list', 'study5_exp_Bestfits.list'
