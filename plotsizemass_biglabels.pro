pro plotsizemass_biglabels, lgnlist, explist, $
                            CIRCULARIZE = circularize

  if NOT keyword_set(NSIG) then nsig = 1
  if NOT keyword_set(CIRCULARIZE) then circularize = 0 else circularize = 1
    
  summaryData = mrdfits('objectSummaryData.fits', 1)
  finalMasses = mrdfits('finalMasses.fits', 1)
  
  readcol, explist, expfiles, f = 'A'
  readcol, lgnlist, lgnfiles, f = 'A'

  ngals = n_elements(lgnfiles)

  test     = mrdfits(lgnfiles[0], 1)

  regions  = test.REGION
  nregions = n_elements(regions)

  sizes   = fltarr(ngals, 2)
  esizes  = fltarr(ngals, 2)
  ssizes  = fltarr(ngals, 2)
  essizes = fltarr(ngals, 2)
  areas   = fltarr(ngals, nregions, 2)
  eareas  = fltarr(ngals, nregions, 2)
  masses  = fltarr(ngals, nregions, 2)
  emasses = fltarr(ngals, nregions, 2)
  imasses = fltarr(ngals, 2)
  eimasses = fltarr(ngals, 2)
  
  ;; Masses are consistent w/in the errors including magnification
  
  mags = fltarr(ngals)
  emags = fltarr(ngals)
  sfrs = fltarr(ngals,2)
  bts  = fltarr(ngals, 2)
  ebts = fltarr(ngals,2)
  files = strarr(ngals)

  sbts = fltarr(ngals, 5)
  
  savedata = {ID: string(0), $
              OMASS: fltarr(2), $
              EOMASS: fltarr(2), $
              INTMASS: fltarr(2), $
              EINTMASS: fltarr(2), $
              BT_O: fltarr(2), $
              EBT_O: fltarr(2), $
              BT_INT: fltarr(2), $
              EBT_INT: fltarr(2), $
              SBTS: fltarr(5)}
  savedata = replicate(savedata, ngals)

  names = []
  
  dr = 0.02
  rad = findgen(2./dr+1) * dr; - 2.
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
                
        data = mrdfits(file, 1)
        sfrs[ii,jj] = data[where(data.REGION eq 'OPTFL')].SFR[0]

        q = where(data.REGION ne 'OPTFL' AND $
                  data.REGION ne 'INTER' AND $
                  data.REGION ne 'OUTER')

        uq = where(data.REGION eq 'INNFL' OR $
                   data.REGION eq 'INTER' OR $
                   data.REGION eq 'OUTER')

        master = where(data.REGION eq 'OPTFL')
        
        ro = data[0].RE_OBS
        rp = data[0].RE_PHYS
;        print, ro, rp
        r  = (mrdfits(data[0].FITSFILE, 1, /silent)).RE
        rs = (mrdfits(data[0].FITSFILE, 1, /silent)).RE_SEX

        pscale = ro / r
        kscale = rp / r
        
        ;; Correct for magnification
        data.LMASS -= alog10(mag)
        data.AREA_PHYS /= mag
        
        masses[ii,*,jj]  = data.LMASS[0]
        emasses[ii,*,jj] = sqrt((0.5 * (data.LMASS[2] - data.LMASS[1]))^2 + emag^2)
        
        areas[ii,*,jj]  = alog10(data.AREA_PHYS)
        eareas[ii,*,jj] = emag
        
        sizes[ii,jj]  = alog10(rp / sqrt(mag))
        esizes[ii,jj] = 0.5 * emag

        ssizes[ii,jj]  = alog10(rs * kscale / sqrt(mag))
        essizes[ii,jj] = 0.5 * emag

        ;; Estimate B/T
        ;; Mass w/in r_e over mass w/in 2 * re
        tdat     = data[uq]
        mmaster  = masses[ii,master[0],jj]
        emmaster = emasses[ii,master[0],jj]

        imasses[ii,jj]  = alog10(total(10.^data[q].LMASS[0]))
        eimasses[ii,jj] = sqrt(total(emasses[ii,q,jj]^2))
        
        tr   = tdat.RNORM        
;        tr = [-tdat[2].RNORM, -tdat[1].RNORM, tdat[0].RNORM, tdat[1].RNORM, tdat[2].RNORM]
        
        niter = 100
        sprofs = fltarr(n_elements(rad), niter)
        mprofs = fltarr(n_elements(rad), niter)
        tbts   = fltarr(niter,2)

        if NOT circularize then begin
           m    = tdat.LMASS[0] + [0, alog10(2), alog10(2)]                                     ;; FOLDED -- must multiply inter and outer by 2  ; - areas[ii,uq,jj]
           em   = 0.5 * (tdat.LMASS[2] - tdat.LMASS[1]) + [0, alog10(sqrt(2)), alog10(sqrt(2))] / alog(10) ;)^2) + eareas[ii,uq,jj]^2)
        endif else begin
           m    = tdat.LMASS[0] - areas[ii,uq,jj]
           em   = sqrt((0.5 * (tdat.LMASS[2] - tdat.LMASS[1])^2) + eareas[ii,uq,jj]^2)
        endelse

        for ll = 0, niter - 1 do begin
           tm = m + randomn(seed, n_elements(em)) * em
           if NOT circularize then begin
              t = total(10.^tm, /cum)                                                     ; 
              mprofs[*,ll] = interpol(t, tr, rad)                                         ;
              tmp1 = mprofs[value_locate(rad,0),ll] / mprofs[value_locate(rad,tr[-1]),ll] ;t[0] / t[-1];mprofs[*,ll] / mprofs[-1,ll]                           ;10.^mmaster ;
              tmp2 = mprofs[value_locate(rad,0),ll] / 10.^mmaster                         ;t[0] / 10.^(mmaster + randomn(seed, 1) * emmaster);mprofs[*,ll] / 10.^mmaster
           endif else begin
              t = interpol(tm, tr, rad)
              sprofs[*,ll] = total(2 * !pi * rad * dr * rp^2 * 10.^t, /cum) ;; MAG INDEPENDENT
              mprofs[*,ll] = sprofs[*,ll] / sprofs[-1,ll] * 10.^(mmaster + randomn(seed, 1) * emmaster)[0]
              tmp1 = mprofs[value_locate(rad * rp, 2.5),ll] / mprofs[-1,ll]
              tmp2 = mprofs[value_locate(rad * rp, 2.5),ll] / 10.^(mmaster + randomn(seed, 1) * emmaster)
           endelse
           tbts[ll,*] = [tmp1, tmp2] ;[tmp1[value_locate(rad, tr[0])], tmp2[value_locate(rad, tr[0])]]
        endfor
        maxlike = median(mprofs, dim = 2)
        maxmass = median(sprofs[-1,*])
        emaxmass = stddev(sprofs[-1,*], /nan)
        mwrfits, savedata, 'intOptMassComp.fits', /create
        perr = 1./alog(10) * stddev(mprofs, dim = 2, /nan) / maxlike
        plot, rad, alog10(maxlike), /ynoz, $
              yran = [8,11.5], linesty = 2
        oplot, replicate(2.5/rp, 2), !Y.CRANGE, col = '00ff00'x
        oplot, rad, alog10(maxlike) + perr
        oplot, rad, alog10(maxlike) - perr
        qui = where(data.REGION eq 'OPTFL')
        oplot, !X.CRANGE, replicate(masses[ii,qui,jj] - emasses[ii,qui,jj], 2), col = 255
        oplot, !X.CRANGE, replicate(masses[ii,qui,jj] + emasses[ii,qui,jj], 2), col = 255
        print, alog10(maxmass), 1./alog(10) * emaxmass/maxmass
        print, mmaster, emmaster ;masses[ii,qui,jj], emasses[ii,qui,jj]
        legend, /top, /left, [files[ii], string(jj, f = '(I)'), string(median(tbts[*,0]))]
;        stop
        
        bts[ii,jj]  = median(tbts[*,0]); 10.^data[uq[0]].LMASS[0] / 10.^imasses[ii,jj];total(10.^data[uq].LMASS[0]);;10.^data[uq[0]].LMASS[0] / total(10.^data[uq].LMASS[0]);median(tbts[*,0]); + median(tbts[*,1])); 
        ebts[ii,jj] = stddev(tbts[*,0], /nan);reform(sqrt(total(emasses[ii,q[1:-1],jj]^2)) * alog(10) * bts[ii,jj]); ;sqrt((0.5 * abs(median(tbts[*,0]) - median(tbts[*,1])))^2 + $
;                           0.5 * (variance(tbts[*,0], /nan) + variance(tbts[*,1], /nan)))

;        imasses[ii,jj] = alog10(total(10.^
;        eimasses[ii,jj] = 1./alog(10) * emaxmass/maxmass

        tbt = sprofs[value_locate(rad*rp,2.5),*] / mprofs[-1,*]
        ttbt = median(tbt)
        etbt = stddev(tbt, /nan)

        savedata[ii].ID = files[ii]
        savedata[ii].OMASS[jj] = mmaster
        savedata[ii].EOMASS[jj] = emmaster
        savedata[ii].INTMASS[jj] = maxmass
        savedata[ii].EINTMASS[jj] = emaxmass
        savedata[ii].BT_O[jj] = bts[ii,jj]
        savedata[ii].EBT_O[jj] = ebts[ii,jj]
        savedata[ii].BT_INT[jj] = ttbt
        savedata[ii].EBT_INT[jj] = etbt

        if jj eq 0 then begin
           sbts[ii,*] = getfluxfrac(files[ii]+'_f160w.fits', data[qui].Z[0], mag)
           savedata[ii].SBTS = sbts[ii,*]
        endif
        
     endfor
     
  endfor
  masses  = mean(masses, dim = 3)
  emasses = mean(emasses, dim = 3)
;  sfrs    = mean(sfrs, dim = 2, /nan)
  
  pregions = ['OPTFL', 'INNFL', 'INTER', 'OUTER']
  extract = []
  for ii = 0, n_elements(pregions) - 1 do $
     extract = [extract, where(regions eq pregions[ii])]
  resreg = extract[1:*]
;  cols = ['777777'x, '0000ff'x, '00a500'x, 'ff5500'x]

;  stop
  
  ;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;
  ;; PLOT STUFF ;;
  ;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;
  
  set_plot, 'PS'
  xsize = 7.5;7
  ysize = 6;7
  ar = ysize / xsize
  device, filename = 'sizeMassPlus_biglabels.eps', $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = xsize, ysize = ysize, /in

                                ; -- Re -- 
  
  plot, masses[*,extract[0]], alog10(ssizes[*,1]), /nodat, $
        psym = 8, $
        xran = $;minmax(masses[*,extract[0]])+$
        [10.4,11.49], /xsty, $
        xtickint = 0.5, $
        xminor = 5, $
        yran = [-0.3,0.8], $
        xthick = 4, ythick = 4, $
        charsize = 1.6, charthick = 5, $
;        xtitle = 'log !18M!X!D*!N [M!D'+sunsymbol()+'!N]', $
;        ytitle = 'log !18r!X!De, intrinsic!N [kpc]', $
        pos = [0.15,0.20,0.50,0.9] + [0,0,0.025,0]
  cgtext, 0.06, mean(!Y.WINDOW), $
          'log !18r!X!De, intrinsic!N [kpc]', $
          charsize = 1.7, charthick = 5, align = 0.5, orien = 90, /norm
  ;; Add size-mass fits from vdW14 @ z = 1.25 + resolution limits
  res = alog10(0.09 / zang(1.0, mean(summaryData.Z), $
                           h0 = 73., omega_m = 0.27, lambda0 = 0.73))

  mmm = findgen((!X.CRANGE[1] - !X.CRANGE[0]) / 0.1 + 2) * 0.1 + !X.CRANGE[0]
  xxx    = [mmm, reverse(mmm)]
  yyyred = [0.22 + 0.76 * (mmm - 10.7) - 0.2, $
            reverse(0.22 + 0.76 * (mmm - 10.7) + 0.2)]
  yyyblue= [0.70 + 0.22 * (mmm - 10.7) - 0.17, $
            reverse(0.70 + 0.22 * (mmm - 10.7) + 0.17)]
  polyfill, xxx, yyyblue < !Y.CRANGE[1], col = 'ffa500'x, $
            /line_fill, thick = 1, spacing = 0.05, orien = 45
  polyfill, xxx, yyyred < !Y.CRANGE[1], col = '00a5ff'x, $
            /line_fill, thick = 1, spacing = 0.05, orien = -45
  xxx = [!X.CRANGE, reverse(!X.CRANGE)]
  yyy = [!Y.CRANGE[0], !Y.CRANGE[0], res, res]
  polyfill, xxx, yyy, /line_fill, $
            col = '555555'x, orien = 45, thick = 1, spacing = 0.05

  cgloadct, 33, /rev, ncolors = 9, $
            clip = [10,240]
  ssfrs = alog10(sfrs/10.^masses[*,extract[0]])
  minssfr = -11
  maxssfr = -8
  dssfr   = (maxssfr - minssfr) / 8.
  tssfrs  = findgen(9)*dssfr + minssfr
  cols  = []
  for ii = 0, ngals - 1 do begin
     foo = min(abs(tssfrs - ssfrs[ii]), hit);value_locate(tssfrs, ssfrs[ii]);
     cols = [cols, hit[0]]
  endfor

;  print, masses[*,extract[0]]
;  print, ssfrs
  for jj = 0, 0 do begin
     case jj of
        0: begin
           plotsym, 0, /fill
           lsty = 0
        end
        1: begin
           plotsym, 0, thick = 2
           lsty = 1
        end
     endcase
     
     q  = where(regions eq 'OPTFL')
     tsizes   = sizes[*,jj]
     tssizes  = ssizes[*,jj]
     tesizes  = esizes[*,jj]
     tessizes = essizes[*,jj]
     
     tmasses  = masses[*,q]
     temasses = emasses[*,q]

     oploterror, tmasses, tsizes, $
                 temasses, tesizes, $
                 psym = 8, $
;                 col = '777777'x, $
;                 errcol = '777777'x, $
                 symsize = 1.75, /nohat, $
                 thick = 4, errthick = 4
     oploterror, tmasses, tssizes, $
                 temasses, tessizes, $
                 psym = 1, $
;                 col = '777777'x, $
;                 errcol = '777777'x, $
                 symsize = 1.75, $
                 errsty = 1, /nohat, $
                 thick = 4, errthick = 4

     plotsym, 0, /fill          ;
     for ii = 0, ngals - 1 do begin
        oplot, replicate(tmasses[ii],2), [tsizes[ii],tssizes[ii]], $
               col = '777777'x, linesty = 1, thick = 4
        oplot, [tmasses[ii],tmasses[ii]+alog10(mags[ii])], $
               !Y.CRANGE[0]+0.1+[0, 0.5 * alog10(mags[ii])], $
               thick = 6, col = '770077'x, linesty = 1
        oplot, tmasses[ii]+emags[ii]*[-1,1], $
               (!Y.CRANGE[0]+0.1) + emags[ii]*[-1,1]/2., thick = 2, col = '0055ff'x
        oplot, [tmasses[ii]+alog10(mags[ii])], $
               [!Y.CRANGE[0]+0.1+0.5 * alog10(mags[ii])], $
               psym = 8, symsize = 1, col = cgcolor(string(cols[ii]))
        cgtext, tmasses[ii], tsizes[ii] - 0.1, names[ii], $
                charsize = 1.3, charthick = 4, align = 0.5
        if tmasses[ii] eq min(tmasses) then $
           cgtext, mean([tmasses[ii]-emags[ii],tmasses[ii]+alog10(mags[ii])]), $
                   mean(!Y.CRANGE[0]+0.1+[0, 0.5 * alog10(mags[ii])])+0.05, $
                   greek('mu'), col = '770077'x, charthick = 3, orien = atan(ar) * 180/!pi, $
                   charsize = 1.3

        oplot, [tmasses[ii]], [tsizes[ii]], psym = 8, $
               symsize = 1.25, col = cgcolor(string(cols[ii]))
        oplot, [tmasses[ii]], [tssizes[ii]], psym = 1, $
               symsize = 1.25, col = cgcolor(string(cols[ii]))

     endfor
  
  endfor
  cgcolorbar, pos = [0.175,0.825,0.375,0.85], $
              charsize = 1.2, charthick = 4, $
              minrange = minssfr, maxrange = maxssfr, $
              title = 'log !18sSFR!X', thick = 3, /discrete, $
              xtickint = 1, textthick = 3, ncolors = 9;, $
;              oob_low = 240, oob_hi = 10
  legend, pos = [10.95,0.05], /data, box = 0, $;/horiz, $
          reverse(['1D (used)', '2D (SE)']), $
          psym = reverse([8, 1]), symsize = [1.5,1.5], thick = [4,1], $;col = '777777'x, $
          charsize = 1.3, charthick = 4, pspacing = 0.5, spacing = 1.5

                                ; -- B/T -- 

  if NOT circularize then $
     yttl = '!18M!X!D*,Inner!N!18/M!X!D*,sum!N' $
  else $
     yttl = '!18M!X(!18r!X<2.5 kpc)!18/M!X!D*!N'

  xw = !X.WINDOW[1] - !X.WINDOW[0]
  plotsym, 0, /fill
  plot, masses[*,extract[0]], bts[*,1], /nodat, $
        psym = 8, $
        xran = !X.CRANGE + [1d-6,0], /xsty, $
        yran = [(min(bts)-0.2) > 0,max(bts)+0.2], $
        xtickint = 0.5, $
        xminor = 5, $
        ysty = 8+1, $
        ytickname = replicate(' ',60), $
        xthick = 4, ythick = 4, $
        charsize = 1.6, charthick = 5, $
;        xtitle = 'log !18M!X!D*!N [M!D'+sunsymbol()+'!N]', $
        pos = [!X.WINDOW[1],!Y.WINDOW[0],!X.WINDOW[1] + xw,!Y.WINDOW[1]], /noer, $
        yminor = 2, ytickint = 0.1
  axis, yaxis = 1, $;ytitle = yttl, $; ["!18B/T!X"]', $
        ythick = 4, charsize = 1.5, charthick = 5, $
        yminor = 2, ytickint = 0.1, yran = !Y.CRANGE, /ysty
  cgtext, 0.985, mean(!Y.WINDOW), $
          yttl, $
          charsize = 1.7, charthick = 5, align = 0.5, orien = 90, /norm
  
  ;; Read your SDSS shit
  local = mrdfits('bt_stuff_for_dan.fits', 1)
  lm = local.LOC + 0.15
  lllx = ([lm, reverse(lm)] > !X.CRANGE[0]) < !X.CRANGE[1]
  llly = ([local.P25, reverse(local.P75)] > !Y.CRANGE[0]) < !Y.CRANGE[1]
  polyfill, lllx, llly, col = cgcolor('pink');, /line_fill, $
;            thick = 1, orien = -45, spacing = 0.05
  
  q  = where(regions eq 'OPTFL')
  tmasses  = reform(masses[*,q])
  temasses = reform(emasses[*,q])
  oploterror, tmasses, bts[*,1], $
              temasses, ebts[*,1], $
              psym = 8, symsize = 1.25, $
              col = '777777'x, errcol = '777777'x, $
              /nohat, errthick = 4
  oploterror, tmasses, bts[*,0], $
              temasses, ebts[*,0], $
              psym = 8, symsize = 1.75, $
              /nohat, errthick = 4
  for ii = 0, ngals - 1 do begin
;     oplot, replicate(tmasses[ii], 5), sbts[ii,*], psym = 1
     oplot, [tmasses[ii]], [bts[ii,1]], psym = 8, $
            symsize = 1, col = cgcolor(string(cols[ii]))
     oplot, [tmasses[ii]], [bts[ii,0]], psym = 8, $
            symsize = 1.25, col = cgcolor(string(cols[ii]))
     oploterror, tmasses[ii], reform(sbts[ii,2]), $
                 temasses[ii], reform([sbts[ii,3] - sbts[ii,2]]), $
                 psym = 1, errsty = 1, symsize = 1.75, errthick = 1, thick = 6, /nohat
     oplot, [tmasses[ii]], reform([sbts[ii,2]]), $
            psym = 1, symsize = 1.25, thick = 2, $
            col = cgcolor(string(cols[ii]))
     plotsym, 0, /fill          ;
;     oploterror, [imasses[ii,0]], [bts[ii,0]], $
;                 [eimasses[ii,0]], [ebts[ii,0]], /nohat, $
;                 psym = 8, symsize = 1, $
;                 col = cgcolor(string(cols[ii])), errcol = cgcolor(string(cols[ii]))
;     oplot, [imasses[*,1]], [bts[*,1]], psym = 8, symsize = 1, col = cgcolor(string(cols[ii]))
  endfor
  for ii = 0, ngals - 1 do $
     oplot, [tmasses[ii],tmasses[ii]+alog10(mags[ii])], $
            !Y.CRANGE[0]+0.015+[0,0] + 0.01 * ii, $
            thick = 6, col = '770077'x, linesty = 1
  for ii = 0, ngals - 1 do $
     oplot, tmasses[ii]+emags[ii]*[-1,1], $
            !Y.CRANGE[0]+0.015+[0,0] + 0.01 * ii, thick = 2, col = '0055ff'x
  for ii = 0, ngals - 1 do $
     oplot, [tmasses[ii]+alog10(mags[ii])], [!Y.CRANGE[0]+0.015] + 0.01 * ii, $
            psym = 8, symsize = 0.75, col = cgcolor(string(cols[ii]))
  plot, masses[*,extract[0]], bts[*,1], /nodat, $
        psym = 8, $
        xran = !X.CRANGE + [1d-6,0], /xsty, $
        yran = [(min(bts)-0.2) > 0,max(bts)+0.2], $
        xtickint = 0.5, $
        xminor = 5, $
        ysty = 8+1, $
        xtickname = replicate(' ',60), $
        ytickname = replicate(' ',60), $
        xthick = 4, ythick = 4, $
        pos = [!X.WINDOW[0], !Y.WINDOW[0], !X.WINDOW[1], !Y.WINDOW[1]], /noer, yminor = 2, ytickint = 0.1
  axis, yaxis = 1, $
        ythick = 4, ytickname = replicate(' ', 60), $
        yminor = 2, ytickint = 0.1, yran = !Y.CRANGE, /ysty

  plotsym, 0, /fill
  legend, pos = [!X.crange[0], !Y.crange[1]-0.02], /data, box = 0, $;/top, /left, box = 0, $
          ['LogNormal SFH', $
           'DelayedExp SFH', $
           'SE F160W flux ratio', $
           '!18z!X=0 !18B/T!X (Abramson+14)'], $
          charsize = 1.1, charthick = 4, $
          col = [0,'777777'x, 0, cgcolor('pink')], psym = [8,8,1,8], $
          pspacing = 0.5, spacing = 1.5, $
          symsize = [1.75, 1.25, 1.25, 1.75], thick = [1,1,4,1]

  cgtext, mean([0.15,!X.WINDOW[1]]), !Y.WINDOW[0] - 0.125, $
          'log !18M!X!D*!N [M!D'+sunsymbol()+'!N]', $
          charsize = 1.7, charthick = 5, align = 0.5, /norm
  
  device, /close
  set_plot, 'X'
;  spawn, 'gv sizeMassPlus_biglabels.eps &'

  mwrfits, savedata, 'sizeMassSummaryData.fits', /create
  
end

;plotsizemass_biglabels, 'study5_lgn_Bestfits.list', 'study5_exp_Bestfits.list', /circularize

