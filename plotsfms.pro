pro plotsfms, lgnlist, explist, $
              NSIG = nsig

  if NOT keyword_set(NSIG) then nsig = 1

  summaryData = mrdfits('objectSummaryData.fits', 1)
  finalMasses = mrdfits('finalMasses.fits', 1)
;  magerr = 2./alog(10) * summaryData.EMU/summaryData.MU
  
  readcol, explist, expfiles, f = 'A'
  readcol, lgnlist, lgnfiles, f = 'A'

  ngals = n_elements(lgnfiles)

  test     = mrdfits(lgnfiles[0], 1)

  regions  = test.REGION
  nregions = n_elements(regions)
  
  sfrs    = fltarr(ngals, nregions, 2) ;; Lgn, Exp
  esfrs   = fltarr(ngals, nregions, 2)
  lsfrs   = fltarr(ngals, nregions, 2) ;; for Upper limits
  
  masses  = fltarr(ngals, nregions, 2)
  emasses = fltarr(ngals, nregions, 2)

  files = strarr(ngals)
  for ii = 0, ngals - 1 do begin
     for jj = 0, 1 do begin
        case jj of
           0: begin
              file = lgnfiles[ii]
              print, file
              files[ii] = strmid(file, 0, 7)
              qui = where(summaryData.ID eq files[ii])
              mag = summaryData[qui].MU
              mag = mag[0]
              emag = summaryData[qui].EMU
              emag = emag[0]
              print, mag
           end
           1: file = expfiles[ii]
        endcase
                
        data = mrdfits(file, 1)

        ;; Correct for magnification
        ;; Typical error in mag is <0.1 dex so forget it
        data.LMASS -= alog10(mag)
        data.SFR   /= mag
        
        masses[ii,*,jj]  = data.LMASS[0]
        emasses[ii,*,jj] = 0.5 * (data.LMASS[2] - data.LMASS[1]) ;; mag errs added below
        
        sfrs[ii,*,jj]    = data.SFR[0]
        esfrs[ii,*,jj]   = 0.5 * (data.SFR[2] - data.SFR[1]) ;; mag errs added below

;        uls = where(alog10(sfrs[ii,*,jj]) lt nsig * esfrs[ii,*,jj] /
;        alog(10) / sfrs[ii,*,jj], nuls)
        uls = where(sfrs[ii,*,jj] lt nsig * (data.SFR[0] - data.SFR[1]), nuls) ;; lower error bar
        if nuls gt 0 then begin
           sfrs[ii,uls,jj] = 0
           lsfrs[ii,uls,jj] = (data.SFR[0] - data.SFR[1])[uls]
           esfrs[ii,uls,jj] = (data.SFR[0] - data.SFR[1])[uls]
        endif
           
     endfor
  endfor


  pregions = ['OPTFL', 'INNFL', 'INTER', 'OUTER']
  extract = []
  for ii = 0, n_elements(pregions) - 1 do $
     extract = [extract, where(regions eq pregions[ii])]
  resreg = extract[1:*]
  cols = ['777777'x, '0000ff'x, '00a500'x, 'ff5500'x]

  fmasses  = finalMasses.LMASS  ;; from getfinalmasses.pro
  efmasses = finalMasses.ELMASS

  masses[*,extract[0],0] = fmasses
  masses[*,extract[0],1] = fmasses
  emasses[*,extract[0],0] = efmasses
  emasses[*,extract[0],1] = efmasses

  ;; SFMS from whitaker+14 @ 1 < z < 1.5
  a  = -26.03 - 0.15 ;; to Salpeter
  ae = 1.69
  b  = 4.62
  be = 0.34
  c  = -0.19
  ce = 0.02

  dm = 0.1
  m1 = 11.5
  m0 = 9
  m = findgen((m1-m0)/dm+1) * dm + m0 - 0.15
  sfms = a + b * m + c * m^2
  
  plot, fmasses, alog10(sfrs[*,*,1]), /nodat, $
        psym = 8, $
        xran = [9.5,11.5], $
        yran = [-1,2], $
        xtitle = 'log !18M!X!D*!N [M!D'+sunsymbol()+'!N]', $
        ytitle = 'log !18SFR!X [M!D'+sunsymbol()+'!N yr!E-1!N]'
  q = where(m gt !X.CRANGE[0] and m lt !X.crange[1]) ; AND sfms gt !Y.CRANGE[0] AND sfms lt !Y.CRANGE[1])
  polyfill, [m[q], reverse(m[q])], [sfms[q]-0.3, reverse(sfms[q]+0.3)] < !Y.CRANGE[1], $ ;col = '555555'x
            /line_fill, orien = -45, thick = 1, spacing = 0.025, col = sfmscol
  for ii = 0, n_elements(pregions) - 1 do begin

     for jj = 0, 1 do begin
        case jj of
           0: begin
              plotsym, 0, /fill
              lsty = 0
              ec = 0
           end
           1: begin
              plotsym, 0, thick = 2
              lsty = 2
              ec = '777777'x
           end
        endcase

        q  = where(regions eq pregions[ii])
        tsfrs    = sfrs[*,q,jj]
        tesfrs   = esfrs[*,q,jj]
        tmasses  = masses[*,q,jj]
        temasses = emasses[*,q,jj]
        
        uq   = where(tsfrs gt 0) 
        lims = where(tsfrs eq 0, nlims) 
;        print, alog10(nsig*tesfrs[lims])

        for kk = 0, ngals - 1 do begin
           oplot, masses[kk,extract,0], alog10(sfrs[kk,extract,0] > nsig * lsfrs[kk,extract,0]), linesty = 0
           oplot, masses[kk,extract,1], alog10(sfrs[kk,extract,1] > nsig * lsfrs[kk,extract,1]), linesty = 2
        endfor

        oplot, [tmasses[uq]], [alog10(tsfrs[uq])], psym= 8, symsize = 1.25, col = ec
        oploterror, [tmasses[uq]], [alog10(tsfrs[uq])], $
                    [temasses[uq]], [1./alog(10) * tesfrs[uq] / tsfrs[uq]], $
                    psym = 8, col = cols[ii], errcol = ec, symsize = 1., /nohat
        if nlims gt 0 then $
           for kk = 0, nlims - 1 do $
              oplot, tmasses[lims[kk]] + temasses[lims[kk]] * [-1,1], $
                     replicate(alog10(nsig * tesfrs[lims[kk]]), 2), $
                     col = cols[ii], linesty = lsty
     endfor
  endfor
  
;  !p.multi = [0,ngals,0]
;  window, 0, xsize = 1200, ysize = 300
;  xs = 0.1
;  ys = 0.25
;  yf = 0.95
;  xw = (1 - 2*xs) / ngals

  xm = 0.15
  xw = (1 - 1.5 * xm) / 2.
  ym = xm
  yw = (1 - 1.5 * ym) / 2.
  pos1 = [xm,ym+yw,xm+xw,ym+2*yw]
  pos2 = [xm+xw,ym+yw,xm+2*xw,ym+2*yw]
  pos3 = [xm,ym,xm+xw,ym+yw]
  pos4 = [xm+xw,ym,xm+2*xw,ym+yw]
  
  set_plot, 'PS'
;  xsize = 8.5
  xsize = 5.
  ysize = 5.
  device, filename = 'resolvedSFMS.eps', $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = xsize, ysize = ysize, /in
  sfmscol = 'aaaaaa'x
  yr = [-0.3,2]
  for kk = 0, ngals - 1 do begin
     
     case kk of
        0: plot, masses[kk,*,*], alog10(sfrs[kk,*,*] > esfrs[kk,*,*]) > 0, /nodat, $
                 xran = [9.5,11.5-1d-6], psym = 8, yran = yr, $ ;yr = [minsfr,maxsfr], $
                 /xsty, $
;              xtitle = 'log !18M!X!D*!N [M!D'+sunsymbol()+'!N]', $
;              ytitle = 'log !18SFR!X [M!D'+sunsymbol()+'!N yr!E-1!N]', $
                 xthick = 3, ythick = 3, charsize = 1.25, charthick = 4, $
                 pos = pos1, $
                 xtickname = replicate(' ', 60), /ysty
        1: plot, masses[kk,*,*], alog10(sfrs[kk,*,*] > esfrs[kk,*,*]) > 0, /nodat, $
                 xran = [9.5+1d-6,11.5], psym = 8, yran = yr, $ ;yr = [minsfr,maxsfr], $
                 /xsty, $
;              xtitle = 'log !18M!X!D*!N [M!D'+sunsymbol()+'!N]', $
;              ytitle = 'log !18SFR!X [M!D'+sunsymbol()+'!N yr!E-1!N]', $
                 ytickname = replicate(' ', 60), $
                 xtickname = replicate(' ', 60), $
                 xthick = 3, ythick = 3, charsize = 1.25, charthick = 4, $
                 pos = pos2, /noer
        2: plot, masses[kk,*,*], alog10(sfrs[kk,*,*] > esfrs[kk,*,*]) > 0, /nodat, $
                 xran = [9.5,11.5-1d-6], psym = 8, yran = yr, $ ;yr = [minsfr,maxsfr], $
                 /xsty, $
;              xtitle = 'log !18M!X!D*!N [M!D'+sunsymbol()+'!N]', $
;              ytitle = 'log !18SFR!X [M!D'+sunsymbol()+'!N yr!E-1!N]', $
;                 ytickname = replicate(' ', 60), $
                 xthick = 3, ythick = 3, charsize = 1.25, charthick = 4, $
                 pos = pos3, /noer, /ysty
        3: plot, masses[kk,*,*], alog10(sfrs[kk,*,*] > esfrs[kk,*,*]) > 0, /nodat, $
                 xran = [9.5-1d-6,11.5], psym = 8, yran = yr, $ ;yr = [minsfr,maxsfr], $
                 /xsty, $
;              xtitle = 'log !18M!X!D*!N [M!D'+sunsymbol()+'!N]', $
;              ytitle = 'log !18SFR!X [M!D'+sunsymbol()+'!N yr!E-1!N]', $
                 ytickname = replicate(' ', 60), $
                 xthick = 3, ythick = 3, charsize = 1.15, charthick = 4, $
                 pos = pos4, /noer, /ysty
     endcase
     q = where(m gt !X.CRANGE[0] and m lt !X.crange[1]) ; AND sfms gt !Y.CRANGE[0] AND sfms lt !Y.CRANGE[1])
     polyfill, [m[q], reverse(m[q])], [sfms[q]-0.3, reverse(sfms[q]+0.3)] < !Y.CRANGE[1], $;col = '555555'x
               /line_fill, orien = -45, thick = 1, spacing = 0.025, col = sfmscol
;     oplot, m[q], sfms[q], thick = 6, col = 'ffa500'x, linesty = 2
     oplot, masses[kk,extract[1:*],0], alog10(sfrs[kk,extract[1:*],0] > nsig * lsfrs[kk,extract[1:*],0]), $
            linesty = 0, thick = 6
     oplot, masses[kk,extract[1:*],1], alog10(sfrs[kk,extract[1:*],1] > nsig * lsfrs[kk,extract[1:*],1]), $
            linesty = 2, thick = 6

     q = where(summaryData.ID eq files[kk])
     mu  = summaryData[q].MU
     emu = 1./alog(10) * summaryData[q].EMU / mu
     xctr = 11.
     yctr = 0.3
     sx = xctr - emu
     ex = xctr + emu
     sy = yctr - emu
     ey = yctr + emu
     oplot, [xctr,xctr+alog10(mu)], [yctr, yctr+alog10(mu)], thick = 6, col = '770077'x, linesty = 1
     oplot, [sx,ex], [sy,ey], thick = 2, col = '0055ff'x

;     oploterror, [alog10(total(10.^masses[kk,resreg,0]))], [alog10(total(sfrs[kk,resreg,0]))], $
;                 sqrt(total(emasses[kk,resreg,0]^2)), sqrt(total((0.5/alog(10) * esfrs[kk,resreg,0]/sfrs[kk,resreg,0])^2)), $
;                 psym = 2, col = 'ffffff'x, errcol = 'ffffff'x
;     oploterror, [alog10(total(10.^masses[kk,resreg,1]))], [alog10(total(sfrs[kk,resreg,1]))], $
;                 sqrt(total(emasses[kk,resreg,1]^2)), sqrt(total((0.5/alog(10) * esfrs[kk,resreg,1]/sfrs[kk,resreg,1])^2)), $
;                 psym = 4, col = 'ffffff'x, errcol = 'ffffff'x
     for ii = 0, n_elements(pregions) - 1 do begin
        for jj = 1,0,-1 do begin
           case jj of
              0: begin
                 plotsym, 0, /fill
                 lsty = 0
                 ec = 0
              end
              1: begin
                 plotsym, 0, thick = 6
                 lsty = 2
                 ec = '777777'x
              end
           endcase
           q  = where(regions eq pregions[ii])
           tsfrs    = (sfrs[kk,q,jj])[0]
           tesfrs   = sqrt((esfrs[kk,q,jj])[0]^2)
           tlsfrs   = (lsfrs[kk,q,jj])[0]
           tmasses  = (masses[kk,q,jj])[0]
           temasses = sqrt((emasses[kk,q,jj])[0]^2 + emu^2)

           if tsfrs gt 0 then begin
              oploterror, [tmasses], [alog10(tsfrs)], $
                          [temasses], [sqrt((0.5/alog(10) * tesfrs / tsfrs)^2 + emu^2)], $
                          psym = 8, col = ec, errcol = ec, symsize = 1.25, $
                          errthick = 4, /nohat
              oplot, [tmasses], [alog10(tsfrs)], psym = 8, symsize = 1, col = cols[ii]
           endif else begin
              oplot, tmasses + temasses * [-1,1], $
                     replicate(alog10(nsig * tlsfrs) > (!Y.CRANGE[0] + 0.2), 2), $
                     col = cols[ii], linesty = lsty, thick = 6
              plotsym, 1, thick = 4
              oplot, [tmasses], $
                     [alog10(nsig * tlsfrs) > (!Y.CRANGE[0] + 0.2)], $
                     col = cols[ii], psym = 8, symsize = 1.5, thick = 4
           endelse
        endfor
        cgtext, !X.CRANGE[0]+0.1, 0.9*!Y.CRANGE[1], /data, $
                files[kk], charsize = 1, charthick = 3;, align = 1
        if kk eq 0 then begin
           legend, pos = [!X.CRANGE[0], !Y.CRANGE[0]+0.75], /data, /left, box = 0, $
                   ['LogNormal', 'DelayedExp'], $
                   linesty = [0,2], thick = [4,4], $
                   pspacing = 1, charsize = 1, charthick = 3
        endif else if kk eq 1 then begin
           ar = $;((!X.CRANGE[1] - !X.CRANGE[0]) / (!Y.CRANGE[1] - !Y.CRANGE[0])) * $
                xw / (0.25 * (!Y.WINDOW[1] - !Y.WINDOW[0]))
           cgtext, charsize = 1, charthick = 2, greek('mu'), $
                   align = 0.5, xctr + 0.025, yctr+0.2, $
                   /data, orien = atan(1./ar) * 180. / !pi, col = '770077'x
        endif else if kk eq 2 then begin
           x0 = 10.6;!X.CRANGE[1] - 0.5
           x1 = 11.0;!X.CRANGE[1] - 0.1
           y0 = !Y.CRANGE[0] + 0.15
           y1 = !Y.CRANGE[0] + 0.25
           polyfill, [x0,x1,x1,x0], [y0,y0,y1,y1], $
                     /line_fill, orien = 45, thick = 1, spacing = 0.025, col = sfmscol
;           oplot, [x0,x1], replicate(mean([y0,y1]), 2), thick = 4, linesty = 2, col = 'ffa500'x
           cgtext, charsize = 1, charthick = 3, 'W+14 SFMS', $
                   align = 1, x0-0.1, y0, /data
        endif else if kk eq 3 then begin
           plotsym, 0, /fill
           legend, /bottom, /left, box = 0, $ ;pos = [!X.CRANGE[0], !Y.CRANGE[0] + 1]
                   ['Opt', 'Inn', 'Mid', 'Out'], $
                   psym = 8, col = cols, $
                   pspacing = 0.5, charsize = 1, charthick = 3
        endif

     endfor

  endfor
  cgtext, 0.05, ym+yw, 'log !18SFR!X [M!D'+sunsymbol()+'!N yr!E-1!N]', $
          charsize = 1.25, charthick = 4, align = 0.5, /norm, orien = 90
  cgtext, xm+xw, 0.05, 'log !18M!X!D*!N [M!D'+sunsymbol()+'!N]', $
          charsize = 1.25, charthick = 4, align = 0.5, /norm
  device, /close

  ;;;;;;;;;;
  ;;;;;;;;;;
  ;; SSFR ;;
  ;;;;;;;;;;
  ;;;;;;;;;;
  
  xm = 0.15
  xw = (1 - 1.5 * xm) / 2.
  ym = 0.15
  yw = (1 - 1.5 * ym) / 2.
  pos1 = [xm,ym+yw,xm+xw,ym+2*yw]      + [0.025,0,0.025,0]
  pos2 = [xm+xw,ym+yw,xm+2*xw,ym+2*yw] + [0.025,0,0.025,0]
  pos3 = [xm,ym,xm+xw,ym+yw]           + [0.025,0,0.025,0]
  pos4 = [xm+xw,ym,xm+2*xw,ym+yw]      + [0.025,0,0.025,0]
  
  device, filename = 'resolvedSSFMS.eps', $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = xsize, ysize = ysize, /in
  sfmscol = 'aaaaaa'x
  yr = -9 - [2,-0.3]
  for kk = 0, ngals - 1 do begin
     
     case kk of
        0: plot, masses[kk,*,*], (alog10(sfrs[kk,*,*] > esfrs[kk,*,*]) > 0) - masses[kk,*,*], /nodat, $
                 xran = [9.5,11.5-1d-6], psym = 8, yran = yr, $ ;yr = [minsfr,maxsfr], $
                 /xsty, $
;              xtitle = 'log !18M!X!D*!N [M!D'+sunsymbol()+'!N]', $
;              ytitle = 'log !18SFR/M!X!D*!N [yr!E-1!N]', $
                 xthick = 3, ythick = 3, charsize = 1.25, charthick = 4, $
                 pos = pos1, /ysty, $
                 xtickname = replicate(' ', 60)
        1: plot, masses[kk,*,*], (alog10(sfrs[kk,*,*] > esfrs[kk,*,*]) > 0) - masses[kk,*,*], /nodat, $
                 xran = [9.5,11.5-1d-6], psym = 8, yran = yr, $ ;yr = [minsfr,maxsfr], $
                 /xsty, $
;              xtitle = 'log !18M!X!D*!N [M!D'+sunsymbol()+'!N]', $
;              ytitle = 'log !18SFR!X [M!D'+sunsymbol()+'!N yr!E-1!N]', $
                 ytickname = replicate(' ', 60), $
                 xthick = 3, ythick = 3, charsize = 1.25, charthick = 4, $
                 pos = pos2, /ysty, $
                 xtickname = replicate(' ', 60), /noer
        2: plot, masses[kk,*,*], (alog10(sfrs[kk,*,*] > esfrs[kk,*,*]) > 0) - masses[kk,*,*], /nodat, $
                 xran = [9.5,11.5-1d-6], psym = 8, yran = yr, $ ;yr = [minsfr,maxsfr], $
                 /xsty, $
;              xtitle = 'log !18M!X!D*!N [M!D'+sunsymbol()+'!N]', $
;              ytitle = 'log !18SFR!X [M!D'+sunsymbol()+'!N yr!E-1!N]', $
                 xthick = 3, ythick = 3, charsize = 1.25, charthick = 4, $
                 pos = pos3, /noer, /ysty
        3: plot, masses[kk,*,*], (alog10(sfrs[kk,*,*] > esfrs[kk,*,*]) > 0) - masses[kk,*,*], /nodat, $
                 xran = [9.5,11.5], psym = 8, yran = yr, $ ;yr = [minsfr,maxsfr], $
                 /xsty, $
;              xtitle = 'log !18M!X!D*!N [M!D'+sunsymbol()+'!N]', $
;              ytitle = 'log !18SFR!X [M!D'+sunsymbol()+'!N yr!E-1!N]', $
                 ytickname = replicate(' ', 60), $
                 xthick = 3, ythick = 3, charsize = 1.25, charthick = 4, $
                 pos = pos4, /noer, /ysty
     endcase
     
     q = where(m gt !X.CRANGE[0] and m lt !X.crange[1]) ; AND sfms gt !Y.CRANGE[0] AND sfms lt !Y.CRANGE[1])
     polyfill, [m[q], reverse(m[q])], [sfms[q]-0.3-m[q], reverse(sfms[q]+0.3-m[q])] < !Y.CRANGE[1], $;col = '555555'x
               /line_fill, orien = -45, thick = 1, spacing = 0.025, col = sfmscol
;     oplot, m[q], sfms[q], thick = 6, col = 'ffa500'x, linesty = 2
     oplot, masses[kk,extract[1:*],0], alog10(sfrs[kk,extract[1:*],0] > nsig * lsfrs[kk,extract[1:*],0]) - masses[kk,extract[1:*],0], $
            linesty = 0, thick = 6
     oplot, masses[kk,extract[1:*],1], alog10(sfrs[kk,extract[1:*],1] > nsig * lsfrs[kk,extract[1:*],1]) - masses[kk,extract[1:*],0], $
            linesty = 2, thick = 6

     q = where(summaryData.ID eq files[kk])
     mu  = summaryData[q].MU
     emu = 1./alog(10) * summaryData[q].EMU / mu
     xctr = 11.
     yctr = -9 - 1.5
     sx = xctr - emu
     ex = xctr + emu
     sy = yctr
     ey = yctr
     oplot, [xctr,xctr+alog10(mu)], [yctr, yctr], thick = 6, col = '770077'x, linesty = 1
     oplot, [sx,ex], [sy,ey], thick = 2, col = '0055ff'x

     for ii = 0, n_elements(pregions) - 1 do begin
        for jj = 1,0,-1 do begin
           case jj of
              0: begin
                 plotsym, 0, /fill
                 lsty = 0
                 ec = 0
              end
              1: begin
                 plotsym, 0, thick = 6
                 lsty = 2
                 ec = '777777'x
              end
           endcase
           q  = where(regions eq pregions[ii])
           tsfrs    = (sfrs[kk,q,jj])[0]
           tesfrs   = sqrt((esfrs[kk,q,jj])[0]^2 + emu^2)
           tlsfrs   = (lsfrs[kk,q,jj])[0]
           tmasses  = (masses[kk,q,jj])[0]
           temasses = sqrt((emasses[kk,q,jj])[0]^2 + emu^2)

           if tsfrs gt 0 then begin
              oploterror, [tmasses], [alog10(tsfrs) - tmasses], $
                          [temasses], [sqrt((0.5/alog(10) * tesfrs / tsfrs)^2 + temasses^2 + 2*emu^2)], $
                          psym = 8, col = ec, errcol = ec, symsize = 1.25, $
                          errthick = 4, /nohat
              oplot, [tmasses], [alog10(tsfrs) - tmasses], psym = 8, symsize = 1, col = cols[ii]
           endif else begin
              oplot, tmasses + temasses * [-1,1], $
                     replicate((alog10(nsig * tlsfrs) - tmasses) > (!Y.CRANGE[0] + 0.2), 2), $
                     col = cols[ii], linesty = lsty, thick = 6
              plotsym, 1, thick = 4
              oplot, [tmasses], $
                     [(alog10(nsig * tlsfrs) - tmasses) > (!Y.CRANGE[0] + 0.2)], $
                     col = cols[ii], psym = 8, symsize = 1.5, thick = 4
           endelse
        endfor
        cgtext, !X.CRANGE[1]-0.1, !Y.CRANGE[1]-0.25, /data, $
                files[kk], charsize = 1, charthick = 3, align = 1
        if kk eq 0 then begin
           legend, pos = [!X.CRANGE[0], !Y.CRANGE[0]+0.75], /data, /left, box = 0, $
                   ['LogNormal', 'DelayedExp'], $
                   linesty = [0,2], thick = [4,4], $
                   pspacing = 1, charsize = 1, charthick = 3
        endif else if kk eq 1 then begin
           cgtext, charsize = 1, charthick = 2, greek('mu'), $
                   align = 0.5, ex, yctr+0.1, $
                   /data, col = '770077'x
        endif else if kk eq 2 then begin
           x0 = 10.6;!X.CRANGE[1] - 0.5
           x1 = 11.0;!X.CRANGE[1] - 0.1
           y0 = !Y.CRANGE[0] + 0.15
           y1 = !Y.CRANGE[0] + 0.25
           polyfill, [x0,x1,x1,x0], [y0,y0,y1,y1], $
                     /line_fill, orien = 45, thick = 1, spacing = 0.025, col = sfmscol
;           oplot, [x0,x1], replicate(mean([y0,y1]), 2), thick = 4, linesty = 2, col = 'ffa500'x
           cgtext, charsize = 1, charthick = 3, 'W+14 SFMS', $
                   align = 1, x0-0.1, y0, /data
        endif else if kk eq 3 then begin
           plotsym, 0, /fill
           legend, /bottom, /left, box = 0, $ ;pos = [!X.CRANGE[0], !Y.CRANGE[0] + 1]
                   ['Opt', 'Inn', 'Mid', 'Out'], $
                   psym = 8, col = cols, $
                   pspacing = 0.5, charsize = 1, charthick = 3
        endif

     endfor

  endfor
  cgtext, 0.025, ym+yw, 'log !18SFR/M!X!D*!N [yr!E-1!N]', $
          charsize = 1.25, charthick = 4, align = 0.5, /norm, orien = 90
  cgtext, xm+xw, 0.05, 'log !18M!X!D*!N [M!D'+sunsymbol()+'!N]', $
          charsize = 1.25, charthick = 4, align = 0.5, /norm
  
  
  device, /close
  set_plot, 'X'
  !P.multi = 0
  spawn, 'gv resolvedSFMS.eps &'
  spawn, 'gv resolvedSSFMS.eps &'
  
  stop
  
end

;plotsfms, 'study5_lgn_Bestfits.list', 'study5_exp_Bestfits.list', nsig = 2




































;; Plot the observed integrated SFMS and then that fixed to a common
;; reference redshift
pro plotsfmsBAD, inlist

  readcol, inlist, files, f = 'A'
  nfiles = n_elements(files)
  
  masses   = fltarr(nfiles)
  emasses  = fltarr(nfiles)
  sfrs     = fltarr(nfiles)
  esfrs    = fltarr(nfiles)
  sfrs_ha  = fltarr(nfiles)
  esfrs_ha = fltarr(nfiles)
  zs       = fltarr(nfiles)
  ezs      = fltarr(nfiles)
  
  for ii = 0, nfiles - 1 do begin
     data = mrdfits(files[ii], 1)
     opt = where(data.REGION eq 'OPTFL')
     data = data[opt]

     masses[ii]  = data.LMASS[0]
     emasses[ii] = 0.5 * (data.LMASS[2] - data.LMASS[1])
     sfrs[ii]    = data.SFR[0]
     esfrs[ii]   = 0.5 * (data.SFR[2] - data.SFR[1])
;     sfrs_ha[ii,*] = data.HA_FLUX * 1d-18 * dluminosity(zs[ii,0], /cm)^2 * 4 *!pi / 1.26d41
  endfor

  
  rmasses   = fltarr(nfiles, n_elements(data))
  remasses  = fltarr(nfiles, n_elements(data))
  rsfrs     = fltarr(nfiles, n_elements(data))
  for ii = 0, nfiles - 1 do begin
     data = mrdfits(files[ii], 1)
     opt = where(data.REGION ne 'OPTFL')
     data = data[opt]

     if ii eq 0 then regions = data.REGION
     
     rmasses[ii,*]  = data.LMASS[0]
     remasses[ii,*] = 0.5 * (data.LMASS[2] - data.LMASS[1])

     tes = 0.5 * (data.SFR[2] - data.SFR[1])
     ul = where(data.SFR[0] lt 2 * tes, nul, compl = meas)
     rsfrs[ii,meas]  = data[meas].SFR[0]
     resfrs[ii,meas] = tes[meas]

     ;; Add upper limits
     rsfrs[ii,ul] = 2 * tes[ul]
     resfrs[ii,ul] = 2 * tes[ul]
     
     q = where(data.REGION eq 'INNFL')
     zs[ii]  = data[q].Z[0]
     ezs[ii] = 0.5 * (data[q].Z[2] - data[q].Z[1]) 
     

  endfor

;  stop
  
;  z_common = median(zs[*,0]);
;
;  readcol, 'allRecons.list', files, f = 'A'
;  nfiles = n_elements(files)
;  cmasses  = fltarr(nfiles, 3)
;  csfrs    = fltarr(nfiles, 3)
;  for ii = 0, nfiles - 1 do begin
;     data = mrdfits(files[ii], 1)
;     hit = value_locate(data.OPTFL_REDSHIFT, z_common)
;          
;     cmasses[ii,*] = alog10(data.OPTFL_MGH[hit,*])
;     csfrs[ii,*]    = data.OPTFL_SFH[hit,*]
;  endfor
    
  plotsym, 0, /fill
;  !p.multi = [0,2,0]
  plot, masses, alog10(sfrs), /nodat, $
        xran = [9.5,11.5], yran = [-1,2], $
        xtitle = 'log !18M!X!D*!N [M!D'+sunsymbol()+'!N]', $
        ytitle = 'log !18SFR!X [M!D'+sunsymbol()+'!N yr!E-1!N]', /iso
  oploterror, masses, alog10(sfrs), $
              emasses, $
              1./alog(10) * esfrs/sfrs, $
              psym = 8, col = '777777'x, errcol = '777777'x
  cols = [255, '00a500'x, '00ff00'x, 'ff0000'x, 'ffa500'x]
  for ii = 0, 4 do begin
     q = where(rsfrs[*,ii] ne resfrs[*,ii], compl = ul)
     oploterror, [rmasses[q,ii]], [alog10(rsfrs[q,ii])], $
                 [remasses[q,ii]], 1./alog(10) * [resfrs[q,ii]/rsfrs[q,ii]], $
                 psym = 8, col = cols[ii], errcol = cols[ii]
;     if n_elements(ul) gt 0 then stop
     for jj = 0, n_elements(ul) - 1 do begin
        oplot, rmasses[ul[jj],ii] + remasses[ul[jj],ii] * [-1,1], replicate(alog10(resfrs[ul[jj],ii]), 2), $
               col = cols[ii]
        oplot, [rmasses[ul[jj],ii]], [alog10(resfrs[ul[jj],ii])], psym = 6, col = cols[ii]
     endfor
  endfor
  for ii = 0, nfiles - 1 do $
     oplot, rmasses[ii,*], alog10(rsfrs[ii,*])
     
  stop

  !p.multi =[0,3,nfiles/3+1]
  regions = ['INNFL', 'INTUP', 'INTDN', 'OUTUP', 'OUTDN', 'OPTFL']
  cols = [255, '00ff00'x, '007700'x, 'ffa500'x, 'ff0000'x, '777777'x]
  for ii = 0, nfiles - 1 do begin
     data = mrdfits(files[ii], 1)     
     qui = where(data.REGION eq 'OPTFL')
     tsfrl = data.HA_FLUX[0] * 1d-18 * dluminosity(zs[ii,0], /cm)^2 * 4 *!pi / 1.26d41
     tssfrl = alog10(tsfrl) - data.LMASS[0]
     plot, [data[qui].LMASS[0]], alog10([data[qui].SFR[0]]) - data[qui].LMASS, $
;           xran = minmax(data.LMASS), yran = alog10(minmax(data.SFR[0]/10.^data.LMASS[0])), $
;           xran = minmax(data.LMASS), yran = minmax(tssfrl), $
           xran = [9,11.5], yran = [-12,-8.5], $
           xtitle = 'log !18M!X!D*!N [M!D'+sunsymbol()+'!N]', $
           ytitle = 'log !18sSFR!X [yr!E-1!N]', $
           /nodat
     for jj = 0, n_elements(regions) - 1 do begin
        qui = where(data.REGION eq REGIONS[JJ])
        tdata = data[qui]
        mass = tdata.LMASS[0]
        sfr = tdata.SFR[0]
        ssfr    = alog10(sfr) - mass
        mass_err = 0.5 * (tdata.LMASS[2] - tdata.LMASS[1])
        sfr_err  = 0.5/alog(10) * (tdata.SFR[2] - tdata.SFR[1])/tdata.SFR[0]
        ssfr_err = sqrt(mass_err^2 + sfr_err^2)
        hasfr = tdata.HA_FLUX * 1d-18 * dluminosity(zs[ii,0], /cm)^2 * 4 *!pi / 1.26d41
        hasfr_err = 0.5/alog(10) * (hasfr[2] - hasfr[1]) / hasfr[0]
        hassfr_err = sqrt(hasfr_err^2 + mass_err^2)
        hassfr = alog10(hasfr) - mass        
        oploterror, [mass], [hassfr], mass_err, hassfr_err, $
                    psym = 8, col = cols[jj], errcol = cols[jj]
        oplot, [mass], [ssfr], psym = 8, symsize = 1.3
        oploterror, [mass], [ssfr], mass_err, ssfr_err, $
                    psym = 8, col = cols[jj], errcol = cols[jj]
        legend, /bottom, /right, box = 0, $
                repstr(files[ii], '_lgn_bestfit.fits', ''), $
                charsize = 1
     endfor
  endfor

;  sfrs    = fltarr(n_elements(regions), nfiles, 3)
;  sfrs_Ha = fltarr(n_elements(regions), nfiles, 3)
;  for ii = 0, nfiles - 1 do begin
;
;     data = mrdfits(files[ii], 1)     
;     for jj = 0, n_elements(regions) - 1 do begin
;        tdata = data[where(data.REGION eq regions[jj])]
;        sfrs[jj,ii,*] = tdata.SFR
;        sfrs_Ha[jj,ii,*] = tdata.HA_FLUX * 1d-18 * dluminosity(zs[ii,0], /cm)^2 * 4 *!pi / 1.26d41
;     endfor
;     
;  endfor
;
;  !p.multi = 0
;  
;  plot, alog10(sfrs), alog10(sfrs_Ha), $
;        /iso, /nodat, xran = [-4,4], yran = [-4,4], /xsty, /ysty
;  one_one
;  for jj = 0, n_elements(regions) - 1 do $
;     oploterror, alog10(sfrs[jj,*,0]), alog10(sfrs_ha[jj,*,0]), $
;                 0.5/alog(10) * (sfrs[jj,*,2] - sfrs[jj,*,1])/sfrs[jj,*,0], $
;                 0.5/alog(10) * (sfrs_ha[jj,*,2] - sfrs_ha[jj,*,1])/sfrs_Ha[jj,*,0], $
;                 col = cols[jj], errcol = cols[jj], psym = 8
  
;  oploterror, masses[*,0], alog10(sfrs_ha[*,0]), $
;              0.5 * (masses[*,2] - masses[*,1]), $
;              0.5/alog(10) * (sfrs_ha[*,2] - sfrs_ha[*,1])/sfrs_ha[*,0], $
;              psym = 8, col = '0000ff'x, errcol = '0000ff'x

  
  
  stop
  
end
;plotsfms, 'allLgnBestfits.list'
;plotsfms, 'study3_lgn_Bestfits.list'
;plotsfms, 'study5_lgn_Bestfits.list'
