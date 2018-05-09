pro plotsfms, lgnlist, explist, $
              NSIG = nsig

  if NOT keyword_set(NSIG) then nsig = 1

  summaryData = mrdfits('objectSummaryData.fits', 1)
  finalMasses = mrdfits('finalMasses.fits', 1)
;  magerr = 2./alog(10) * summaryData.EMU/summaryData.MU

  ;; Get MW values from G13
  samplez = mean(summaryData.Z)
  g13 = mrdfits('mikeGenBasic.fits', 1)
  mws = where(alog10(g13.MSTEL_OBS) ge 10.6 AND alog10(g13.MSTEL_OBS) le 10.8)
  hit = value_locate(g13.REDSHIFT, samplez)
  sf = mws[where(alog10(g13.SFR_OBS[MWS]/g13.MSTEL_OBS[MWS]) ge -11)]
  pas = mws[where(alog10(g13.SFR_OBS[MWS]/g13.MSTEL_OBS[MWS]) lt -11 OR g13.SFR_OBS[MWS] eq 0)]
  mw_sf_m1 = median(alog10(g13.MSTEL_T[hit,sf]))
  mw_sf_s1 = median(alog10(g13.SFR_T[hit,sf])) - 9
  emw_sf_m1 = stddev(g13.MSTEL_T[hit,sf], /nan) / $
              median(g13.MSTEL_T[hit,sf])/alog(10)
  emw_sf_s1 = stddev(g13.SFR_T[hit,sf], /nan) / $
              median(g13.SFR_T[hit,sf])/alog(10)
  mw_pas_m1 = median(alog10(g13.MSTEL_T[hit,pas]))
  mw_pas_s1 = median(alog10(g13.SFR_T[hit,pas])) - 9
  emw_pas_m1 = stddev(g13.MSTEL_T[hit,pas], /nan) / $
               median(g13.MSTEL_T[hit,pas])/alog(10)
  emw_pas_s1 = stddev(g13.SFR_T[hit,pas], /nan) / $
               median(g13.SFR_T[hit,pas])/alog(10)

;  stop
  
  readcol, explist, expfiles, f = 'A'
  readcol, lgnlist, lgnfiles, f = 'A'

  ;; Read evo tracks just to put arrows on the points
  readcol, repstr(lgnlist, 'Bestfit', 'Recon'), lgnrecons, f = 'A'
  
  ngals = n_elements(lgnfiles)

  test     = mrdfits(lgnfiles[0], 1)

  regions  = test.REGION
  nregions = n_elements(regions)
  
  sfrs    = fltarr(ngals, nregions, 2) ;; Lgn, Exp
  esfrs   = fltarr(ngals, nregions, 2)
  lsfrs   = fltarr(ngals, nregions, 2) ;; for Upper limits
  
  masses  = fltarr(ngals, nregions, 2)
  emasses = fltarr(ngals, nregions, 2)

  dms = fltarr(ngals, nregions)
  dss = fltarr(ngals, nregions)
  dnorms = fltarr(ngals, nregions)
  
  files = strarr(ngals)
  names = []
  for ii = 0, ngals - 1 do begin
     for jj = 0, 1 do begin
        case jj of
           0: begin
              file = lgnfiles[ii]
              print, file
              files[ii] = strmid(file, 0, 7)
              case files[ii] of
                 '00900_1': nn = 'SSF'
                 '00451_2': nn = 'CSF'
                 '01916_2': nn = 'PSB'
                 '00660_2': nn = 'PAS'
              endcase
              names = [names, nn]
              qui = where(summaryData.ID eq files[ii])
              mag = summaryData[qui].MU
              mag = mag[0]
              emag = summaryData[qui].EMU
              emag = emag[0]
              print, mag
              evo = mrdfits(lgnrecons[ii], 1, /silent)
              hit = value_locate(evo[0].time, $
                                 getage(summaryData[qui].Z)+[-1,0,1])
              for kk = 0, nregions - 1 do begin
                 thit = where(evo.REGION eq regions[kk])
                 dms[ii,kk] = alog10(evo[thit].MGH[hit[2],1]/evo[thit].MGH[hit[0],1])
                 dss[ii,kk] = alog10(evo[thit].SFH[hit[2],1]/evo[thit].SFH[hit[0],1])
                 dnorms[ii,kk] = sqrt(dms[ii,kk]^2 + dss[ii,kk]^2)
                 print, regions[kk], ' ', evo[thit].REGION, dms[ii,kk], dss[ii,kk]
              endfor
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

  imasses = fltarr(ngals, 2)
  eimasses = fltarr(ngals, 2)
  isfrs = fltarr(ngals,2)
  eisfrs = fltarr(ngals,2)
  for ii = 0, ngals - 1 do begin
     for jj = 0, 1 do begin
        imasses[ii,jj]  = alog10(10.^masses[ii,extract[1],jj] + 2 * 10.^masses[ii,extract[2],jj] + 2 * 10.^masses[ii,extract[3],jj])
        eimasses[ii,jj] = sqrt(emasses[ii,extract[1],jj]^2 + 2 * emasses[ii,extract[2],jj]^2 + 2 * emasses[ii,extract[3],jj]^2 + emag^2)
        isfrs[ii,jj] = sfrs[ii,extract[1],jj] + 2 * sfrs[ii,extract[2],jj] + 2 * sfrs[ii,extract[3],jj]
        eisfrs[ii,jj] = sqrt(esfrs[ii,extract[1],jj]^2 + 2 * esfrs[ii,extract[2],jj]^2 + 2 * esfrs[ii,extract[3],jj]^2 + emag^2)
     endfor
  endfor

  plotsym, 0, /fill
  plot, masses[*,extract[0],0], imasses[*,0], $
        xran = [10.2,11.7], yran = [10.2,11.7], /iso, /nodat
  oploterror, masses[*,extract[0],0], imasses[*,0], $
              emasses[*,extract[0],0], eimasses[*,0], psym = 8
  oploterror, masses[*,extract[0],1], imasses[*,1], $
              emasses[*,extract[0],1], eimasses[*,1], errcol = '777777'x, psym = 8
  one_one
  
;  stop
  
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
     if names[kk] eq 'CSF' or names[kk] eq 'SSF' then $
        oplot, m[q], 0.99 * (m[q]-10.2) + 1.31 - 0.31, linesty = 1, thick = 2 ;; W+14 fixed slope z=1-1.5 low-mass
     polyfill, (mw_pas_m1 + emw_pas_m1 * [-1,-1,1,1]) > !X.CRANGE[0], $
               mw_pas_s1 + emw_pas_s1 * [-1,1,1,-1], $
               /line_fill, col = '00a5ff'x, thick = 1, spacing = 0.1, orien = -60
;     polyfill, (10.7 + 0.1 * [-1,-1,1,1]) > !X.CRANGE[0], $
;               [!Y.CRANGE[0], 0.5, 0.5, !Y.CRANGE[0]], $ ;mw_pas_s1 - emw_pas_s1 * [1,1]
;               /line_fill, col = '00a5ff'x, thick = 1, spacing = 0.1, orien = -60
     if kk eq 1 then  $
        cgtext, mw_pas_m1, mw_pas_s1, /data, $ ; + emw_pas_M1 + 0.15
                'PAS MWs', col = long('0055ff'x), charthick = 4, charsize = 0.9, align = 0.5, orien = 50
     polyfill, (mw_sf_m1 + emw_sf_m1 * [-1,-1,1,1]) > !X.CRANGE[0], $
               mw_sf_s1 + emw_sf_s1 * [-1,1,1,-1], $
               /line_fill, col = 'ffa500'x, thick = 1, spacing = 0.1, orien = 60
     if kk eq 1 then $
        cgtext, mean([mw_sf_m1+emw_sf_m1, !X.CRANGE[0]]), mw_sf_s1 - 0.05, /data, $;
                'SF MWs', col = long('ff5500'x), charthick = 4, charsize = 0.9, align = 0.5
;     oplot, mw_sf_m1 + emw_sf_m1 * [-1,1], replicate(mw_sf_s1 - emw_sf_s1, 2), $
;            col = 'ffa500'x, thick = 10
;     oplot, mw_sf_m1 + emw_sf_m1 * [-1,1], replicate(mw_sf_s1 + emw_sf_s1, 2), $
;            col = 'ffa500'x, thick = 10
;     oplot, replicate(mw_sf_m1 - emw_sf_m1, 2), mw_sf_s1 + emw_sf_s1 * [-1,1], $
;            col = 'ffa500'x, thick = 10
;     oplot, replicate(mw_sf_m1 + emw_sf_m1, 2), mw_sf_s1 + emw_sf_s1 * [-1,1], $
;            col = 'ffa500'x, thick = 10
;     oplot, mw_pas_m1 + emw_pas_m1 * [-1,1], replicate(mw_pas_s1 - emw_pas_s1, 2), $
;            col = '00a5ff'x, thick = 6
;     oplot, mw_pas_m1 + emw_pas_m1 * [-1,1], replicate(mw_pas_s1 + emw_pas_s1, 2), $
;            col = '00a5ff'x, thick = 6
;     oplot, replicate(mw_pas_m1 - emw_pas_m1, 2), mw_pas_s1 + emw_pas_s1 * [-1,1], $
;            col = '00a5ff'x, thick = 6
;     oplot, replicate(mw_pas_m1 + emw_pas_m1, 2), mw_pas_s1 + emw_pas_s1 * [-1,1], $
;            col = '00a5ff'x, thick = 6
  
     oplot, masses[kk,extract[1:*],0], alog10(sfrs[kk,extract[1:*],0] > nsig * lsfrs[kk,extract[1:*],0]), $
            linesty = 0, thick = 6
     oplot, masses[kk,extract[1:*],1], alog10(sfrs[kk,extract[1:*],1] > nsig * lsfrs[kk,extract[1:*],1]), $
            linesty = 2, thick = 6
     for ll = 0, n_elements(extract) - 1 do $
        if sfrs[kk,extract[ll],0] then $
           arrow, masses[kk,extract[ll],0], $
                  alog10(sfrs[kk, extract[ll], 0]), $
                  masses[kk,extract[ll],0] + dms[kk, extract[ll]] / dnorms[kk,extract[ll]] / 4., $
                  alog10(sfrs[kk, extract[ll], 0]) + dss[kk, extract[ll]] / dnorms[kk, extract[ll]] / 4., $
                  col = cols[ll], thick = 4, hsize = 200, /data
     
     q = where(summaryData.ID eq files[kk])
     mu  = summaryData[q].MU
     emu = 1./alog(10) * summaryData[q].EMU / mu
     xctr = 11.
     yctr = 0.2
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
                          psym = 8, col = ec, errcol = ec, symsize = 1.5, $
                          errthick = 4, /nohat
              oplot, [tmasses], [alog10(tsfrs)], psym = 8, symsize = 1.2, col = cols[ii]
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
                names[kk], charsize = 1, charthick = 3;, align = 1
        if kk eq 0 then begin
           legend, pos = [!X.CRANGE[0], !Y.CRANGE[0]+0.6], /data, /left, box = 0, $
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
           x1 = 10.8;!X.CRANGE[1] - 0.1
           y0 = !Y.CRANGE[0] + 0.15
           y1 = !Y.CRANGE[0] + 0.25
           polyfill, [x0,x1,x1,x0], [y0,y0,y1,y1] + 0.25, $
                     /line_fill, orien = 45, thick = 1, spacing = 0.025, col = sfmscol
;           oplot, [x0,x1], replicate(mean([y0,y1]), 2), thick = 4, linesty = 2, col = 'ffa500'x
           cgtext, charsize = 1, charthick = 3, 'W+14 SFMS', $
                   align = 1, x0-0.1, y0 + 0.25, /data
           oplot, [10.8,10.95], [-0.1,-0.1], linesty = 1, thick = 2
           cgtext, charsize = 1, charthick = 3, '1:1 lower-env.', $
                   9.6, -0.15, /data

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
                   align = 1, x0-0.1, y0+0.1, /data
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
;  spawn, 'gv resolvedSSFMS.eps &'
  
;  stop
  
end

;plotsfms, 'study5_lgn_Bestfits.list', 'study5_exp_Bestfits.list', nsig = 2
