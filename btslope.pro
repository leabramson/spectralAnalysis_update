pro btslope, evolist, bestfitlist, $
             DT = dt;, $
;             DOZ = doz

  if not keyword_set(DT) then dt = 1. ; tobs +/- this in Gyr

  summaryData = mrdfits('objectSummaryData.fits', 1)
  finalMasses = mrdfits('finalMasses.fits', 1)

  readcol, evolist, files, f= 'A'
  readcol, bestfitlist, sumfiles, f= 'A'
  ngals = n_elements(files)

  test = mrdfits(files[0], 1)

  regions = test.REGION
  time    = test[0].TIME
  z       = test[0].REDSHIFT

  masses = fltarr(n_elements(time), ngals, 2)
  sfrs   = fltarr(n_elements(time), ngals, 2)
  sizes  = fltarr(n_elements(time), ngals, 2)
  bts    = fltarr(n_elements(time), ngals, 2)  

  robs = fltarr(ngals, 2)
  btobs = fltarr(ngals, 2)
  
  tobs  = fltarr(ngals)
  ssfrs = fltarr(ngals,2)
  
  for ii = 0, ngals - 1 do begin

     id = strmid(files[ii], 0, 7)
     
     data  = mrdfits(files[ii], 1)
     q = where(data.REGION eq 'OPTFL')

     masses[*,ii,*] = data[q].OMASS
     sfrs[*,ii,*]   = data[q].OSFR
     sizes[*,ii,*]  = data[q].HM_RAD
     bts[*,ii,*]    = data[q].BT

     tobs[ii]  = data[q].TOBS

     robs[ii,*] = data[q].HM_RAD_OBS
     btobs[ii,*] = data[q].BT_OBS
     
     tdata = mrdfits(sumfiles[ii], 1)
     q = where(tdata.REGION eq 'OPTFL')
     ssfrs[ii,0] = alog10(tdata[q].SFR[0]) - tdata[q].LMASS[0]
     q2 = where(finalMasses.ID eq id)
     ssfrs[ii,1] = sqrt(finalMasses[q2].ELMASS^2 + $
                        (0.5/alog(10)*(tdata[q].SFR[2] - tdata[q].SFR[1])/tdata[q].SFR[0])^2)
     
  endfor


  ncolors = 12
  cgloadct, 33, ncolors = ncolors, /rev, $;, $
            clip = [10,240];, /brewer
  minssfr = -11
  maxssfr = -8
  dssfr   = (maxssfr - minssfr) / 8.
  tssfrs  = findgen(9)*dssfr + minssfr
  cols  = []
  for ii = 0, ngals - 1 do begin
     foo = min(abs(tssfrs - ssfrs[ii,0]), hit)
     cols = [cols, hit[0]]
  endfor

  sss = reverse(sort(ssfrs[*,0]))
  
  set_plot, 'PS'
  xsize = 8.5
  ysize = 4
  ar = ysize / xsize
  device, filename = 'sizeMassEvo.eps', $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = xsize, ysize = ysize, /in
  plot, masses[*,0,0], alog10(sizes[*,0,0]), /nodat, $
        psym = 8, $
        xran = [9.4,11.4], /xsty, $
        yran = [0,1], $
        xthick = 3, ythick = 3, charsize = 1.25, charthick = 4, $
;        xtitle = 'log !18M!X!D*!N [M!D'+sunsymbol()+'!N]', $
        ytitle = 'log !18r!X!Dhalf-mass!N [kpc]', $
        pos = [0.1,0.15,0.5,0.9]
  x = findgen((!X.CRANGE[1] - !X.CRANGE[0])/0.01+1) * 0.01 + !X.CRANGE[0]
  xxx = [x, reverse(x)]
  yyyred = [0.22 + 0.76 * (x - 10.7) - 0.2, $
            reverse(0.22 + 0.76 * (x - 10.7) + 0.2)]
  yyyblue= [0.70 + 0.22 * (x - 10.7) - 0.17, $
            reverse(0.70 + 0.22 * (x - 10.7) + 0.17)]
  polyfill, xxx, yyyblue < !Y.CRANGE[1], col = 'ffc000'x;, $
;            /line_fill, thick = 1, spacing = 0.05, orien = 45
  polyfill, xxx, (yyyred > !Y.CRANGE[0]) < !Y.CRANGE[1], col = '00c0ff'x;, $
;            /line_fill, thick = 1, spacing = 0.05, orien = -45

  plotsym, 0, /fill
;  pcols = [255,'00a500'x,'00bbff'x,'0000ff'x]
  for ii = 0, ngals - 1 do begin

     jj = sss[ii]
     
     ;toplot = where(time gt tobs[jj] - dt $
     ;               AND $
     ;               time lt tobs[jj] + dt)

     toplot = value_locate(z, [2,1.75,1.5,1.25,getredshift(tobs[jj]),0.75,0.5])
     
     tsizes   = reform(sizes[toplot,jj,0]);, 3, /edge_truncate)
     tesizes  = reform(sizes[toplot,jj,1] / sizes[toplot,jj,0]) / alog(10);, 3, /edge_truncate)
     
     tmasses  = reform(masses[toplot,jj,0]);, 3, /edge_truncate)
     temasses = reform(masses[toplot,jj,1]);, 3, /edge_truncate)

     xxx = [tmasses, reverse(tmasses)] > !X.CRANGE[0]
     yyy = [alog10(tsizes)-tesizes, reverse(alog10(tsizes)+tesizes)]

     polyfill, xxx, yyy, col = '555555'x, $; cgcolor(string(fix(jj)))
               /line_fill, orien = 60, $
               thick = 1, spacing = 0.05
;     oplot, tmasses, alog10(tsizes)-tesizes, thic = 4, col = '777777'x;, col = cgcolor(string(cols[jj])) 
;     oplot, tmasses, alog10(tsizes)+tesizes, thic = 4, col = '777777'x;, col = cgcolor(string(cols[jj])) 
     oplot, tmasses, alog10(tsizes), thick = 4, col = '222222'x

     epochs = value_locate(time, tobs[jj]+[-1,1,0])
     psyms = [2,4,8]
     
;     oplot, masses[epochs,jj,0], alog10(sizes[epochs, jj, 0]), $
;                 masses[epochs,jj,1], sizes[epochs,jj,1] / sizes[epochs,jj,0] / alog(10), $
;                 psym = 8, symsize = 1.25;, errthick = 4, /nohat, 

     tsfrs = sfrs[epochs,jj,0] - masses[epochs,jj, 0]
     for kk = 0, n_elements(epochs) - 1 do begin
        foo = min(abs(tssfrs - tsfrs[kk]), hit)
        oplot, [masses[epochs[kk],jj,0]], [alog10(sizes[epochs[kk], jj, 0])], $
               psym = psyms[kk], symsize = 1.25, thick = 4
        oplot, [masses[epochs[kk],jj,0]], [alog10(sizes[epochs[kk], jj, 0])], $
               psym = psyms[kk], symsize = 1, col = cgcolor(string(hit[0])), thick = 2
     endfor
  endfor
  cgcolorbar, pos = [0.125,0.825,0.25,0.85], $
              charsize = 1, charthick = 3, $
              minrange = minssfr, maxrange = maxssfr, $
              title = 'log !18sSFR!X', thick = 3, /discrete, $
              xtickint = 1, textthick = 3, ncolors = ncolors;, $
;              oob_low = 240, oob_hi = 10
  plotsym, 8, /fill
  legend, /bottom, /left, /clear, $
          ['vdW+14 SF', $
           'vdW+14 Pas'], $
          charsize = 1, charthick = 4, $
          col = ['ffc000'x, '00c0ff'x], psym = 8, $
          pspacing = 0.5, spacing = 1.2, symsize = [1, 1]+0.2
  plotsym, 0, /fill
;  legend, pos = [11,0.1], /data, box = 0, $;/horiz, $
;          reverse(['1D (used)', '2D (SE)']), $
;          psym = reverse([8, 1]), symsize = [1.25,1.25], $;col = '777777'x, $
;          charsize = 1, charthick = 3, pspacing = 0, spacing = 1.35
  
  plot, masses[*,0,0], bts[*,0,0], /nodat, $
        psym = 8, $
        xran = !X.CRANGE + [1d-6,0], /xsty, $
        yran = [(min(btobs)-0.2) > 0,max(btobs)+0.3],  ytickint = 0.1, $
        ysty = 8+1, $
        ytickname = replicate(' ',60), $
        xthick = 3, ythick = 3, $
        charsize = 1.25, charthick = 4, $
;        xtitle = 'log !18M!X!D*!N [M!D'+sunsymbol()+'!N]', $
        pos = [0.5,0.15,0.9,0.9], /noer, yminor = 2
  axis, yaxis = 1, ytitle = '!18B/T!X', $
        ythick = 3, charsize = 1.25, charthick = 4, $
        yminor = 2, ytickint = 0.1, yran = !Y.CRANGE, /ysty

  ;; Read your SDSS shit
  local = mrdfits('bt_stuff_for_dan.fits', 1)
  lm = local.LOC + 0.15
  lllx = ([lm, reverse(lm)] > !X.CRANGE[0]) < !X.CRANGE[1]
  llly = ([local.P25, reverse(local.P75)] > !Y.CRANGE[0]) < !Y.CRANGE[1]
  polyfill, lllx, llly, col = cgcolor('pink');, /line_fill, $
;            thick = 1, orien = -45, spacing = 0.05

  for ii = 0, ngals- 1 do begin

     jj = sss[ii]
     
;     toplot = where(time gt getage(2) $
;                    AND $
;                    time lt getage(0.5))

     toplot = value_locate(z, [2,1.75,1.5,1.25,getredshift(tobs[jj]),0.75,0.50]);getredshift(value_locate(masses[*,jj,0], !X.CRANGE[0])),
     
     tsizes   = reform(bts[toplot,jj,0]);, 3, /edge_truncate)
     tesizes  = reform(bts[toplot,jj,1]);, 3, /edge_truncate)

     tmasses  = reform(masses[toplot,jj,0]) > !X.CRANGE[0];, 3, /edge_truncate)
     temasses = reform(masses[toplot,jj,1]);, 3, /edge_truncate)

     xxx = [tmasses, reverse(tmasses)] > !X.CRANGE[0]
     yyy = ([tsizes-tesizes, reverse(tsizes+tesizes)] > !Y.CRANGE[0]) < !Y.CRANGE[1]

     polyfill, xxx, yyy, col = '555555'x, $;cgcolor(string(cols[jj])), $
               /line_fill, orien = 60, $
               thick = 1, spacing = 0.05
;     oplot, tmasses, tsizes-tesizes, thic = 4, col = cgcolor(string(cols[jj])) 
;     oplot, tmasses, tsizes+tesizes, thic = 4, col = cgcolor(string(cols[jj])) 
     oplot, tmasses, tsizes, thick = 4, col = '222222'x

     epochs = value_locate(time, tobs[jj]+[-1,1,0])
     
;     oplot, masses[epochs,jj,0], bts[epochs,jj,0], $
;                 masses[epochs,jj,1], sizes[epochs,jj,1] / sizes[epochs,jj,0] / alog(10), $
;                 psym = 8, symsize = 1.25;, errthick = 4, /nohat, 
                 
     tsfrs = sfrs[epochs,jj,0] - masses[epochs,jj, 0]
     for kk = 0, n_elements(epochs) - 1 do begin
        foo = min(abs(tssfrs - tsfrs[kk]), hit)
        oplot, [masses[epochs[kk],jj,0]], [bts[epochs[kk], jj, 0]], $
               psym = psyms[kk], symsize = 1.25, thick = 4
        oplot, [masses[epochs[kk],jj,0]], [bts[epochs[kk], jj, 0]], $
               psym = psyms[kk], symsize = 1, col = cgcolor(string(hit[0])), thick = 2
;        oplot, [masses[epochs[kk],jj,0]], [bts[epochs[kk], jj, 0]], $
;               psym = 8, symsize = 1, col = cgcolor(string(hit[0]))
     endfor
     
  endfor
  legend, /bottom, /right, /clear, $
          ['!18t!X=!18t!X!Dobs!N-1 Gyr', $
           '!18t!X=!18t!X!Dobs!N+1 Gyr', $
           '!18t!X=!18t!X!Dobs!N'], $
          charsize = 1, charthick = 4, $
          col = 0, psym = psyms, $
          pspacing = 0.5, spacing = 1.2, symsize = [1, 1, 1]+0.2

  cgtext, mean([0.1,0.9]), 0.03, 'log !18M!X!D*!N [M!D'+sunsymbol()+'!N]', $
          charsize = 1.25, charthick = 4, align = 0.5, /norm

  device, /close
  spawn, 'gv sizeMassEvo.eps &'

  stop
  
end
;btslope, 'study5_lgn_Structs.list', 'study5_lgn_Bestfits.list', dt = 1.5
