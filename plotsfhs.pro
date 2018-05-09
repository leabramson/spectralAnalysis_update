pro plotsfhs

  master = loadsummarydata('..')
  master = master[sort(master.ID)] ;; put it inline w/ the 'ls' call below
  ngals = n_elements(master)
  
  spawn, 'ls *lgn_chains.fits > lgnChains.list'
  spawn, 'ls *exp_chains.fits > expChains.list'

  readcol, 'lgnChains.list', lgnfiles, f = 'A'
  readcol, 'expChains.list', expfiles, f = 'A'

  readcol, 'study5_lgn_Bestfits.list', lgnsum, f = 'A'
  readcol, 'study5_exp_Bestfits.list', expsum, f = 'A'
  
  ;; set clocks
  t0 = getage(10)
  t1 = getage(0)
  dt = 0.05

  time     = findgen((t1-t0)/dt+1) * dt + t0
  redshift = getredshift(time)

  nregions = 4

  lsfhs = fltarr(n_elements(time), ngals, 4)
  elsfhs = fltarr(n_elements(time), ngals, 4)

  ;; add for ref just to see if it works 20180429
  lssfrs = fltarr(n_elements(time) * dt/0.2 + 1, ngals)
  
  esfhs = fltarr(n_elements(time), ngals, 4)
  eesfhs = fltarr(n_elements(time), ngals, 4)
  
  tdexes = lonarr(ngals)
  pdexes = fltarr(ngals)

  psfrs = fltarr(ngals)
  epsfrs = fltarr(ngals)
  pms = fltarr(ngals)

  areas = []
  mmag = mean([1.87,1.67,1.52,1.66]) ;; mean magnification to convert to rough SFR
;  sum = fltarr(n_elements(time), ngals)
  
  ;; Reconstruct the SFHs
  for ii = 0, ngals - 1 do begin

     data = mrdfits(lgnfiles[ii], 1)
     edata = mrdfits(expfiles[ii], 1)
     
     hit = value_locate(redshift, master[ii].Z)
     tdexes[ii] = hit

     lsum = mrdfits(lgnsum[ii], 1)
     esum = mrdfits(expsum[ii], 1)    
     for jj = 0, nregions - 1 do begin
        case jj of
           0: begin
              region = 'OPTFL'
              scl = 1.
           end
           1: begin
              region = 'INNFL'
              scl = 1.
           end
           2: begin
              region = 'INTER'
              scl = 2.
           end
           3: begin
              region = 'OUTER'
              scl = 2.
           end
        endcase
        use = where(data.REGION eq region)
        tdat = data[use]
        tedat = edata[use]
        
        if region eq 'OPTFL' then region = 'INNFL'

        use = where(lsum.REGION eq region)
        tlsum = lsum[use]
        ina = tlsum.AREA_PHYS
        if jj eq 0 then areas = [areas, ina]
        
        use = where(esum.REGION eq region)
        tesum = esum[use]
                
        niter = n_elements(tdat.T0)
        tsfr  = fltarr(n_elements(time), niter)
        if jj eq 1 then tssfrs = fltarr(n_elements(time), niter)
;        print, median(tdat.T0), median(tdat.TAU)
        for kk = 0, niter - 1 do begin
           tsfr[*,kk] = 1./sqrt(2 * !pi * tdat.TAU[kk]^2) $
                        * exp(-0.5 * (alog(time) - tdat.T0[kk])^2 / tdat.TAU[kk]^2) $
                        / time ;; bug caught; lea 20170824
;           sum[*,ii] += scl * tsfr[*,kk] / tsfr[hit,kk] * tdat.SFR[kk] / niter
           tsfr[*,kk] /= tsfr[hit,kk] / tdat.SFR[kk] * tlsum.AREA_PHYS
           if jj eq 1 then $
              tssfrs[*,kk] = tsfr[*,kk] $
                             / total(tsfr[*,kk] * 0.6 * dt * 1d9, /cum)
        endfor

        if jj eq 1 then begin
           case ii of
              0: ssfrs1 = tssfrs
              1: ssfrs2 = tssfrs
              2: ssfrs3 = tssfrs
              3: ssfrs4 = tssfrs
           endcase
        endif
        
;        stop
        
        if jj eq 1 then begin
           pop = max(median(tsfr, dim = 2), foo)           
           psfrs[ii] = median(tssfrs[foo,*], /even)
           epsfrs[ii] = stddev(tssfrs[foo,*], /nan)
;           pms[ii] = total(median(tsfr[0:foo,*] * ina, dim =2)) * 0.6 * dt * 1d9
           pdexes[ii] = foo
        endif
           
        niter = n_elements(tedat.AGE)
        tesfr = fltarr(n_elements(time), niter)
        for kk = 0, niter - 1 do begin
           t_0 = time[hit] - 10.^(tedat.AGE[kk] - 9)
           qui = where(time gt t_0)
           tesfr[qui,kk] = (time[qui] - t_0) / 10.^(tedat.LTAU[kk] - 9)^2 * exp(-(time[qui] - t_0)/10.^(tedat.LTAU[kk] - 9))
           tesfr[qui,kk] /= tesfr[hit,kk] / tedat.SFR[kk] * tesum.AREA_PHYS
        endfor
           
        lsfhs[*,ii,jj]  = median(tsfr, dim = 2)
        elsfhs[*,ii,jj] = stddev(tsfr, dim = 2, /nan) 

        esfhs[*,ii,jj]  = median(tesfr, dim = 2)
        eesfhs[*,ii,jj] = stddev(tesfr, dim = 2, /nan) 

     endfor
  endfor

  print, mean(areas) / mmag

  mmm = mean(areas) / mmag
  
;  stop
  
  cols = ['777777'x, 255, '00a500'x, 'ff5500'x]
;  !P.MULTI = [0,ngals,0]
  xs = 0.1
  ys = 0.25
  xw = (0.975 - 1.5 * xs) / 4
  set_plot, 'PS'
  device, filename = 'logNormalSFHfits.eps', $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = 10.5, ysize = 3.25, /in
  for ii = 0, ngals - 1 do begin
     if ii eq 0 then begin
        plot, time, alog10(lsfhs[*,ii,0]), $ ;/ylog, $
              xran = [0,10], xsty = 8+1, $
              yran = [-2,0.7], $
              xthick = 3, ythick = 3, $
              charthick = 4, charsize = 1.25, /nodat, $
;              xtitle = '!18t!X [Gyr]', $
              ytitle = 'log '+greek('Sigma')+'!D!18SFR!X!N [M!D'+sunsymbol()+'!N yr!E-1!N kpc!E-2!N]', $
              pos = [xs, ys, xs+(ii+1)*xw, 0.85]
     endif else if ii gt 0 AND ii lt ngals - 1 then begin
        plot, time, alog10(lsfhs[*,ii,0]), $ ;/ylog, $
              xran = [1d-6,10], xsty = 8+1, $
              yran = [-2,0.7], $
              xthick = 3, ythick = 3, $
              charthick = 4, charsize = 1.25, /nodat, $
;              xtitle = '!18t!X [Gyr]', $
              ytickname = replicate(' ', 60), $
              pos = [xs+ii*xw, ys, xs+(ii+1)*xw, 0.85], /noer
     endif else if ii eq ngals - 1 then begin
        plot, time, alog10(esfhs[*,ii,0]), $ ;/ylog, $
              xran = [1d-6,10], xsty = 8+1, $
              ysty = 8+1, $
              yran = [-2,0.7], $
              xthick = 3, ythick = 3, $
              charthick = 4, charsize = 1.25, /nodat, $
;              xtitle = '!18t!X [Gyr]', $
              ytickname = replicate(' ', 60), $
              pos = [xs+ii*xw, ys, xs+(ii+1)*xw, 0.85], /noer
        axis, yaxis = 1, yran = !Y.CRANGE + alog10(mmm), /ysty, $
              col = 255, ythick = 3, ytitle = $
              '!18<!Xlog!18 SFR!X!DInner!N!18>!X [M!D'+sunsymbol()+'!N yr!E-1!N]', $
              charthick = 4, charsize = 1.25
     endif

     polyfill, time[tdexes[ii]] + [-1,-1,1,1], $
               [!Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1], !Y.CRANGE[0]], $
               col = 'cccccc'x, /line_fill, thick = 2, spacing = 0.075, orien = 45
     
;     oplot, replicate(time[tdexes[ii]], 2), !Y.CRANGE, $
;            col = '555555'x, thick = 3
     oplot, replicate(time[tdexes[ii]], 2) + 1, !Y.CRANGE, $
            col = '555555'x, thick = 2, linesty = 1
     oplot, replicate(time[tdexes[ii]], 2) - 1, !Y.CRANGE, $
            col = '555555'x, thick = 2, linesty = 1

     for jj = 0, nregions - 1 do begin
        s1 = reform(lsfhs[0:tdexes[ii],ii,jj])
        s2 = reform(lsfhs[tdexes[ii]:*,ii,jj])
        e1 = reform(elsfhs[0:tdexes[ii],ii,jj])
        e2 = reform(elsfhs[tdexes[ii]:*,ii,jj])
        e1[where(~finite(e1))] = 0
        e2[where(~finite(e2))] = 0
        xxx1 = [time[0:tdexes[ii]], reverse(time[0:tdexes[ii]])]
        xxx2 = [time[tdexes[ii]:*] < !X.CRANGE[1], reverse(time[tdexes[ii]:*]) < !X.CRANGE[1]]
        yyy1 = alog10([(s1+e1) < 10.^!Y.CRANGE[1], reverse(s1-e1) > 10.^!Y.CRANGE[0]])
        yyy2 = alog10([(s2+e2) < 10.^!Y.CRANGE[1], reverse(s2-e2) > 10.^!Y.CRANGE[0]])
        polyfill, xxx1, yyy1 > !Y.CRANGE[0], col = cols[jj], $
                  /line_fill, orien = 45, thick = 1, spacing = 0.025
        polyfill, xxx2 < !X.CRANGE[1], yyy2 > !Y.CRANGE[0], col = cols[jj], $
                  /line_fill, orien = -45 + randomn(seed, 1) * 5, thick = 1, spacing = 0.075
        oplot, time[0:tdexes[ii]], alog10(s1), $
               col = cols[jj], thick = 8
        oplot, time[tdexes[ii]:*], alog10(s2), $
               col = cols[jj], thick = 8, linesty = 2
;        oplot, time, sum[*,ii] / ina
     endfor
     axis, xaxis = 1, xtitle = '!18z!X', $
           xtickname = ['4', '3', '2', '1', '0.5'], $
           xtickv = time[value_locate(redshift, [4,3,2,1,0.5])], $
           xticks = 4, $
           charsize = 1.25, charthick = 4
     case master[ii].ID of
        '00900_1': nn = 'SSF'
        '00451_2': nn = 'CSF'
        '01916_2': nn = 'PSB'
        '00660_2': nn = 'PAS'
     endcase
     cgtext, !X.CRANGE[1]-0.25, !Y.CRANGE[1]-0.25, $
             nn, charsize = 0.9, charthick = 4, align = 1
  endfor
  cgtext, (2*xs+ngals*xw)/2., 0.07, /norm, $
          '!18t!X [Gyr]', charsize = 1.25, charthick = 4, align = 0.5
  device, /close
;  spawn, 'gv logNormalSFHfits.eps &'

  set_plot, 'PS'
  device, filename = 'delExponentSFHfits.eps', $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = 10.5, ysize = 3, /in
  for ii = 0, ngals - 1 do begin
     if ii eq 0 then begin
        plot, time, alog10(esfhs[*,ii,0]), $ ;/ylog, $
              xran = [0,10], xsty = 8+1, $
              yran = [-2,0.7], $
              xthick = 3, ythick = 3, $
              charthick = 4, charsize = 1.25, /nodat, $
;              xtitle = '!18t!X [Gyr]', $
              ytitle = 'log '+greek('Sigma')+'!D!18SFR!X!N [M!D'+sunsymbol()+'!N yr!E-1!N kpc!E-2!N]', $
              pos = [xs, ys, xs+(ii+1)*xw, 0.85]
     endif else if ii gt 0 and ii ne ngals - 1 then begin
        plot, time, alog10(esfhs[*,ii,0]), $ ;/ylog, $
              xran = [1d-6,10], xsty = 8+1, $
              yran = [-2,0.7], $
              xthick = 3, ythick = 3, $
              charthick = 4, charsize = 1.25, /nodat, $
;              xtitle = '!18t!X [Gyr]', $
              ytickname = replicate(' ', 60), $
              pos = [xs+ii*xw, ys, xs+(ii+1)*xw, 0.85], /noer
     endif else if ii eq ngals - 1 then begin
        plot, time, alog10(esfhs[*,ii,0]), $ ;/ylog, $
              xran = [1d-6,10], xsty = 8+1, $
              ysty = 8+1, $
              yran = [-2,0.7], $
              xthick = 3, ythick = 3, $
              charthick = 4, charsize = 1.25, /nodat, $
;              xtitle = '!18t!X [Gyr]', $
              ytickname = replicate(' ', 60), $
              pos = [xs+ii*xw, ys, xs+(ii+1)*xw, 0.85], /noer
        axis, yaxis = 1, yran = !Y.CRANGE + alog10(mmm), /ysty, $
              col = 255, ythick = 3, ytitle = $
              '!18<!Xlog!18 SFR!X!DInner!N!18>!X [M!D'+sunsymbol()+'!N yr!E-1!N]', $
              charthick = 4, charsize = 1.25
     endif
     
     polyfill, time[tdexes[ii]] + [-1,-1,1,1], $
               [!Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1], !Y.CRANGE[0]], $
               col = 'cccccc'x, /line_fill, thick = 2, spacing = 0.075, orien = 45
     
;     oplot, replicate(time[tdexes[ii]], 2), !Y.CRANGE, $
;            col = '555555'x, thick = 3
     oplot, replicate(time[tdexes[ii]], 2) + 1, !Y.CRANGE, $
            col = '555555'x, thick = 2, linesty = 1
     oplot, replicate(time[tdexes[ii]], 2) - 1, !Y.CRANGE, $
            col = '555555'x, thick = 2, linesty = 1

     for jj = 0, nregions - 1 do begin
        s1 = reform(esfhs[0:tdexes[ii],ii,jj])
        s2 = reform(esfhs[tdexes[ii]:*,ii,jj])
        e1 = reform(eesfhs[0:tdexes[ii],ii,jj])
        e2 = reform(eesfhs[tdexes[ii]:*,ii,jj])
        e1[where(~finite(e1))] = 0
        e2[where(~finite(e2))] = 0
        xxx1 = [time[0:tdexes[ii]], reverse(time[0:tdexes[ii]])]
        xxx2 = [time[tdexes[ii]:*] < !X.CRANGE[1], reverse(time[tdexes[ii]:*]) < !X.CRANGE[1]]
        yyy1 = alog10([(s1+e1) < 10.^!Y.CRANGE[1], reverse(s1-e1) > 10.^!Y.CRANGE[0]])
        yyy2 = alog10([(s2+e2) < 10.^!Y.CRANGE[1], reverse(s2-e2) > 10.^!Y.CRANGE[0]])
        polyfill, xxx1, yyy1 > !Y.CRANGE[0], col = cols[jj], $
                  /line_fill, orien = 45, thick = 1, spacing = 0.025
        polyfill, xxx2 < !X.CRANGE[1], yyy2 > !Y.CRANGE[0], col = cols[jj], $
                  /line_fill, orien = -45 + randomn(seed, 1) * 5, thick = 1, spacing = 0.075
        oplot, time[0:tdexes[ii]], alog10(s1), $
               col = cols[jj], thick = 8
        oplot, time[tdexes[ii]:*], alog10(s2), $
               col = cols[jj], thick = 8, linesty = 2
     endfor
     axis, xaxis = 1, xtitle = '!18z!X', $
           xtickname = ['4', '3', '2', '1', '0.5'], $
           xtickv = time[value_locate(redshift, [4,3,2,1,0.5])], $
           xticks = 4, $
           charsize = 1.25, charthick = 4
     case master[ii].ID of
        '00900_1': nn = 'SSF'
        '00451_2': nn = 'CSF'
        '01916_2': nn = 'PSB'
        '00660_2': nn = 'PAS'
     endcase
     cgtext, !X.CRANGE[1]-0.25, !Y.CRANGE[1]-0.25, $
             nn, charsize = 0.9, charthick = 4, align = 1
  endfor
  cgtext, (2*xs+ngals*xw)/2., 0.05, /norm, $
          '!18t!X [Gyr]', charsize = 1.25, charthick = 4, align = 0.5
  device, /close
;  spawn, 'gv delExponentSFHfits.eps &'

;  stop
  
  bt = [0.28,0.36,0.37,0.55]
  ncolors = 4
  cgloadct, 3, ncolors = ncolors, /rev
  minbt = 0.2
  maxbt = 0.6
  dbt   = (maxbt - minbt) / ncolors
  tbts  = findgen(ncolors)*dbt + minbt
  cols  = []
  for ii = 0, 3 do begin
     foo = min(abs(tbts - bt[ii]), hit)
     cols = [cols, hit[0]]
  endfor
  cols = [0,1,1,3]
  
  set_plot, 'PS'
  device, filename = 'feedbackTest.eps', $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = 5, ysize = 4, /in
  plotsym, 0, /fill
  plot, time[pdexes], alog10(psfrs), psym = 1, $
        xran = [0,6], $
        yran = [-9,-8], $;[-11,-7], $;
        xtitle = '!18t!X!Dpeak!N [Gyr]', $
        ytitle = '!18sSFR!X!Dinner!N(!18t!X!Dpeak!N) [yr!E-1!N]', $
        xthick = 3, ythick = 3, charsize = 1.25, charthick = 4, $
        xsty = 8+1, xminor = 2, yminor = 2, $
        pos = [0.2,0.15,0.95,0.85]
  axis, xaxis = 1, xran = !X.CRANGE, /xsty, $
        charsize = 1.25, charthick = 4, xtitle =' !18z!X!Dpeak!N', $
        xtickv = time[value_locate(redshift, [5,4,3,2,1])], $
        xtickname = ['5','4','3','2','1'], xticks = 4
;  polyfill, [time, reverse(time)], $
;            [alog10(4./time/1d9), replicate(!Y.CRANGE[0], n_elements(time))] > !Y.CRANGE[0], ;$
;            /line_fill, col = 'ff5500'x, thick = 1, spacing = 0.05, orien = -45
  polyfill, [time, reverse(time)] < !X.CRANGE[1], $
            ([alog10(4./time/1d9)-0.3, alog10(4./reverse(time)/1d9)+0.3] > !Y.CRANGE[0])<!Y.CRANGE[1], $
            /line_fill, col = 'ff5500'x, thick = 1, spacing = 0.05, orien = -45
  oplot, time, alog10(4./time) - 9, col = 'ff0000'x, thick = 8
  for ii = 0, 3 do begin
     foom = total(lsfhs[*,ii,1] * 0.6 * dt * 1d9, /cum)
     foossfr = lsfhs[*,ii,1] / foom
     oplot, time, alog10(foossfr), thick = 4, linesty = 1;, $
;            col = '777777'x
;     oplot, time, alog10(foossfr), thick = 2, $
;            col = cgcolor(string(fix(cols(ii)))), linesty = 4
  endfor
  oploterror, time[pdexes], alog10(psfrs), epsfrs/psfrs/alog(10), $
              /nohat, psym = 8, errthick = 2, symsize = 2
  for ii = 0, 3 do begin
     oplot, [time[pdexes[ii]]], [alog10(psfrs[ii])], psym = 8, $
            col = cgcolor(string(fix(cols[ii]))), symsize = 1.5
     case strmid(master[ii].ID,0,5) of
        '00900': nn = 'SSF'
        '00451': nn = 'CSF'
        '01916': nn = 'PSB'
        '00660': nn = 'PAS'
     endcase
     cgtext, time[pdexes[ii]]+0.2, alog10(psfrs[ii])+0.05, $
             nn, charsize = 0.9, charthick = 3
  endfor
  cgtext, 0.25, -8.57, 'No outflows', $
          charsize = 1, charthick = 3, col = '777777'x, $
          orien = -60, /data
  cgtext, 4.25, -8.675, 'Outflows win', $
          charsize = 1, charthick = 3, col = '777777'x, $
          orien = -28, /data, align = 0.5
  cgtext, 3, -8.8, 'Equilibrium', $
          charsize = 1, charthick = 3, col = '777777'x, $
          orien = -35, /data, align = 0.5
  legend, /left, pos = [!X.CRANGE[0]+0.125,!Y.CRANGE[0]+0.12], box = 1, /clear, $
          '!18<sSFR>!X!DSFMS!N=4!18/t!X!N', $
          linesty = 0, col = 'ff0000'x, thick = 6, pspacing  = 1, $
          charsize = 1, charthick = 4
  cgcolorbar, $;/vert, pos = [!X.WINDOW[1]-0.075,0.7*!Y.WINDOW[1],!X.WINDOW[1]-0.05,$
              ;              !Y.WINDOW[1]-0.05], /norm, $
     pos = [!X.WINDOW[1]-0.3,!Y.WINDOW[1] - 0.06,$
            !X.WINDOW[1]-0.05,!Y.WINDOW[1] - 0.035], $
     charsize = 1, charthick = 4, $
     minrange = minbt, maxrange = maxbt, $
     title = '!18B/T!X!Dobs!N', xthick = 3, /discrete, $
     tickint = 0.1, textthick = 3, ncolors = ncolors
  device, /close
  spawn, 'gv feedbackTest.eps &'
  
  set_plot, 'X'
  
end
