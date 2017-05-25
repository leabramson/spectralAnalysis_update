pro plot1resps, blueRes, redRes, $
                redshift = redshift, $
                photfile = photfile, $
                output   = output, $
                maxL = maxL, $
                minL = minL

  if not keyword_set(REDSHIFT) then begin
     spawn, 'cat '+repstr(blueRes, 'bestspec0', 'analyze')+' | grep "z" > tz.txt'
     readcol, 'tz.txt', redshift, f = 'X,F', /silent
  endif
  redshift = redshift[0]
  
  readcol, blueRes, $
           lambdaB, innTraceB, innEerrB, innPolyB, innModelB, innUsedB
  readcol, redRes, $
           lambdaR, innTraceR, innEerrR, innPolyR, innModelR, innUsedR

  readcol, repstr(blueRes, 'bestspec0', 'bestmodel'), tl, tm
  
  innSpecB = innTraceB; * innPolyB
  innSpecB[where(~innUsedB)] = !values.F_NAN

  innErrB = innEerrB; * innPolyB
  innErrB[where(~innUsedB)] = !values.F_NAN
  
;  innModelB *= innPolyB
  innModelB[where(~innUsedB)] = !values.F_NAN

  innSpecR = innTraceR; * innPolyR
  innSpecR[where(~innUsedR)] = !values.F_NAN

  innErrR = innEerrR; * innPolyR
  innErrR[where(~innUsedR)] = !values.F_NAN
  
;  innModelR *= innPolyR
  innModelR[where(~innUsedR)] = !values.F_NAN

;  lambda = [lambdaB, lambdaR] / (1 + redshift)
;  inner  = [innSpecB, innSpecR]
;  innerErr = [innErrB, innErrR]
;  innMod = [innModelB, innModelR]

  if keyword_set(photfile) then begin
     readcol, photfile, l, d, e, m, $
              f= 'F,F,F,F'
     lmin = min(l)
     lmax = max(l)
  endif else begin
     lmin = min(lambdaB)
     lmax = max(lambdaR)
  endelse

  if keyword_set(maxL) then lmax = maxL * (1+redshift)
  if keyword_set(minL) then lmin = minL * (1+redshift)

  rwaves = findgen(15) * 2000 + 2000 
  
  set_plot, 'PS'
  device, filename = output, $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = 7, ysize = 5, /in
  plotsym, 0, /fill
  plot, lambdaB/(1 + redshift), innSpecB, yran = [0,2. * median(innModelR)], $        
        xran = [0.9 * lmin, 1.1 * lmax]/(1 + redshift), /nodat, $
        xtitle = greek('lambda')+'!Drest!N ['+texToIDL('\AA')+']', $
        ytitle = '!18f!X!D'+greek('lambda')+'!N [10!E-18!N erg s!E-1!N cm!E-2!N '+$
        textoIDl('\AA')+'!E-1!N]', xsty=8+1, $
        pos = [0.15,0.15,0.95,0.9], $
        xthick = 3, ythick = 3, charsize = 1.25, charthick = 4, /ysty
  xr = !X.CRANGE[1] - !X.CRANGE[0]
  xr = findgen(1d3) * xr / (1d3-1) + !X.CRANGE[0]
  ints = value_locate(xr, rwaves/(1+redshift))  
  iuse = where(ints ge 0 and ints lt n_elements(xr) - 1, niuse)
  lab = []
  for ii = 0, niuse - 1 do $
     lab = [lab, string(rwaves[iuse[ii]], f = '(E7.1)')]
  axis, xaxis = 1, $
        xtickv = xr[ints[iuse]], $
        xtickname = lab, $
        xticks = niuse - 1, $
        xtitle = greek('lambda')+'!Dobs!N [Ang]', $
        xthick = 3, charsize = 1.25, charthick = 4, xminor = 4
;  e1 = [innSpecB+innErrB, reverse(innSpecB-innErrB)]
;  e2 = [innSpecR+innErrR, reverse(innspecR-innErrR)]
;  x1 = [lambdaB, reverse(lambdaB)] / (1+redshift)
;  x2 = [lambdaR, reverse(lambdaR)] / (1+redshift)
;  oplot, tl/(1+redshift), tm, thick = 1, col = '555555'x
;  polyfill, x1, e1, col = '777777'x
;  polyfill, x2, e2, col = '777777'x
  oploterror, lambdaB/(1+redshift), innSpecB, innErrB, /nohat, $
              errcol = 'cccccc'x, errthick = 1, psym = 3, linesty = -1, symsize = 0.5
  oploterror, lambdaR/(1+redshift), innSpecR, innErrR, /nohat, $
              errcol = 'cccccc'x, errthick = 1, psym = 3, linesty = -1, symsize = 0.5
  oplot, lambdaB / (1 + redshift), innSpecB, psym = 8, col = '777777'x, symsize = 0.5
  oplot, lambdaR / (1 + redshift), innSpecR, psym = 8, col = '777777'x, symsize = 0.5
  oplot, lambdaB / (1 + redshift), smooth(innSpecB, 5, /nan, /edge_truncate), thick = 4
  oplot, lambdaR / (1 + redshift), smooth(innSpecR, 5, /nan, /edge_truncate), thick = 4
    if keyword_set(photfile) then begin
     oploterror, l/(1+redshift), d/1d-18, e/1d-18, $
                 psym = 8, col = long('ff0000'x), errcol = long('ff0000'x), errthick = 3
     plotsym, 0, thick = 4
     oplot, l/(1+redshift), m/1d-18, col = long('0055ff'x), psym = 8, symsize = 2
     plotsym, 0, /fill
  endif
;  oplot, lambdaB / (1 + redshift), innErrB, col = '777777'x, thick = 1
;  oplot, lambdaR / (1 + redshift), innErrR, col = '777777'x, thick = 1
  oplot, lambdaB / (1 + redshift), innModelB, col = 255, thick = 4
  oplot, lambdaR / (1 + redshift), innModelR, col = 255, thick = 4
  oplot, replicate(4000, 2), !Y.CRANGE, col = '555555'x, linesty = 1, thick = 2
  oplot, replicate(4863, 2), !Y.CRANGE, col = '555555'x, linesty = 1, thick = 2
  oplot, replicate(5007, 2), !Y.CRANGE, col = '555555'x, linesty = 1, thick = 2
  oplot, replicate(6563, 2), !Y.CRANGE, col = '555555'x, linesty = 1, thick = 2
  oplot, replicate(6720, 2), !Y.CRANGE, col = '555555'x, linesty = 1, thick = 2
  legend, /top, /left, box = 0, $
          [strcompress('!18z!X=!N'+string(redshift, f = '(F6.3)'), /rem), $
           'spectrum', 'smoothed', 'photometry', 'model spec.', 'model photom'], $
          psym = [0,1,0,8,0,8], linesty = [-1,0,0,0,0,0], $
          col = [0,'777777'x, 0, 0, '0000ff'x, '0055ff'x], $
          charsize = 1.1, charthick = 3, symsize = [0,0.5,0,1,0,2], pspacing = 1
  device, /close
  set_plot, 'X'
  
;  stop
  
end
