pro plotfittingregion, pyspecin, z_guess

  z_guess = z_guess > 0
  z_guess = z_guess[0]
  
;  readcol, pyspecin, $
;           lambda, f, err, mask, $
;           f = 'F,F,F,F', /silent

  data = mrdfits(pyspecin, 1, /silent)
  lambda = data.LAMBDA
  f      = data.FLUX
  err    = data.STDDEV
  mask   = data.MASK
  
  use = where(mask, compl=masked)

  window, 1, xsize = 700, ysize = 500, retain = 2
  plot, lambda, f, /nodat, $
        xran = [0.9 * min(lambda) > 7000, 1.1 * max(lambda)], $
        yran = [0, 2.5 * median(f[use])], $
        xtitle = greek('lambda')+'!Dobs!N [Ang]', ytitle = '!18f!X', xsty=8+1, $
        pos = [0.15,0.15,0.95,0.90], yminor = 4
  axis, xaxis = 1, xran = [0.9 * min(lambda) > 7000, 1.1 * max(lambda)] / (1+z_guess), $
        xtitle = greek('lambda')+'!S!Drest!N!R!Uguess!N [Ang]', xsty = 1
  xxx = [!X.CRANGE[0], !X.CRANGE[0], lambda[use[0]], lambda[use[0]]]
  yyy = [!Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1], !Y.CRANGE[0]]
  polyfill, xxx, yyy, col = '777777'x, $
            /line_fill, spacing = 0.1, thick = 1, orien = 45
  xxx = [lambda[use[-1]], lambda[use[-1]], !X.CRANGE[1], !X.CRANGE[1]]
  yyy = [!Y.CRANGE[0], !Y.CRANGE[1], !Y.CRANGE[1], !Y.CRANGE[0]]
  polyfill, xxx, yyy, col = '777777'x, $
            /line_fill, spacing = 0.1, thick = 1, orien = -45
  oplot, lambda, f, thick = 1
  oplot, lambda, err, col = 120, thick = 1
  oplot, lambda[use], f[use], thick = 2, col = 'ff5500'x, psym = 1
  oplot, lambda[use], err[use], thick = 2, col = '0055ff'x, psym = 1
  oplot, lambda, mask * median(f[use]), col = '00ff00'x
  oplot, replicate(4000, 2) * (1+z_guess), !Y.CRANGE, col = '555555'x, linesty = 1
  oplot, replicate(4863, 2) * (1+z_guess), !Y.CRANGE, col = '555555'x, linesty = 1
  oplot, replicate(5007, 2) * (1+z_guess), !Y.CRANGE, col = '555555'x, linesty = 1
  oplot, replicate(6563, 2) * (1+z_guess), !Y.CRANGE, col = '555555'x, linesty = 1
  oplot, replicate(6720, 2) * (1+z_guess), !Y.CRANGE, col = '555555'x, linesty = 1
  window, 0

end
