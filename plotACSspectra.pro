pro plotACSspectra, infits, output, $
                    STRETCH = stretch, $
                    CROP = crop, $
                    NEG = neg

  if NOT keyword_set(STRETCH) then stretch = 7
;  if NOT keyword_set(HEIGHT) then height = 2. ;; in asec
  if NOT keyword_set(CROP) then crop = 10 ;; in pix
  if NOT keyword_set(NEG) then neg = 0 else neg = 1

  data = mrdfits(infits, 1, /silent)
  di   = data.DIRECT_IM
  s2d  = data.SPEC2D
  msk  = data.MASK
  lambda = data.LAMBDA
  
  s = size(di, /dim)
  
  foo = mrdfits(data.FILENAME, 'sci', head, /silent)
  pscale = sxpar(head, 'CD2_2') ;; It's downsampled by ~2x from WFC3 plate scale
  lscale = sxpar(head, 'CD1_1')
  pivot  = sxpar(head, 'CRVAL1')
;  tpixz  = lscale / pivot * 10 * (1+zpick) ;; 10 pix offset in z space
;  pscale = 0.13 ;; WFC3IR "/pix    
  run = (findgen(s[0]) - (s[0]+1)/2) * pscale
  run += pscale/2

  innOff = median(data.INNERUP - data.TRACE)
  intOff = median(data.INTERUP - data.TRACE)
  outOff = median(data.OUTERUP - data.TRACE)

;  kpc       = 1. / zang(1.0, zpick, /silent) ;optFlRes.ZFIT) ;; kpc / asec
;  kpcPix    = kpc * pscale                  ;; kpc / pix
;  hlrKpc    = data.RE * kpcPix            ;; re in Kpc
  hlrObs    = data.RE * pscale            ;; re in asec
;  normRadii = [0, mean([innOff,intOff]),mean([intOff,outOff])] / data.RE
  normRadii = [0, mean([data.EXTR_INTDN_OFF, data.EXTR_INTUP_OFF]), $
               mean([data.EXTR_OUTDN_OFF, data.EXTR_OUTUP_OFF])] / data.RE
  
  bdots    = [-reverse(normRadii[-2:-1]), 0, normRadii[1:2]] * data.RE
  bctroids = [data.LSF_OUTER_DN_PARS[1], $
              data.LSF_INTER_DN_PARS[1], $
              data.LSF_INNER_PARS[1], $
              data.LSF_INTER_UP_PARS[1], $
              data.LSF_OUTER_UP_PARS[1]]
  xs1 = crop
  xs2 = n_elements(run) - crop

  regions  = ['OPTFL', $
              'INNFL', $
              'INTUP', $
              'INTDN', $
              'OUTUP', $
              'OUTDN']
  nregions = n_elements(regions)
  
  set_plot, 'PS'
  device, filename = output, $
          /col, /encap, /decomp, bits_per_pixel = 8, $
          xsize = 7, ysize = 5, /in    

  scols = [255, '00aa00'x, '005500'x, 'ffa500'x, 'ff5500'x]
  plot, run[xs1:xs2], run[xs1:xs2], /xsty, /ysty, /nodat, $
        xtitle = greek('delta')+'!18x!X ["]', $
        ytitle = greek('delta')+'!18y!X ["]', /iso, $
        pos = [0.1, 0.7, 0.25, 0.90], $
        xthick = 4, ythick = 4, charthick = 3, charsize = 1, $
        xtickint = 1, ytickint = 1, title = 'direct image'
  cgimage, di[xs1:xs2,xs1:xs2], stretch = stretch, /over, /neg
  x = findgen(101) * (2*!pi) / 100.
  oplot, hlrobs * cos(x), hlrobs * sin(x), col = '00a5ff'x
  oplot, data.RE_SEX * cos(x) * pscale, data.RE_SEX * sin(x) * pscale, col = '00a5ff'x, linesty = 1
  oplot, run, replicate(innOff, n_elements(run))  * pscale, col = scols[0], linesty = 0, thick = 1
  oplot, run, replicate(-innOff, n_elements(run)) * pscale, col = scols[0], linesty = 0, thick = 1
  oplot, run, replicate(intOff, n_elements(run))  * pscale, col = scols[1], linesty = 0, thick = 1
  oplot, run, replicate(-intOff, n_elements(run)) * pscale, col = scols[2], linesty = 0, thick = 1
  oplot, run, replicate(outOff, n_elements(run))  * pscale, col = scols[3], linesty = 0, thick = 1
  oplot, run, replicate(-outOff, n_elements(run)) * pscale, col = scols[4], linesty = 0, thick = 1
  oplot, (bctroids - bctroids[2]) * pscale, bdots * pscale, psym = 1, col = '0055ff'x, $
         symsize = 0.5
  plot, run[xs1:xs2], run[xs1:xs2], /nodat, $
        pos = [!X.WINDOW[0], !Y.WINDOW[0], !X.WINDOW[1], !Y.WINDOW[1]], $
        xtickname = replicate(' ', 60), ytickname = replicate(' ',60), $
        xthick = 4, ythick = 4, charthick = 3, charsize = 1, $
        xtickint = 1, ytickint = 1, /noer, /iso, yminor = 2, $
        xminor = 2, yticklen = 0.075, xticklen = 0.075
  
  st = !X.WINDOW[1]+0.05
  plot, lambda, data.TRACE, /nodat, $
        xtitle = greek('lambda')+'!Dobs!N ['+texToIDL('\AA')+']', $
        xthick = 5, ythick = 5, charthick = 3, charsize = 1, $
        pos = [st,!Y.WINDOW[0], 0.90,!Y.WINDOW[1]], /noer, $
        ytickname = replicate(' ',60), yran = [xs1,xs2], title = infits, /xsty
;  spec2d = s2d[*,xs1:xs2] * (1 - msk[*,xs1:xs2])
  spec2d = data.SPECORIG[*,xs1:xs2]; * (1 - msk[*,xs1:xs2])
  if neg then begin
     cgimage, spec2d, stretch = stretch, /over, /neg
     cgtext, lambda[0] + 200., 0.9 * xs2, /data, align = 0, $
             'G800L', charthick = 4, charsize = 1.2
  endif else begin
     cgimage, spec2d, stretch = stretch, /over
     cgtext, lambda[0] + 200., 0.9 * xs2, /data, align = 0, $
             'G800L', charthick = 4, col = 'ffffff'x, charsize = 1.2
  endelse
  scols = [255, '00aa00'x, '005500'x, 'ffa500'x, 'ff5500'x]
  oplot, lambda, data.innerUP, col = scols[0]      , linesty = 0, thick = 1
  oplot, lambda, data.innerDN, col = scols[0]      , linesty = 0, thick = 1
  oplot, lambda, data.interUp, col = scols[1], linesty = 0, thick = 1
  oplot, lambda, data.interDN, col = scols[2], linesty = 0, thick = 1
  oplot, lambda, data.outerUP, col = scols[3], linesty = 0, thick = 1
  oplot, lambda, data.outerDN, col = scols[4], linesty = 0, thick = 1

  norm = data.F_INNER / data.SENSITIVITY
  use  = where(lambda ge 6000 AND lambda le 9000)
  height = max(norm[use])
  
  plot, lambda, data.F_OPTI / data.SENSITIVITY, $
        xran = [6000,9000], yran = [0, 1.5 * height], xsty = 1, $ ;;1.1 * max(optFlRes.SPEC * norms[0])
        xtitle = greek('lambda')+'!Dobs!N ['+texToIDL('\AA')+']', $
        ytitle = '!18f!X!D'+greek('lambda')+'!N + offset', /nodat, pos = [0.1,0.1,0.90,0.6], /noer, $
        xthick = 4, ythick = 4, charthick = 3, charsize = 1, yminor = 2, /ysty
;  axis, xaxis = 1, xran = [8000,18000], xtitle = greek('lambda')+'!Dobs!N', $
;        xthick = 4, ythick = 4, charthick = 3, charsize = 1
  oplot, replicate(4000, 2), !Y.CRANGE, col = '555555'x, linesty = 1
  oplot, replicate(4863, 2), !Y.CRANGE, col = '555555'x, linesty = 1
  oplot, replicate(5007, 2), !Y.CRANGE, col = '555555'x, linesty = 1
  oplot, replicate(6563, 2), !Y.CRANGE, col = '555555'x, linesty = 1
  oplot, replicate(6720, 2), !Y.CRANGE, col = '555555'x, linesty = 1
  scols = ['777777'x, 255, '00aa00'x, '005500'x, 'ffa500'x, 'ff5500'x]
;  mcols = [0,         120, '005500'x, '002200'x, 'ff5500'x, 'ff0000'x]
  offs  = [6, 5, 4, 3, 2, 1]
  for ii = 0, nregions - 1 do begin
     case ii of
        0: begin
           spec = data.F_OPTI     / data.SENSITIVITY
           spec /= max(spec[use]) / (1.2 * height)
        end
        1: spec = data.F_INNER    / data.SENSITIVITY  
        2: spec = data.F_INTER_UP / data.SENSITIVITY  
        3: spec = data.F_INTER_DN / data.SENSITIVITY  
        4: spec = data.F_OUTER_UP / data.SENSITIVITY  
        5: spec = data.F_OUTER_DN / data.SENSITIVITY  
     endcase
;     lambda = results.LAMBDA
;     norm = mean(spec[where(lambda ge 8000 AND lambda le 9000)], /nan)
;     oplot, lambda, smooth(spec / norm * offs[ii], 2, /nan), col = scols[ii], thick = 3
     oplot, lambda, spec, col = scols[ii], thick = 3
;     oplot, lambda, smooth(results.MODEL / norm * offs[ii], 2, /nan), col = mcols[ii], thick = 4
;     endif else begin
;        oplot, lambda, smooth(results.SPEC, 2, /nan), col = scols[ii], thick = 3
;        oplot, lambda, smooth(results.MODEL, 2, /nan), col = mcols[ii], thick = 4
 ;    endelse
  endfor
  io = innOff / data.RE
  to = intOff / data.RE
  oo = outOff / data.RE
  legend, /top, /right, box = 0, $
          ['!18R/R!X!De!N!X<'+string(io,f='(F4.2)'), 'fit', $
           string(io,f='(F4.2)')+'<!18R/R!X!De!N!X<'+string(to,f='(F4.2)'), 'fit', $
           string(to,f='(F4.2)')+'<!18R/R!X!De!N<'+string(oo,f='(F4.2)'), 'fit'], $
          col = [255,100,'00a500'x,'004400'x,'ffa500'x,'ff0000'x], $
          pspacing = 1, thick = [2,4,2,4,2,4], $
          charsize = 0.7, charthick = 2, linesty = replicate(0,6)
  legend, /bottom, /right, box = 0, $ ; pos = [0.17,8.4], /data
          ['m!DF814W!N='+string(data.MAG, f = '(F5.2)'), $
          '!18R!X!De!N='+string(hlrObs, f = '(F4.2)')+'"'], $
          charsize = 1, charthick = 3
  
  
;  stop
  
  device, /close
  set_plot, 'X'
  spawn, 'gv '+output+' &'
  
end
