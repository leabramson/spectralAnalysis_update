pro plot2pa, pa1, pa2, $
             z, $
             domasked = domasked

  if NOT keyword_set(domasked) then msksfx = '' else msksfx = '.masked'
  if N_ELEMENTS(z) eq 0 then z = 0
  
  regions = ['INNFL', 'INTUP', 'INTDN', 'OUTUP', 'OUTDN', 'OPTFL']+msksfx
  nregions = n_elements(regions)

  dl = 2.
  for ii = 0, nregions - 1 do begin
     d1 = mrdfits(pa1+'_'+regions[ii]+'.unified.dat', 1)
     d2 = mrdfits(pa2+'_'+regions[ii]+'.unified.dat', 1)

     l = d1.lambda
     l = findgen((max(l) - min(l)) / dl + 1) * dl + min(l)

     f1 = interpol(d1.FLUX, d1.LAMBDA, l)
     e1 = interpol(d1.STDDEV, d1.LAMBDA, l)

     f2 = interpol(d2.FLUX, d2.LAMBDA, l)
     e2 = interpol(d2.STDDEV, d2.LAMBDA, l)

     spec = (f1/e1^2 > 0 + f2/e2^2 > 0) / (1/e1^2 + 1/e2^2)

     case ii of
        0: begin
           inner = spec
           err = sqrt(1./ (1/e1^2 + 1/e2^2))
        end
        1: interUP = spec
        2: interDN = spec
        3: outerUP = spec
        4: outerDN = spec
        5: opt = spec / 8.
     endcase
  endfor
  
  set_plot, 'PS'
  device, filename = strmid(pa1, '_1')+'_combinedPAextractions.eps', $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = 7, ysize = 5, /in
  plot, l / (1+z), inner, /nodat, $
        xtitle = greek('lambda')+'!Drest!N ['+texTOIdl('\AA')+']', $
        ytitle = '!18f!X!D', $
        yran = [0,max(inner[where(l gt 10000 AND l lt 16000)]) * 1.1], $
        charsize = 1.25, charthick = 3, xthick = 4, ythick = 4
  scols = reverse(['555555'x, 255, '00aa00'x, '005500'x, 'ffa500'x, 'ff5500'x])
  for ii = 0, nregions - 1 do begin
     case ii of
        5: spec = opt
        4: spec = inner  
        3: spec = interUP
        2: spec = interDN
        1: spec = outerUP
        0: spec = outerDN
     endcase
     use = where(spec / err gt 0)
     oplot, l[use] / (1+z), spec[use], col = scols[ii], thick = 4
  endfor
  oplot, l / (1+z), err, thick = 2, col = '777777'x
  oplot, replicate(4000, 2), !Y.CRANGE, col = '555555'x, linesty = 1
  oplot, replicate(4863, 2), !Y.CRANGE, col = '555555'x, linesty = 1
  oplot, replicate(5007, 2), !Y.CRANGE, col = '555555'x, linesty = 1
  oplot, replicate(6563, 2), !Y.CRANGE, col = '555555'x, linesty = 1
  oplot, replicate(6720, 2), !Y.CRANGE, col = '555555'x, linesty = 1

  device, /close
  set_plot, 'X'
  spawn, 'gv '+strmid(pa1, '_1')+'_combinedPAextractions.eps &'
  
end
