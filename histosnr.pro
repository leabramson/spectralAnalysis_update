pro histosnr, inners, inters, outers, output

  readcol, inners, innfiles, f = 'A'
  readcol, inters, intfiles, f = 'A'
  readcol, outers, outfiles, f = 'A'

  innSNR = []
  intSNR = []
  outSNR = []
  for jj = 0, 2 do begin

     case jj of
        0: list = innfiles
        1: list = intfiles
        2: list = outfiles
     endcase
     
     for ii = 0, n_elements(innfiles) - 1 do begin
        readcol, list[ii], l, f, e, m
        q = where(m)
        snr = [f[q]/e[q]]        
        
        case jj of
           0: innSNR = [innSNR, snr]
           1: intSNR = [intSNR, snr]
           2: outSNR = [outSNR, snr]
        endcase

     endfor
     
  endfor

  innH = histogram(innsnr, min = 1, max = 15, bins = 1, loc = bins)
  intH = histogram(intsnr, min = 1, max = 15, bins = 1)
  outH = histogram(outsnr, min = 1, max = 15, bins = 1)

  innh /= total(innh)
  inth /= total(inth)
  outh /= total(outh)

  inh = histofill(bins, innh, /pad)
  ith = histofill(bins, inth, /pad)
  oth = histofill(bins, outh, /pad)

  innc = innsnr[sort(innsnr)]
  intc = intsnr[sort(intsnr)]
  outc = outsnr[sort(outsnr)]

  innx = findgen(n_elements(innsnr)) / (n_elements(innsnr)+1)
  intx = findgen(n_elements(intsnr)) / (n_elements(intsnr)+1)
  outx = findgen(n_elements(outsnr)) / (n_elements(outsnr)+1)
  
  set_plot, 'PS'
  device, filename = output, $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = 5, ysize = 5, /in
  plot, bins, outh, /nodat, $
        xran = [1,15], /xsty, $
        xtitle = '!18S/N!X [pix!E-1!N]', $
        ytitle = '!18f!X!Dpix!N', $
        ysty = 8+1, $
        charsize = 1.25, charthick = 4, xthick = 3, ythick = 3, $
        pos = [0.15,0.15,0.85,0.9]
  axis, yaxis = 1, ytitle = '!18F!X(<!18X!X)', $
        ythick = 3, charsize = 1.25, charthick = 4, yran = [0,1], /ysty
  polyfill, inh.BINS > !X.CRANGE[0], inh.HIST, $
            /line_fill, spacing = 0.05, thick = 1, col = 255
  polyfill, ith.BINS > !X.CRANGE[0], ith.HIST, col = '00a500'x, $
            /line_fill, orien = 45, spacing = 0.05, thick = 1
  polyfill, oth.BINS > !X.CRANGE[0], oth.HIST, col = 'ff5500'x, $
            /line_fill, orien = -45, spacing = 0.05, thick = 1
  oplot, bins, innh, psym = 10, col = 200, thick = 6
  oplot, bins, inth, psym = 10, col = '00a500'x, thick = 6
  oplot, bins, outh, psym = 10, col = 'ff0000'x, thick = 6
  oplot, replicate(median(innsnr), 2), !Y.CRANGE, linesty = 5, thick = 6
  oplot, replicate(median(intsnr), 2), !Y.CRANGE, linesty = 5, thick = 6
  oplot, replicate(median(outsnr), 2), !Y.CRANGE, linesty = 5, thick = 6
  oplot, replicate(median(innsnr), 2), !Y.CRANGE, col = 255, linesty = 5, thick = 3
  oplot, replicate(median(intsnr), 2), !Y.CRANGE, col = '00a500'x, linesty = 5, thick = 3
  oplot, replicate(median(outsnr), 2), !Y.CRANGE, col = 'ff5500'x, linesty = 5, thick = 3
  oplot, innc, innx * !Y.CRANGE[1], thick = 6
  oplot, innc, innx * !Y.CRANGE[1], col = 255, thick = 3
  oplot, intc, intx * !Y.CRANGE[1], thick = 6
  oplot, intc, intx * !Y.CRANGE[1], col = '00a500'x, thick = 3
  oplot, outc, outx * !Y.CRANGE[1], thick = 6
  oplot, outc, outx * !Y.CRANGE[1], col = 'ff5500'x, thick = 3
  plotsym, 8, /fill
  legend, /bottom, /right, /clear, $
          ['Inner', 'Middle', 'Outer'], $
          col = [255, '00a500'x, 'ff5500'x], $
          psym = 8, charsize = 1, charthick = 3, pspacing = 0.5
  
  device, /close
  set_plot, 'X'
  spawn , 'gv '+output

  stop
 
  
;  plot, innsnr[sort(innsnr)], total(innsnr[sort(innsnr)], /cum) / total(innsnr);, /nodat, xran = [1,15], /xsty
;  
;  
;  
  
end

pro doit

  spawn, 'ls MACS1149/00900_1_?_INNFL.masked.pyspec > innlist.list'
  spawn, 'ls MACS0744/00660_2_?_INNFL.masked.pyspec > innlist.list'
  spawn, 'ls MACS2129/00451_2_?_INNFL.masked.pyspec > innlist.list'

  spawn, 'ls MACS1149/00900_1_?_INT??.masked.pyspec > intlist.list'
  spawn, 'ls MACS0744/00660_2_?_INT??.masked.pyspec > intlist.list'
  spawn, 'ls MACS2129/00451_2_?_INT??.masked.pyspec > intlist.list'

  spawn, 'ls MACS1149/00900_1_?_OUT??.masked.pyspec > outlist.list'
  spawn, 'ls MACS0744/00660_2_?_OUT??.masked.pyspec > outlist.list'
  spawn, 'ls MACS2129/00451_2_?_OUT??.masked.pyspec > outlist.list'

  histosnr, 'innlist.list', 'intlist.list', 'outlist.list', $
            'SNRs.eps'
  
end
