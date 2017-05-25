pro plotCovar, inlist, output, $
               METAL = metal

  readcol, inlist, files, f = 'A'
  nfiles = n_elements(files)

  set_plot, 'PS'
  device, filename = output, $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = 10, ysize = 4, /in
  plotsym, 0, thick = 1.5
  !p.multi = [0,3,0]
  if keyword_set(metal) then begin
     for ii = 0, nfiles - 1 do begin
        readcol, files[ii], $
                 sii, ha, oiii, sii2, av, age, metal, ltau, lmass, lssfr, sfr, haf, hbf, oiif, siif, sii2f, /silent
     endif else begin
        readcol, files[ii], $
                 sii, ha, oiii, sii2, av, age, ltau, lmass, lssfr, sfr, haf, hbf, oiif, siif, sii2f, /silent
     endelse

     if ii eq 0 then $
        plot, alog10(age), ltau, psym = 1, $
              xran = [7,10], yran = [7,10], $
              xtitle = 'log age [yr]', ytitle = 'log tau [yr]', /nodat, $
              xthick = 3, ythick = 3, charsize = 1.5, charthick = 4
     case files[ii] of
        'INNFLpost_equal_weights.dat': col = '0000ff'x
        'INTUPpost_equal_weights.dat': col = '00a5ff'x
        'INTDNpost_equal_weights.dat': col = '0055ff'x
        'OUTUPpost_equal_weights.dat': col = 'ffa500'x
        'OUTDNpost_equal_weights.dat': col = 'ff0000'x
     endcase
     oplot, alog10(age), ltau, psym = 8, col = col, symsize = 0.1
     legend, /top, /left, $
             ['INNER', 'INTER_UP', 'INTER_DN', 'OUTER_UP', 'OUTER_DN'], $
             psym = 1, col = ['0000ff'x, '00a5ff'x, '0055ff'x, 'ffa500'x, 'ff0000'x], $
             charsize = 1, charthick = 4
  endfor
  
  for ii = 0, nfiles - 1 do begin
     readcol, files[ii], $
              sii, ha, oiii, sii2, av, age, metal, ltau, lmass, lssfr, sfr, haf, hbf, oiif, siif, sii2f, /silent

     if ii eq 0 then $
        plot, alog10(age), alog10(metal / 0.02), psym = 1, $
;              xran = [9,10], yran = [-0.5,0.5], $
              xtitle = 'log age [yr]', ytitle = 'log Z/Solar', /nodat, $
              xthick = 3, ythick = 3, charsize = 1.5, charthick = 4, /ysty
     case files[ii] of
        'INNFLpost_equal_weights.dat': col = '0000ff'x
        'INTUPpost_equal_weights.dat': col = '00a5ff'x
        'INTDNpost_equal_weights.dat': col = '0055ff'x
        'OUTUPpost_equal_weights.dat': col = 'ffa500'x
        'OUTDNpost_equal_weights.dat': col = 'ff0000'x
     endcase
     oplot, !X.CRANGE, replicate(alog10(0.01/0.02), 2), thick = 6
     oplot, !X.CRANGE, replicate(alog10(0.04/0.02), 2), thick = 6
     oplot, alog10(age), alog10(metal/0.02), psym = 8, col = col, symsize = 0.1
  endfor

  for ii = 0, nfiles - 1 do begin
     readcol, files[ii], $
              sii, ha, oiii, sii2, av, age, metal, ltau, lmass, lssfr, sfr, haf, hbf, oiif, siif, sii2f, /silent

     if ii eq 0 then $
        plot, age / 10.^ltau, alog10(metal / 0.02), psym = 1, $
;              xran = [2,5], yran = [-0.5,0.5], $
              xtitle = '!18e!X-foldings', ytitle = 'log Z/Solar', /nodat, $
              xthick = 3, ythick = 3, charsize = 1.5, charthick = 4, /ysty
     case files[ii] of
        'INNFLpost_equal_weights.dat': col = '0000ff'x
        'INTUPpost_equal_weights.dat': col = '00a5ff'x
        'INTDNpost_equal_weights.dat': col = '0055ff'x
        'OUTUPpost_equal_weights.dat': col = 'ffa500'x
        'OUTDNpost_equal_weights.dat': col = 'ff0000'x
     endcase
     oplot, !X.CRANGE, replicate(alog10(0.01/0.02), 2), thick = 6
     oplot, !X.CRANGE, replicate(alog10(0.04/0.02), 2), thick = 6
     oplot, age/10.^ltau, alog10(metal/0.02), psym = 8, col = col, symsize = 0.1
  endfor
  device, /close
  set_plot, 'X'
   
  !p.multi = 0

end

ctroids = [data.LSF_INNER_PARS[1], $
           data.LSF_INTER_DN_PARS[1], $
           data.LSF_INTER_UP_PARS[1], $
           data.LSF_OUTER_DN_PARS[1], $
           data.LSF_OUTER_UP_PARS[1]]
