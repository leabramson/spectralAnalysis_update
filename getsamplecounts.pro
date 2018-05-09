pro getsamplecounts

  readcol, 'allGlassCats.list', cats, $
           f = 'A'
  ncats = n_elements(cats)

  cts = lonarr(ncats)
  scts = lonarr(ncats)
  zcts = lonarr(ncats)
  mcts = lonarr(ncats)
  for ii = 0, ncats - 1 do begin
     data = mrdfits(cats[ii], 1)
     cts[ii] = max(data.NUMBER)
     if strpos(cats[ii], 'ABEL0370') eq -1 AND strpos(cats[ii], 'MACS0416') eq -1 then $
        use = where(data.Z_GLASS ge 1.0 AND data.Z_GLASS le 1.8 $
                    OR $
                    data.KUANG_PZ_A ge 1.0 AND data.KUANG_PZ_A le 1.8, nuse) $
     else if strpos(cats[ii], 'ABEL0370') ge 0 OR strpos(cats[ii], 'MACS0416') ge 0 then $
        use = where(data.Z_GLASS ge 1.0 AND data.Z_GLASS le 1.8 $
                    OR $
                    data.Z_PEAK ge 1.0 AND data.Z_PEAK le 1.8, nuse)
     zcts[ii] = nuse
     scts[ii] = total(data.Z_GLASS gt -99)
     mcts[ii] = total(data[use].MAG_AUTO gt 14 AND data[use].MAG_AUTO le 21.8)
  endfor
  ntot = total(cts)
  ns = total(scts)
  nz = total(zcts)
  nin = total(mcts)
  
  readcol, 'sourceLists/1.0_1.8_14.0_21.8_c99.9_pacut0_ALLSRCS.list', $
           field, id ,ra, dec, z, zq, kz, mag, $
           f = 'A,I,D,D,F,F,F,F', $        
           comment = '#'
  nobj = n_elements(field)

  readcol, '1.0_1.8_14.0_21.8_c99.9_pacut0_ALLSRCS_QUANTITATIVE_3031.list', $
           field, id ,ra, dec, z, zq, kz, mag, $
           f = 'A,I,D,D,F,F,F,F', $        
           comment = '#'
  ngood = n_elements(field)

  print, f = '(%"Total objects in GLASS : %i")', $
         ntot
  print, f = '(%"Total objects with GLASS redshift : %i")', $
         ns
  print, f = '(%"Total objects meeting GLASS or phot z criteria : %i")', $
         nz
  print, f = '(%"Total objects meeting z + mag criteria : %i")', $
         nobj;, nin off by one...?
  print, f = '(%"Total objects meeting quant contam criteria : %i")', $
         ngood
    
end
