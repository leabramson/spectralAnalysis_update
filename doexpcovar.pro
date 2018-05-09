pro doexpcovar, dir, output, $
                FOLDED = folded

  if NOT keyword_set(FOLDED) then folded = 0 else folded = 1

  optf  = dir+'/OPTFL.maskedpost_equal_weights.dat'
  innf  = dir+'/INNFL.maskedpost_equal_weights.dat'
  if NOT folded then begin
     intfu = dir+'/INTUP.maskedpost_equal_weights.dat'
     intfd = dir+'/INTDN.maskedpost_equal_weights.dat'
     outfu = dir+'/OUTUP.maskedpost_equal_weights.dat'
     outfd = dir+'/OUTDN.maskedpost_equal_weights.dat'
     files = [optf, innf, intfu, intfd, outfu, outfd]
  endif else begin
     intf  = dir+'/INTER.maskedpost_equal_weights.dat'
     outf  = dir+'/OUTER.maskedpost_equal_weights.dat'
     files = [optf, innf, intf, outf]
  endelse
   
  maxl = []
  for ii = 0, n_elements(files) - 1 do $
     maxl = [maxl, file_lines(files[ii])]
  len = max(maxl)
   
  chains = {REGION:    string(0), $
            Z:         replicate(!VALUES.F_NaN, len),  $
            SII_EW1:   replicate(!VALUES.F_NaN, len),  $
            SII_EW2:   replicate(!VALUES.F_NaN, len),  $
            H_EW:      replicate(!VALUES.F_NaN, len),  $
            OIII_EW:   replicate(!VALUES.F_NaN, len),  $
            AGE:       replicate(!VALUES.F_NaN, len),  $
            LTAU:      replicate(!VALUES.F_NaN, len),  $
            AV:        replicate(!VALUES.F_NaN, len),  $
            LMASS:     replicate(!VALUES.F_NaN, len),  $
            SFR:       replicate(!VALUES.F_NaN, len),  $
            LSSFR:     replicate(!VALUES.F_NaN, len),  $
            SII_FLUX1: replicate(!VALUES.F_NaN, len),  $
            SII_FLUX2: replicate(!VALUES.F_NaN, len),  $
            HA_FLUX:   replicate(!VALUES.F_NaN, len),  $
            HB_FLUX:   replicate(!VALUES.F_NaN, len),  $
            OIII_FLUX: replicate(!VALUES.F_NaN, len),  $
            CHI_SPEC:  replicate(!VALUES.F_NaN, len),  $
            CHI_PHOT:  replicate(!VALUES.F_NaN, len)}
  chains = replicate(chains, n_elements(files))

  for ii = 0, n_elements(chains) - 1 do begin
     if NOT folded then begin
        case ii of
           0: begin
              region = 'OPTFL'
              file = optf
           end
           1: begin
              region = 'INNFL'
              file = innf
           end
           2: begin
              region = 'INTUP'
              file = intfu
           end
           3: begin
              region = 'INTDN'
              file = intfd
           end
           4: begin
              region = 'OUTUP'
              file = outfu
           end
           5: begin
              region = 'OUTDN'
              file = outfd
           end
        endcase
     endif else begin
       case ii of
           0: begin
              region = 'OPTFL'
              file = optf
           end
           1: begin
              region = 'INNFL'
              file = innf
           end
           2: begin
              region = 'INTER'
              file = intf
           end
           3: begin
              region = 'OUTER'
              file = outf
           end
        endcase
     endelse

     spawn, 'cat '+repstr(file, 'post_equal_weights.dat', '.analyze')+' | grep z > foo.txt'
     t = file_lines('foo.txt')
     
     if t gt 0 then begin
        readcol, file, $
                 siiew, z, siiew2, oiiiew, hew, $
                 av, age, ltau, chis, chip, $
                 lmass, lssfr, sfr, $
                 haf, hbf, $
                 oiiif, siif, siif2
     endif else begin
        readcol, file, $
                 siiew, siiew2, oiiiew, hew, $
                 av, age, ltau, chis, chip, $
                 lmass, lssfr, sfr, $
                 oiiif, siif, siif2, haf, hbf
     endelse
     tlen = n_elements(siiew)
     
     chains[ii].REGION              = region
     chains[ii].Z                   = z
     chains[ii].SII_EW1[0:tlen-1]   = siiew   
     chains[ii].SII_EW2[0:tlen-1]   = siiew2   
     chains[ii].H_EW[0:tlen-1]      = hew     
     chains[ii].OIII_EW[0:tlen-1]   = oiiiew   
     chains[ii].AGE[0:tlen-1]       = alog10(age)
     chains[ii].LTAU[0:tlen-1]      = ltau       
     chains[ii].AV[0:tlen-1]        = av        
     chains[ii].LMASS[0:tlen-1]     = lmass     
     chains[ii].SFR[0:tlen-1]       = sfr       
     chains[ii].LSSFR[0:tlen-1]     = lssfr     
     chains[ii].SII_FLUX1[0:tlen-1] = siif 
     chains[ii].SII_FLUX2[0:tlen-1] = siif2 
     chains[ii].HA_FLUX[0:tlen-1]   = haf   
     chains[ii].HB_FLUX[0:tlen-1]   = hbf   
     chains[ii].OIII_FLUX[0:tlen-1] = oiiif 
     chains[ii].CHI_SPEC[0:tlen-1]  = chis      
     chains[ii].CHI_PHOT[0:tlen-1]  = chip      

  endfor

  cols = ['777777'x, 255, '00a500'x, '00ff00'x, 'ffa500'x, 'ff0000'x]
  plot, chains[indgen(n_elements(chains)-1)+1].AGE, $
        chains[indgen(n_elements(chains)-1)+1].LTAU, psym = 3, $
        /nodat, /ynoz;, xran = [0,2.5], yran = [0,2.5], /iso
  oplot, chains[1].Age, chains[1].LTAU, psym = 3, col = 255
  if NOT folded then begin
     oplot, chains[2].Age, chains[2].LTAU, psym = 3, col = '00a500'x
     oplot, chains[3].Age, chains[3].LTAU, psym = 3, col = '00ff00'x
     oplot, chains[4].Age, chains[4].LTAU, psym = 3, col = 'ffa500'x
     oplot, chains[5].Age, chains[5].LTAU, psym = 3, col = 'ff0000'x
  endif else begin
     oplot, chains[2].Age, chains[2].LTAU, psym = 3, col = '00ff00'x
     oplot, chains[3].Age, chains[3].LTAU, psym = 3, col = 'ffa500'x
  endelse
  oplot, chains[0].Age, chains[0].LTAU, psym = 1

  mwrfits, chains, output, /create
  
end

pro foo

  doexpcovar, 'MACS0744/00660_2_pyspecfitPhotResults', $
              'resultsSummaries/00660_2_exp_chains.fits', /folded

  doexpcovar, 'MACS1149/00900_1_pyspecfitPhotResults', $
              'resultsSummaries/00900_1_exp_chains.fits', /folded

  doexpcovar, 'MACS2129/00451_2_pyspecfitPhotResults', $
              'resultsSummaries/00451_2_exp_chains.fits', /folded

  doexpcovar, 'MACS1423/01916_2_pyspecfitPhotResults', $
              'resultsSummaries/01916_2_exp_chains.fits', /folded

end
