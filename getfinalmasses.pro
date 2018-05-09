pro getfinalmasses, lgnlist, explist

  summaryData = mrdfits('objectSummaryData.fits', 1)
  
  readcol, explist, expfiles, f = 'A'
  readcol, lgnlist, lgnfiles, f = 'A'
  
  ngals = n_elements(lgnfiles)

  test     = mrdfits(lgnfiles[0], 1)

  regions  = test.REGION
  nregions = n_elements(regions)
  
  masses  = fltarr(ngals, nregions, 2)
  emasses = fltarr(ngals, nregions, 2)

  sfrs  = fltarr(ngals, nregions, 2)
  esfrs = fltarr(ngals, nregions, 2)

  files = strarr(ngals)
  
  for ii = 0, ngals - 1 do begin
     for jj = 0, 1 do begin
        case jj of
           0: begin
              file = lgnfiles[ii]
              print, file
              files[ii] = strmid(file, 0, 7)
              qui = where(summaryData.ID eq files[ii])
              mag = summaryData[qui].MU
              mag = mag[0]
              emag = summaryData[qui].EMU
              emag = emag[0]
           end
           1: file = expfiles[ii]
        endcase

        data = mrdfits(file, 1, /silent)
        masses[ii,*,jj]  = data.LMASS[0] - alog10(mag)
        emasses[ii,*,jj] = sqrt((0.5 * (data.LMASS[2] - data.LMASS[1]))^2 + $
                                (1./alog(10) * emag/mag)^2)

        sfrs[ii,*,jj] = alog10(data.SFR[0]) - alog10(mag)
        esfrs[ii,*,jj] = sqrt((0.5 * (data.SFR[2] - data.SFR[1])/data.SFR[0]/alog(10))^2 + (emag/mag/alog(10))^2)
        
     endfor
  endfor

  !p.multi = [0,2,0]
  
  plot, masses, masses, /nodat, yran = [-1,1]/2.
  qui = where(regions eq 'OPTFL')
  qui2 = where(regions eq 'INNFL' $
               OR $
               regions eq 'INTER' $
               OR $
               regions eq 'OUTER')
  cols = [255,'00a500'x, 'ff5500'x]
  oploterror, masses[*,*,0], masses[*,*,1] - masses[*,*,0], $
              emasses[*,*,0], sqrt(emasses[*,*,0]^2 + emasses[*,*,1]^2), psym = 1
  oploterror, masses[*,qui,0], masses[*,qui,1] - masses[*,qui,0], $
              emasses[*,qui,0], sqrt(emasses[*,qui,0]^2 + emasses[*,qui,1]^2), $
              col= '777777'x, errcol = '777777'x, psym = 1
  for ii = 0,2 do $
     oploterror, masses[*,qui2[ii],0], masses[*,qui2[ii],1] - masses[*,qui2[ii],0], $
                 emasses[*,qui2[ii],0], sqrt(emasses[*,qui2[ii],0]^2 + emasses[*,qui2[ii],1]^2), $
                 col= cols[ii], errcol = cols[ii], psym = 1
  oplot, !X.CRANGE, [0,0], col = 255
  fmasses = mean(masses, dim = 3)
  efmasses = mean(emasses, dim = 3) ;; not independent measurements, so can't do the root 2 thing

  oploterror, fmasses[*,qui], alog10(total(10.^(masses[*,qui2,0] + fan(alog10([1,2,2]),4)), 2)) - fmasses[*,qui], $
              efmasses[*,qui], sqrt(total(emasses[*,qui2,0]^2,2) + efmasses[*,qui]^2), $
              psym = 1, col = long('00a5ff'x), errcol = long('00a5ff'x)
  oploterror, fmasses[*,qui], alog10(total(10.^(masses[*,qui2,1] + fan(alog10([1,2,2]),4)), 2)) - fmasses[*,qui], $
              efmasses[*,qui], sqrt(total(emasses[*,qui2,1]^2,2) + efmasses[*,qui]^2), $
              psym = 1, col = long('0055ff'x), errcol = long('0055ff'x)
    
  
  print, mean(masses[*,*,1] - masses[*,*,0])
  print, stddev(masses[*,*,1] - masses[*,*,0])
  print, minmax(masses[*,*,1] - masses[*,*,0])
  print, median((masses[*,*,1] - masses[*,*,0])/sqrt(emasses[*,*,0]^2 + emasses[*,*,1]^2))
  print, minmax((masses[*,*,1] - masses[*,*,0])/sqrt(emasses[*,*,0]^2 + emasses[*,*,1]^2))
  print, files
  print, alog10(total(10.^(masses[*,qui2,0] + fan(alog10([1,2,2]),4)), 2)) - masses[*,qui,0]
  print, (alog10(total(10.^(masses[*,qui2,0] + fan(alog10([1,2,2]),4)), 2)) - masses[*,qui,0]) $
         / sqrt(total(emasses[*,qui2,0]^2,2) + emasses[*,qui,0]^2)
  
  plot, sfrs, sfrs, /nodat, xran = [-2,2], yran = [-2,2];/2.
  qui = where(regions eq 'OPTFL')
  qui2 = where(regions eq 'INNFL' $
               OR $
               regions eq 'INTER' $
               OR $
               regions eq 'OUTER')
  cols = [255,'00a500'x, 'ff5500'x]
  oploterror, sfrs[*,*,0], sfrs[*,*,1] - sfrs[*,*,0], $
              esfrs[*,*,0], sqrt(esfrs[*,*,0]^2 + esfrs[*,*,1]^2), psym = 1
  oploterror, sfrs[*,qui,0], sfrs[*,qui,1] - sfrs[*,qui,0], $
              esfrs[*,qui,0], sqrt(esfrs[*,qui,0]^2 + esfrs[*,qui,1]^2), $
              col= '777777'x, errcol = '777777'x, psym = 1
  for ii = 0,2 do $
     oploterror, sfrs[*,qui2[ii],0], sfrs[*,qui2[ii],1] - sfrs[*,qui2[ii],0], $
                 esfrs[*,qui2[ii],0], sqrt(esfrs[*,qui2[ii],0]^2 + esfrs[*,qui2[ii],1]^2), $
                 col= cols[ii], errcol = cols[ii], psym = 1
  oplot, !X.CRANGE, [0,0], col = 255
  fsfrs = mean(sfrs, dim = 3)
  efsfrs = mean(esfrs, dim = 3) ;; not independent measurements, so can't do the root 2 thing

  oploterror, sfrs[*,qui,0], alog10(total(10.^(sfrs[*,qui2,0] + fan(alog10([1,2,2]),4)), 2)) - sfrs[*,qui,0], $
              efsfrs[*,qui,0], sqrt(total(esfrs[*,qui2,0]^2,2) + esfrs[*,qui,0]^2), $
              psym = 1, col = long('00a5ff'x), errcol = long('00a5ff'x)
  oploterror, sfrs[*,qui,1], alog10(total(10.^(sfrs[*,qui2,1] + fan(alog10([1,2,2]),4)), 2)) - sfrs[*,qui,1], $
              esfrs[*,qui,1], sqrt(total(esfrs[*,qui2,1]^2,2) + esfrs[*,qui,1]^2), $
              psym = 1, col = long('0055ff'x), errcol = long('0055ff'x)

  print, ''
  print, ' - - - - - - - - - - - - - - - - - '
  print, ''
  print, mean(sfrs[*,*,1] - sfrs[*,*,0])
  print, stddev(sfrs[*,*,1] - sfrs[*,*,0])
  print, minmax(sfrs[*,*,1] - sfrs[*,*,0])
  print, median((sfrs[*,*,1] - sfrs[*,*,0])/sqrt(esfrs[*,*,0]^2 + esfrs[*,*,1]^2))
  print, minmax((sfrs[*,*,1] - sfrs[*,*,0])/sqrt(esfrs[*,*,0]^2 + esfrs[*,*,1]^2))
  print, files
  print, alog10(total(10.^(sfrs[*,qui2,0] + fan(alog10([1,2,2]),4)), 2)) - sfrs[*,qui,0]
  print, (alog10(total(10.^(sfrs[*,qui2,0] + fan(alog10([1,2,2]),4)), 2)) - sfrs[*,qui,0]) $
         / sqrt(total(esfrs[*,qui2,0]^2,2) + esfrs[*,qui,0]^2)
  
  
  ;; all the OPTFL measurements agree to within 1-sigma
  ;; Average and print; use for rest of analysis
  ;; no other measurements disagree by more than 2.35 sigma
  ;; with peak-to-peak deltas at any extraction of -0.244 to 0.057 dex
  ;; = 0.3 dex total
  ;;
  ;; This also holds for INNFL and INTER regions. Only in OUTER are
  ;; measurements only in > 1 sigma disagreement (2.35), so you
  ;; *could* average them all, but perhaps just do the OPTFL one so
  ;; you have consistent total masses.
  ;;
  ;; Sum of the INNFL, INTER, OUTER and OPTFL also agree well
  ;; So it's unclear what to use for a "total" mass.
  ;; Print both
  ;;
  ;; Errors are about 2x as big using sum as compared to optimal extraction
  
  savedata = {ID:     string(0), $
              LMASS:  0., $
              ELMASS: 0., $
              TMASS: 0., $
              ETMASS: 0.}
  savedata = replicate(savedata, ngals)
  for ii = 0, ngals - 1 do begin
     savedata[ii].ID = files[ii]
     savedata[ii].LMASS  = fmasses[ii,qui]
     savedata[ii].ELMASS = efmasses[ii,qui]
     savedata[ii].TMASS  = mean([alog10(total(10.^(masses[ii,qui2,0] + alog10([1,2,2])))), alog10(total(10.^(masses[ii,qui2,1] + alog10([1,2,2]))))] )
     savedata[ii].ETMASS = mean([sqrt(total(emasses[ii,qui2,0]^2 + alog10([1,2,2]))), sqrt(total(emasses[ii,qui2,1]^2 + alog10([1,2,2])))])
  endfor
  mwrfits, savedata, 'finalMasses.fits', /create

  
end
;getfinalmasses, 'study5_lgn_Bestfits.list', 'study5_exp_Bestfits.list'
