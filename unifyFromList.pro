pro unifyFromList, inlist

  readcol, inlist, files, f = 'A'
  nfiles = n_elements(files)

  tids = strarr(nfiles)
  for ii = 0, nfiles - 1 do $
     tids[ii] = strmid(files[ii], 0, strpos(files[ii], '_'))

  sids  = sort(tids)
  tids  = tids[sids]
  uids = tids[UNIQ(tids)]
  
  ii = 0
  counter = 0
  while ii le nfiles - 1 do begin

     targ = uids[counter]
     nuse = n_elements(where(tids eq targ))
     
     b1 = targ+'_1_B.fits'
     r1 = targ+'_1_R.fits'
     b2 = targ+'_2_B.fits'
     r2 = targ+'_2_R.fits'
     
     o1 = targ+'_1_unified_1D.fits'
     o2 = targ+'_2_unified_1D.fits'
     
     unifySpectra, b1, r1, output = o1
     unifySpectra, b2, r2, output = o2

     print, ''
     print, 'Unified spectra for obj. '+targ+' PA1 written to : '+o1
     print, 'Unified spectra for obj. '+targ+' PA2 written to : '+o2
     print, ''

     ii += nuse
     counter++
     
  endwhile
  

end





;unifyMultiple, 'faintResultDirs.list'
