pro unifyMultiple, dirlist

  readcol, dirlist, dirs, f = 'A'
  ndirs = n_elements(dirs)

  for ii = 0, ndirs - 1 do begin
     dir  = strmid(dirs[ii], 0, 5)
     targ = strmid(dirs[ii], 6, strpos(dirs[ii], "_") - 6)

     b1 = dir+'/'+targ+'_1_B.fits'
     r1 = dir+'/'+targ+'_1_R.fits'
     b2 = dir+'/'+targ+'_2_B.fits'
     r2 = dir+'/'+targ+'_2_R.fits'

     o1 = dirs[ii]+'/'+targ+'_1_unified_1D.fits'
     o2 = dirs[ii]+'/'+targ+'_2_unified_1D.fits'
     
     unifySpectra, b1, r1, output = o1
     unifySpectra, b2, r2, output = o2

     print, ''
     print, 'Unified spectra for '+dirs[ii]+' PA1 written to : '+o1
     print, 'Unified spectra for '+dirs[ii]+' PA2 written to : '+o2
     print, ''
     
  endfor

end
;unifyMultiple, 'resultDirs.list'

