pro compsfhres, expDirInlist, lgnDirInlist

  readcol, expDirInlist, $
           expdirs, f = 'A'
  readcol, lgnDirInlist, $
           lgndirs, f = 'A'

  ngals = n_elements(lgndirs)
  
  for ii = 0, ngals - 1 do begin
     analyzetofits, expdirs[ii], 'efoo.fits'
     analyzetofits, lgndirs[ii], 'lfoo.fits', $
                    sfh = 'LogNormal'

     stop
     
  endfor

  
  
end
