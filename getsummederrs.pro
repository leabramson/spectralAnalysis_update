pro getsummederrs

  readcol, 'study5_lgn_Bestfits.list', lgnfiles, f= 'A'
  readcol, 'study5_exp_Bestfits.list', expfiles, f= 'A'

  ngals = n_elements(lgnfiles)

  smerrs = fltarr(ngals, 2)
  sserrs = fltarr(ngals, 2)

  for ii = 0, ngals - 1 do begin

     d = mrdfits(lgnfiles[ii], 1)

     inn = where(d.REGION eq 'INNFL')
     int = where(d.REGION eq 'INTER')
     out = where(d.REGION eq 'OUTER')

     innE = (d[inn].LMASS[2] - d[inn].LMASS[1]) / 2
     intE = (d[int].LMASS[2] - d[int].LMASS[1]) / 2
     outE = (d[out].LMASS[2] - d[out].LMASS[1]) / 2

     smerrs[ii,0] = sqrt(innE^2 + 2 * intE^2 + 2 * outE^2)
     
     d = mrdfits(expfiles[ii], 1)

     inn = where(d.REGION eq 'INNFL')
     int = where(d.REGION eq 'INTER')
     out = where(d.REGION eq 'OUTER')

     innE = (d[inn].LMASS[2] - d[inn].LMASS[1]) / 2
     intE = (d[int].LMASS[2] - d[int].LMASS[1]) / 2
     outE = (d[out].LMASS[2] - d[out].LMASS[1]) / 2

     smerrs[ii,1] = sqrt(innE^2 + 2 * intE^2 + 2 * outE^2)

     print, lgnfiles[ii], smerrs[ii,0], smerrs[ii,1]
     
  endfor

end
