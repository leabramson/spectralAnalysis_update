pro fitSpectra, dir, PA, $
                ncores = ncores

  if NOT keyword_set(ncores) then ncores = 12

  if PA eq 1 then $
     spawn, 'ls '+dir+'/*_1_B_inner.pyspec > tmpObjects.list' $
  else $
     spawn, 'ls '+dir+'/*_2_B_inner.pyspec > tmpObjects.list'

  readcol, 'tmpObjects.list', files, f = 'A'
  nfiles = n_elements(files)
  
  base = strarr(nfiles)
  for jj = 0, nfiles - 1 do begin
     base[jj] = strmid(files[jj], strpos(files[jj], '/'), strpos(files[jj],' _B'))

     spawn, 'mkdir '+dir+'/'+base[jj]
     spawn, 'mpirun -np '+string(ncores, f = '(I2)')+' python fitspec.py '+$
            base[jj]+' INNER 1 1.0 '+base[jj]+'/inner'
     spawn, 'cat '+base[jj]+'/inner.analyze | grep z > ztmp.list'
     readcol, 'ztmp.list', zmed, zlow, zhi, f = 'X,F,F,F'

     spawn, 'mpirun -np '+string(ncores, f = '(I2)')+' python fitspec.py '+$
            base[jj]+' INTER '+string(zmed, f = '(F5.3)')+' 0.005 '+base[jj]+'/inter'
     spawn, 'mpirun -np '+string(ncores, f = '(I2)')+' python fitspec.py '+$
            base[jj]+' OUTER '+string(zmed, f = '(F5.3)')+' 0.005 '+base[jj]+'/outer'

  endfor
  

end

