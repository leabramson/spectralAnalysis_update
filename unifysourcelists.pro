pro unifysourcelists, sourcelistdir, output

  spawn, 'ls '+sourcelistdir+' > sourceLists.dat'
  readcol, 'sourceLists.dat', lists, f = 'A'
  nlists = n_elements(lists)

  for ii = 0, nlists - 1 do begin
     if ii eq 0 then $
        spawn, 'cat '+sourcelistdir+'/'+lists[ii]+' > '+sourcelistdir+'/'+output $
     else $
        spawn, 'cat '+sourcelistdir+'/'+lists[ii]+' >> '+sourcelistdir+'/'+output
  endfor

end
;; unifysourcelists, 'sourceLists', '1.0_1.8_14.0_21.8_c99.9_pacut0_ALLSRCS.list'

