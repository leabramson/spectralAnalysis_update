pro organizeonefield, field

  spawn, 'ls '+field+'/?????_?_?.fits > '+field+'/files.list'
  readcol, field+'/files.list', files, f = 'A'
  nfiles = n_elements(files)

  tfiles = strarr(nfiles)
  IDs    = strarr(nfiles)
  grisms = strarr(nfiles)
  PAs    = strarr(nfiles)

  for ii = 0, nfiles - 1 do begin

     tfile = repstr(files[ii], field+'/', '')

     tfiles[ii] = tfile
     IDs[ii]    = strmid(tfile, 0, 5)
     PAs[ii]    = strmid(tfile, 6, 1)
     grisms[ii] = strmid(tfile, 8, 1)
     
  end

  sIDs = IDs[sort(IDs)]
  uIDs = sIDs[UNIQ(sIDs)]
  ngals = n_elements(uIDs)
    
  close, 1
  openw, 1, field+'/fileInfo.list', width = 256
  printf, 1, '#FILENAME GLASS_ID PA GRISM Z_PHOT Z_SPEC Z_SPEC_QUAL'
  for ii = 0, ngals - 1 do begin
     hits = where(IDs eq uIDs[ii], nhits)
     td = mrdfits(field+'/'+tfiles[hits[0]], 1)
     for jj = 0, nhits - 1 do $
        printf, 1, tfiles[hits[jj]], ' ', $
                IDs[hits[jj]],       ' ', $
                PAs[hits[jj]],       ' ', $
                grisms[hits[jj]],    ' ', $
                td.Z_PHOT,           ' ', $
                td.Z,                ' ', $
                td.Z_QUAL
                
  endfor
  close, 1
  
end

;;
;;
;;

pro organizeAll

  organizeonefield, 'MACS0717'
  organizeonefield, 'MACS0744'
  organizeonefield, 'MACS1149'
  organizeonefield, 'MACS1423'
  organizeonefield, 'MACS2129'
  organizeonefield, 'RXJC1347'
  organizeonefield, 'RXJC2248'
  organizeonefield, 'ABEL2744'
  
end
