function retrieve1ID, fieldName, id
  
  ;; Set-up file structure stuff
;  prefix =
;  '/Volumes/LAbramiDrive/GLASS_VIS_INSPECT/'+fieldName+'/IndividualObjects/'
  prefix = '/home/labramson/LocalData/Abramson/GLASS_official_products/'+fieldName+'/IndividualObjects/'
  prelen = strlen(prefix)
  spawn, 'ls '+prefix+'*-???-150515_0????-G102.2D.fits > tmp.list'
  readcol, 'tmp.list', tfiles, f = 'A', /silent
  ntfiles = n_elements(tfiles)
  tfile = tfiles[0]
  tfile = strmid(tfile, prelen)
  vNo = strpos(tfile, '150515')
  pas = strarr(ntfiles)
  for ii = 0, ntfiles - 1 do $
     pas[ii] = strmid(tfiles[ii], vNo-4+prelen, 3)
  pas = pas[sort(pas)]
  pas = pas[UNIQ(pas)]

  fieldFileStart = strmid(tfile, 0, vNo-5)
  toGet = prefix+fieldFileStart+'-PPP-150515_0XXXX-G1YY.2D.fits'
  
  ;; Build the G102 and G141 file lists for each PA
  PA1_filelist = strarr(2)
  PA2_filelist = strarr(2)
  fields       = string(0)

  fname = repstr(toGet, 'XXXX', strcompress(string(id, f= '(I04)'), /rem))

  pa1_fname = repstr(fname, 'PPP', pas[0])
  fname_blue = repstr(pa1_fname, 'YY', '02')
  fname_red  = repstr(pa1_fname, 'YY', '41')
  PA1_filelist = [fname_blue, fname_red]

  pa2_fname = repstr(fname, 'PPP', pas[1])
  fname_blue = repstr(pa2_fname, 'YY', '02')
  fname_red  = repstr(pa2_fname, 'YY', '41')
  PA2_filelist = [fname_blue, fname_red]

;  redshifts[ii] = z[ii]
;  qflags[ii]    = zq[ii]
;  pzs[ii]       = pz[ii]
;  mags[ii]      = mag[ii]
;  ras[ii]       = ra[ii]
;  decs[ii]      = dec[ii]
;  fields[ii]    = field[ii]
  

  savedata = {PA1B  : reform(PA1_filelist[0,*]), $
              PA2B  : reform(PA2_filelist[0,*]), $
              PA1R  : reform(PA1_filelist[1,*]), $
              PA2R  : reform(PA2_filelist[1,*]), $
              FIELD : fieldName, $
              ID    : id, $
              PA1   : pas[0], $
              PA2   : pas[1]}
  
  RETURN, savedata
  
end

