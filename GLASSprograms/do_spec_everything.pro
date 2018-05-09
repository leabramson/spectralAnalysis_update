;; Cull and get the files for 1 field

pro do_spec_everything, field, $
                        ZMIN = zmin, $
                        ZMAX = zmax, $
                        MINMAG = minmag, $
                        MAXMAG = maxmag, $
                        DOCONTAM = docontam, $
                        CLEVEL = clevel, $
                        BOTHPA = bothPa, $
                        OUTPREFIX = outprefix
  
  if NOT keyword_set(ZMIN) then zmin = 1.0
  if NOT keyword_set(ZMAX) then zmax = 1.8
  if NOT keyword_set(MINMAG) then minmag = 14
  if NOT keyword_set(MAXMAG) then maxmag = 21.5
  if NOT keyword_set(DOCONTAM) then docontam = 0 else docontam = 1
  if docontam then begin
     if NOT keyword_set(CLEVEL) then clevel = 1.5
     if NOT keyword_set(BOTHPA) then bothPA = 0 else bothPA = 1
  endif
  
  data_dir   = '../MasterCats_V1'
  cat_suffix = '_allAvailableData.fits' ;'_MASTER_z_GiG_Kpz.fits'


  case field of
     'a370'    : field = 'ABEL0370'
     'a0370'   : field = 'ABEL0370'
     'abel370' : field = 'ABEL0370'
     'abel0370': field = 'ABEL0370'
     'ABEL370' : field = 'ABEL0370'
     'ABEL0370': field = 'ABEL0370'

     'a2744'   : field = 'ABEL2744'
     'abel2744': field = 'ABEL2744'
     'ABEL2744': field = 'ABEL2744'
     
     'm0416'   : field = 'MACS0416'
     'macs0416': field = 'MACS0416'
     'MACS0416': field = 'MACS0416'
     
     'm0717'   : field = 'MACS0717'
     'macs0717': field = 'MACS0717'
     'MACS0717': field = 'MACS0717'
     
     'm0744'   : field = 'MACS0744'
     'macs0744': field = 'MACS0744'
     'MACS0744': field = 'MACS0744'

     'm1149'   : field = 'MACS1149'
     'macs1149': field = 'MACS1149'
     'MACS1149': field = 'MACS1149'

     'm1423'   : field = 'MACS1423'
     'macs1423': field = 'MACS1423'
     'MACS1423': field = 'MACS1423'
     
     'm2129'   : field = 'MACS2129'
     'macs2129': field = 'MACS2129'
     'MACS2129': field = 'MACS2129'

     'r1347'   : field = 'RXJC1347'
     'rxj1347' : field = 'RXJC1347'
     'rxjc1347': field = 'RXJC1347'
     'RXJ1347' : field = 'RXJC1347'
     'RXJC1347': field = 'RXJC1347'

     'r2248'   : field = 'RXJC2248'
     'rxj2248' : field = 'RXJC2248'
     'rxjc2248': field = 'RXJC2248'
     'RXJ2248' : field = 'RXJC2248'
     'RXJC2248': field = 'RXJC2248'
  endcase

  cat = data_dir+'/'+field+cat_suffix

  zzmin = string(zmin  , f = '(F3.1)')
  zzmax = string(zmax  , f = '(F3.1)')
  lomag = string(minmag, f = '(F4.1)')
  himag = string(maxmag, f = '(F4.1)')
  
  if docontam then begin
     cclevel = 'c'+string(clevel, f = '(F3.1)')
     bbpa    = 'pacut'+string(bothPa, f = '(I1)')
  endif else begin
     cclevel = 'c99.9'
     bbpa    = 'pacut0'
  endelse
  
  outsuffix = '/sourceLists/'+field+'_'+zzmin+'_'+zzmax+'_'$
              +lomag+'_'+himag+'_'$
              +cclevel+'_'+bbpa

  if NOT keyword_set(OUTPREFIX) then outprefix = '.'

  cutOnZ, cat, zmin = zmin, zmax = zmax, $
          /photoz, $
          docontam = docontam, clevel = clevel, bothPa = bothpa, $
          output = 'tmpZ.list'
  cutOnMag, 'tmpZ.list', minMag = minmag, maxMag = maxmag, $
            output = outprefix+'/'+outsuffix+'_FILELIST.list'
  
  toGet = dump1field(outprefix+'/'+outsuffix+'_FILELIST.list', field)
  ntoget = n_elements(toGet.Z)

  specOutDir = outprefix+'/'+field
  spawn, 'mkdir '+specOutDir
  if ntoget ge 1 then begin
     for ii = 0, ntoget - 1 do begin
        
        check = file_info(toget.PA1B[ii])
        if check.exists then begin
           toWrite1b = extract_1(toget.PA1B[ii], toGet.Z[ii], $
                                 Z_PHOT = toget.Z_PHOT[ii], PA = toget.PA1)   ;, /refit)
           toWrite2b = extract_1(toget.PA2B[ii], toGet.Z[ii], $
                                 Z_PHOT = toget.Z_PHOT[ii], PA = toget.PA2)   ;, /refit)
           toWrite1r = extract_1(toget.PA1R[ii], toGet.Z[ii], $
                                 Z_PHOT = toget.Z_PHOT[ii], PA = toget.PA1)   ;, /refit)
           toWrite2r = extract_1(toget.PA2R[ii], toGet.Z[ii], $
                                 Z_PHOT = toget.Z_PHOT[ii], PA = toget.PA2) ;, /refit)
           
           mwrfits, toWrite1b, specOutDir+'/'+strcompress(string(toGet.ID[ii], f = '(I05)')+'_1_B.fits', /rem), /create
           mwrfits, toWrite2b, specOutDir+'/'+strcompress(string(toGet.ID[ii], f = '(I05)')+'_2_B.fits', /rem), /create
           mwrfits, toWrite1r, specOutDir+'/'+strcompress(string(toGet.ID[ii], f = '(I05)')+'_1_R.fits', /rem), /create
           mwrfits, toWrite2r, specOutDir+'/'+strcompress(string(toGet.ID[ii], f = '(I05)')+'_2_R.fits', /rem), /create
        endif
        
     endfor
  endif
  
end

pro doclean

  spawn, 'ls ~/Projects/GLASS/MasterCats_V1/*allAvail*.fits > catalogs.list'
  readcol, 'catalogs.list', cats, f = 'A'
  ncats = n_elements(cats)

  cat2 = strarr(ncats)
  for ii = 0, ncats - 1 do $
     cat2[ii] = strmid(cats[ii], strpos(cats[ii], 'allAvailableData')-9, 8)
  
  print, ''
  print, 'Extracing spectra for '+string(ncats, f = '(I02)')+' GLASS fields'
  for ii = 0, ncats - 1 do $
     print, ' >>> '+cat2[ii]
  print, ''

  wait, 1

  for ii = 0, ncats - 1 do $
  doeverything, cat2[ii], $
                zmin = 1.0, zmax = 1.8, $
                minmag = 14, maxmag = 21.8, $
                /docontam, /bothpa, clevel = 2.0, $
                outprefix = '~/Projects/GLASS/spectralAnalysis'

  print, ''
  print, ''
  print, ''
  print, '*** DONE! ***'
  print, ''
  print, ''
  print, ''
    
end
