pro doeverythingRadialExtract_A370_M0416

  doEverything, 'ABEL0370', $
                zmin = 1.0, $
                zmax = 1.8, $
                minmag = 14., $
                maxmag = 21.8, $
                outprefix = '~/Projects/GLASS/spectralAnalysis', $
                /delete_prev

  doEverything, 'MACS0416', $
                zmin = 1.0, $
                zmax = 1.8, $
                minmag = 14., $
                maxmag = 21.8, $
                outprefix = '~/Projects/GLASS/spectralAnalysis', $
                /delete_prev

end

;;
;;
;;

pro doEverything, field, $
                  ZMIN = zmin, $
                  ZMAX = zmax, $
                  MINMAG = minmag, $
                  MAXMAG = maxmag, $
                  BOTHPA = bothPa, $
                  OUTPREFIX = outprefix, $
                  DELETE_PREV = delete_prev

  if NOT keyword_Set(DELETE_PREV) then $
     delete_prev = 0 $
  else $
     delete_prev = 1
  
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
  cat_suffix = '_allAvailableDataTaka.fits' ;; Modified to accomodate cats w/ Taka's stuff, not Kuang's.

  case field of
     'a370'    : field = 'ABEL0370'
     'a0370'   : field = 'ABEL0370'
     'abel370' : field = 'ABEL0370'
     'abel0370': field = 'ABEL0370'
     'ABEL370' : field = 'ABEL0370'
     'ABEL0370': field = 'ABEL0370'
     
     'm0416'   : field = 'MACS0416'
     'macs0416': field = 'MACS0416'
     'MACS0416': field = 'MACS0416'
  endcase

  cat = data_dir+'/'+field+cat_suffix

  zzmin = string(zmin  , f = '(F3.1)')
  zzmax = string(zmax  , f = '(F3.1)')
  lomag = string(minmag, f = '(F4.1)')
  himag = string(maxmag, f = '(F4.1)')

  cclevel = 'c99.9'  ;; MACS0416 is missing a GiG Cat...
  bbpa    = 'pacut0' ;; Too few sources pass both PA cuts
  
  outsuffix = '/sourceLists/'+field+'_'+zzmin+'_'+zzmax+'_'$
              +lomag+'_'+himag+'_'$
              +cclevel+'_'+bbpa

  if NOT keyword_set(OUTPREFIX) then outprefix = '.'

  cutOnZTaka, cat, zmin = zmin, zmax = zmax, $
              /photoz, $
              output = 'tmpZ.list'
  cutOnMagTaka, 'tmpZ.list', minMag = minmag, maxMag = maxmag, $
                output = outprefix+'/'+outsuffix+'_FILELIST.list'
  
  toGet = dump1field(outprefix+'/'+outsuffix+'_FILELIST.list', field)
  ntoget = n_elements(toGet.Z)

  specOutDir = outprefix+'/'+field
  spawn, 'mkdir '+specOutDir
  if delete_prev then $
     spawn, 'rm -r '+specOutDir+'/*'
  
  if ntoget ge 1 then begin
     for ii = 0, ntoget - 1 do begin
        
        check1 = file_info(toget.PA1B[ii])
        if check1.EXISTS then begin
           toWrite1b = extract_1_indep(toget.PA1B[ii], toGet.Z[ii], Z_QUAL = toget.ZQ[ii], $
                                       Z_PHOT = toget.Z_PHOT[ii], PA = toget.PA1) ;, /refit)
           toWrite1r = extract_1_indep(toget.PA1R[ii], toGet.Z[ii], Z_QUAL = toget.ZQ[ii], $
                                       Z_PHOT = toget.Z_PHOT[ii], PA = toget.PA1) ;, /refit)
           mwrfits, toWrite1b, specOutDir+'/'+strcompress(string(toGet.ID[ii], f = '(I05)')+'_1_B.fits', /rem), /create
           mwrfits, toWrite1r, specOutDir+'/'+strcompress(string(toGet.ID[ii], f = '(I05)')+'_1_R.fits', /rem), /create
        endif

        check2 = file_info(toget.PA2B[ii])
        if check2.EXISTS then begin
           toWrite2b = extract_1_indep(toget.PA2B[ii], toGet.Z[ii], Z_QUAL = toget.ZQ[ii], $
                                       Z_PHOT = toget.Z_PHOT[ii], PA = toget.PA2) ;, /refit)
           toWrite2r = extract_1_indep(toget.PA2R[ii], toGet.Z[ii], Z_QUAL = toget.ZQ[ii], $
                                       Z_PHOT = toget.Z_PHOT[ii], PA = toget.PA2) ;, /refit)
           mwrfits, toWrite2b, specOutDir+'/'+strcompress(string(toGet.ID[ii], f = '(I05)')+'_2_B.fits', /rem), /create
           mwrfits, toWrite2r, specOutDir+'/'+strcompress(string(toGet.ID[ii], f = '(I05)')+'_2_R.fits', /rem), /create
        endif
        
     endfor
  endif
  
end
