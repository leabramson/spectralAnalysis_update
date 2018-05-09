pro get_1ID_spectrum, fieldName, ra, dec, $
                      OUTDIR = outdir, $
                      OUTNAME = outname

  if NOT keyword_set(OUTDIR) then outdir = './'

  prefix = '/home/labramson/Projects/GLASS/MasterCats_V1/'
  mastercat = prefix+fieldName+'_allAvailableData.fits'

  if fieldName eq 'ABEL0370' $
     OR  $
     fieldName eq 'MACS0416' then $
        mastercat = prefix+fieldname+'_glass_srclist.fits'
  
  glass_ID = getoneID(ra, dec, mastercat, matchlen = 0.5)
  
  files = retrieve1ID(fieldName, glass_ID)

  if NOT keyword_set(OUTNAME) then outname = string(glass_ID, f = '(I05)')
  outprefix = outdir+outname
  
  check = file_info(files.PA1B)
  if check.exists then begin
     blueSpec_1 = extract_1_indep(files.PA1B, PA = files.PA1)
     redSpec_1  = extract_1_indep(files.PA1R, PA = files.PA1)
     mwrfits, blueSpec_1, outprefix+'_G102_'+string(files.PA1, f = '(I03)')+'.fits', /create
     mwrfits,  redSpec_1, outprefix+'_G141_'+string(files.PA1, f = '(I03)')+'.fits', /create
  endif else $
     print, 'NO EXTRACTION PRESENT'

  check = file_info(files.PA2B)
  if check.exists then begin
     blueSpec_2 = extract_1_indep(files.PA2B, PA = files.PA2)
     redSpec_2  = extract_1_indep(files.PA2R, PA = files.PA2)
     mwrfits, blueSpec_2, outprefix+'_G102_'+string(files.PA2, f = '(I03)')+'.fits', /create
     mwrfits,  redSpec_2, outprefix+'_G141_'+string(files.PA2, f = '(I03)')+'.fits', /create
  endif else $
     print, 'NO EXTRACTION PRESENT'

;  !p.multi = [0,2,0]
;  plot, blueSpec_1.LAMBDA, blueSpec_1.F_OPTI / blueSpec_1.SENSITIVITY, xran = [9000,18000], yran = [0,400]
;  oplot, redSpec_1.LAMBDA,  redSpec_1.F_OPTI /  redSpec_1.SENSITIVITY, col = '777777'x
;
;  plot, blueSpec_2.LAMBDA, blueSpec_2.F_OPTI / blueSpec_2.SENSITIVITY, xran = [9000,18000], yran = [0,400]
;  oplot, redSpec_2.LAMBDA,  redSpec_2.F_OPTI /  redSpec_2.SENSITIVITY, col = '777777'x
;  !p.multi = 0
  
  
end

;;
;;
;;

pro doTaka, field, infile, outdir

  if N_ELEMENTS(infile) eq 0 then $
     infile = '~/Downloads/ps_A2744.txt'

  if N_ELEMENTS(outdir) eq 0 then $
     outdir = '/home/labramson/LocalData/Abramson/A2744specForTaka/'

  if N_ELEMENTS(field) eq 0 then $
     field = 'ABEL2744'
  
;  readcol, infile, ID, RA, DEC, f = 'L, D, D'
  readcol, infile, ID, RA, DEC, f = 'X, L, X, D, D'
  ngals = n_elements(ID)

  for ii = 0, ngals - 1 do $
     get_1ID_spectrum, field, $
                       ra[ii], dec[ii], $
                       outdir = outdir, $
                       outname = string(ID[ii], f = '(I04)')
  
end
;doTaka, 'ABEL0370', '/home/labramson/Downloads/ps_z1.0z6.0_A0370.cat', '/home/labramson/LocalData/Abramson/A0370specForTaka/'
;doTaka, 'ABEL2744', '/home/labramson/Downloads/ps_z1.0z6.0_A2744.cat', '/home/labramson/LocalData/Abramson/A2744specForTaka/'
;doTaka, 'MACS0416', '/home/labramson/Downloads/ps_z1.0z6.0_M0416.cat', '/home/labramson/LocalData/Abramson/M0416specForTaka/'
;doTaka, 'MACS1149', '/home/labramson/Downloads/ps_z1.0z6.0_M1149.cat', '/home/labramson/LocalData/Abramson/M1149specForTaka/'
;doTaka, 'MACS0717', '/home/labramson/Downloads/ps_z1.0z6.0_M0717.cat', '/home/labramson/LocalData/Abramson/M0717specForTaka/'
;doTaka, 'MACS2129', '/home/labramson/Downloads/ps_z1.0z6.0_M2129.cat', '/home/labramson/LocalData/Abramson/M2129specForTaka/'
;doTaka, 'MACS1423', '/home/labramson/Downloads/ps_z1.0z6.0_M1423.cat', '/home/labramson/LocalData/Abramson/M1423specForTaka/'
;doTaka, 'MACS0744', '/home/labramson/Downloads/ps_z1.0z6.0_M0744.cat', '/home/labramson/LocalData/Abramson/M0744specForTaka/'
;doTaka, 'RXJC2248', '/home/labramson/Downloads/ps_z1.0z6.0_A1063.cat', '/home/labramson/LocalData/Abramson/A1063specForTaka/'
;doTaka, 'RXJC1347', '/home/labramson/Downloads/ps_z1.0z6.0_R1347.cat', '/home/labramson/LocalData/Abramson/R1347specForTaka/'
