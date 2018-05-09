;; This will combine all of the catalogs in
;; "$DATA_DIR/GLASS_official_products/..." for release to the rest of
;; the team as an "official master catalog."

pro make_taka_coord_fits, infile, outname

  if n_elements(outname) eq 0 then outname = 'TAKA_ID_TO_RADEC.fits'
  
  readcol, infile, $
           ID, RA, DEC, $
           f = 'L,D,D', /silent
  ngals = n_elements(id)

  savedata = {ID_TM:  0L, $
              RA:     0.d, $
              DEC:    0.d}
  savedata = replicate(savedata, ngals)

  for ii = 0, ngals - 1 do begin

     savedata[ii].ID_TM  = id[ii]
     savedata[ii].RA     = ra[ii]
     savedata[ii].DEC    = dec[ii]
     
  endfor

  mwrfits, savedata, outname, /create

  print, ''
  print, ' >>> TAKA ID to COORDINATE file output to : '+outname
  
end

;;
;;
;;

pro make_taka_gal_fits, infile, outname

  if n_elements(outname) eq 0 then outname = 'TAKA_GALFIT_RESULTS.fits'
  
  readcol, infile, $
           ID, $
           X, Y, $
           GALFIT_mag, r_e, $      
           n, q, $                 
           pa, $ ;ID, $            
           GALFIT_BKG, chi2_nu, $  
           err_x, err_y, $        
           err_mag, err_re, $      
           err_n, err_q, $          
           err_pa, $ ;SEX_mag, $      
;           err_SEX_MAG, sky_sext, $
;           rcut, pa_sext $
;           q_sext, rA_sext $
;           r50_sext, r90_sext, $
           N_neig, $
           f = 'L,'+$
           'D,D,'+$
           'F,F,'+$
           'F,F,'+$
           'F,X,'+$
           'F,F,'+$
           'D,D,'+$
           'F,F,'+$
           'F,F,'+$
           'F,X,'+$
           'X,X,'+$
           'X,X,'+$
           'X,X,'+$
           'X,X,'+$
           'I', $
           comment = '#'
  ngals = n_elements(ID)

  savedata = {ID_TM:           0L,$
              X_GALFIT:        0.d, $
              Y_GALFIT:        0.d, $
              MAG_GALFIT:      0., $
              RE_GALFIT:       0., $
              NSERSIC_GALFIT:  0., $
              AXISRAT_GALFIT:  0., $
              PA_GALFIT:       0., $
              BKG_GALFIT:      0., $
              CHI2_NU_GALFIT:  0., $
              EX_GALFIT:       0., $
              EY_GALFIT:       0., $
              EMAG_GALFIT:     0., $
              ERE_GALFIT:      0., $
              ENSERSIC_GALFIT: 0., $
              EAXISRAT_GALFIT: 0.,$
              EPA_GALFIT:      0., $
;              TM_SEX_MAG:    0., $ ;; Will be handled by FAST input file
;              ETM_SEX_MAG:   0., $
;              TM_SEX_SKY:    0., $
              N_NEIGHBORS:     0.}
  savedata = replicate(savedata, ngals)

  for ii = 0, ngals - 1 do begin

     savedata[ii].ID_TM           = ID[ii]
     savedata[ii].X_GALFIT        = X[ii]     
     savedata[ii].Y_GALFIT        = Y[ii] 
     savedata[ii].MAG_GALFIT      = GALFIT_mag[ii]
     savedata[ii].RE_GALFIT       = r_e[ii]
     savedata[ii].NSERSIC_GALFIT  = n[ii]
     savedata[ii].AXISRAT_GALFIT  = q[ii]
     savedata[ii].PA_GALFIT       = pa[ii]
     savedata[ii].BKG_GALFIT      = GALFIT_BKG[ii]
     savedata[ii].CHI2_NU_GALFIT  = chi2_nu[ii]
     savedata[ii].EX_GALFIT       = err_x[ii]
     savedata[ii].EY_GALFIT       = err_y[ii]
     savedata[ii].EMAG_GALFIT     = err_mag[ii]
     savedata[ii].ERE_GALFIT      = err_re[ii] 
     savedata[ii].ENSERSIC_GALFIT = err_n[ii]
     savedata[ii].EAXISRAT_GALFIT = err_q[ii]
     savedata[ii].EPA_GALFIT      = err_pa[ii]
;     savedata[ii].TM_SEX_MAG      = SEX_mag[ii]
;     savedata[ii].ETM_SEX_MAG     = err_SEX_MAG[ii]
;     savedata[ii].TM_SEX_SKY      = sky_sext[ii]
     savedata[ii].N_NEIGHBORS     = N_neig[ii]
     
  endfor

  mwrfits, savedata, outname, /create

  print, ''
  print, ' >>> TAKA GALFIT RESULT file output to : '+outname
  print, ''
  
end

;;
;;
;;

pro make_taka_ez_fits, infile, outname

  if n_elements(outname) eq 0 then outname = 'TAKA_EAZY_RESULTS.fits'
  
  readcol, infile, $ 
           id, $     
           z_spec, $ 
           z_a, $    
           z_m1, $   
           chi_a, $  
           z_p, $    
           chi_p, $  
           z_m2, $   
           odds, $   
           l68, $    
           u68,  $   
           l95, $    
           u95,  $   
           l99, $    
           u99,  $   
           nfilt, $  
           q_z, $    
           z_peak, $ 
           peak_prob, $
           z_mc, $   
           f = 'L,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'I,'+$
           'F,'+$
           'F,'+$
           'F', $
           comment = '#'
  ngals = n_elements(ID)
  
  savedata = {ID_TM:           0L, $
              TM_Z_SPEC:       0., $
              TM_EZ_Z_A:       0., $    
              TM_EZ_Z_M1:      0., $   
              TM_EZ_CHI_A:     0., $  
              TM_EZ_Z_P:       0., $    
              TM_EZ_CHI_P:     0., $  
              TM_EZ_Z_M2:      0., $   
              TM_EZ_ODDS:      0., $   
              TM_EZ_L68:       0., $    
              TM_EZ_U68:       0., $   
              TM_EZ_L95:       0., $    
              TM_EZ_U95:       0., $   
              TM_EZ_L99:       0., $    
              TM_EZ_U99:       0., $   
              TM_EZ_NFILT:     0., $  
              TM_EZ_Q_Z:       0., $    
              TM_EZ_Z_PEAK:    0., $ 
              TM_EZ_PEAK_PROB: 0., $
              TM_EZ_Z_MC:      0.}
  savedata = replicate(savedata, ngals)

  for ii = 0, ngals - 1 do begin

     savedata[ii].ID_TM           = id[ii]                     
     savedata[ii].TM_Z_SPEC       = z_spec[ii]        
     savedata[ii].TM_EZ_Z_A       = z_a[ii]           
     savedata[ii].TM_EZ_Z_M1      = z_m1[ii]         
     savedata[ii].TM_EZ_CHI_A     = chi_a[ii]       
     savedata[ii].TM_EZ_Z_P       = z_p[ii]           
     savedata[ii].TM_EZ_CHI_P     = chi_p[ii]       
     savedata[ii].TM_EZ_Z_M2      = z_m2[ii]         
     savedata[ii].TM_EZ_ODDS      = odds[ii]         
     savedata[ii].TM_EZ_L68       = l68[ii]           
     savedata[ii].TM_EZ_U68       = u68[ii]
     savedata[ii].TM_EZ_L95       = l95[ii]           
     savedata[ii].TM_EZ_U95       = u95[ii]
     savedata[ii].TM_EZ_L99       = l99[ii]           
     savedata[ii].TM_EZ_U99       = u99[ii]
     savedata[ii].TM_EZ_NFILT     = nfilt[ii]       
     savedata[ii].TM_EZ_Q_Z       = q_z[ii]           
     savedata[ii].TM_EZ_Z_PEAK    = z_peak[ii]     
     savedata[ii].TM_EZ_PEAK_PROB = peak_prob[ii]
     savedata[ii].TM_EZ_Z_MC      = z_mc[ii]         
     
  endfor            

  mwrfits, savedata, outname, /create

  print, ''
  print, ' >>> TAKA EAZY RESULT file output to : '+outname
  print, ''
    
end

;;
;;
;;

pro make_taka_fast_fits, infile, outname

  if n_elements(outname) eq 0 then outname = 'TAKA_FAST_RESULTS.fits'
  
  readcol, infile, $
           id, $
           z, $
;           l68_z, $ ; Z IS FIXED
;           u68_z, $
           ltau, $
           l68_ltau, $
           u68_ltau, $
           metal, $
           l68_metal,$
           u68_metal, $
           lage,$
           l68_lage,$
           u68_lage, $
           Av, $
           l68_Av, $
           u68_Av, $
           lmass, $
           l68_lmass, $
           u68_lmass, $
           lsfr, $
           l68_lsfr, $
           u68_lsfr, $
;           lssfr, $    ; CAN CONSTRUCT FROM SFR AND MSTEL
;           l68_lssfr, $
;           u68_lssfr, $
           la2t, $
           l68_la2t, $
           u68_la2t, $
           chi2, $     
           f = 'L,'+$
           'F,'+$
           'X,'+$
           'X,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'X,'+$
           'X,'+$
           'X,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,', $
           comment = '#'
  ngals = n_elements(ID)

  ;; Now read in the full photometry
  photfile = repstr(infile, '.output', '.input') ;; Switch to FAST input
  readcol, photfile, $
;           id, $
           F105, E105, $
           F125, E125, $
           F140, E140, $
           F160, E160, $
           F435, E435, $
           F606, E606, $
           F814, E814, $
           z_spec, $
           f = 'X,'+$
           'F,F,'+$
           'F,F,'+$
           'F,F,'+$
           'F,F,'+$
           'F,F,'+$
           'F,F,'+$
           'F,F,'+$
           'F', $
           comment = '#'
  ngals = n_elements(ID)
  
  savedata = {ID_TM:             0L, $
              TM_F435:           0., $
              TM_E435:           0., $
              TM_F606:           0., $
              TM_E606:           0., $
              TM_F814:           0., $
              TM_E814:           0., $
              TM_F105:           0., $
              TM_E105:           0., $
              TM_F125:           0., $
              TM_E125:           0., $
              TM_F140:           0., $
              TM_E140:           0., $
              TM_F160:           0., $
              TM_E160:           0., $
              TM_Z_USED:         0., $
              TM_FAST_LTAU:      0., $
              TM_FAST_LTAU_L68:  0., $
              TM_FAST_LTAU_U68:  0., $
              TM_FAST_METAL:     0., $
              TM_FAST_METAL_L68: 0., $
              TM_FAST_METAL_U68: 0., $
              TM_FAST_LAGE:      0., $
              TM_FAST_LAGE_L68:  0., $
              TM_FAST_LAGE_U68:  0., $
              TM_FAST_AV:        0., $
              TM_FAST_AV_L68:    0., $
              TM_FAST_AV_U68:    0., $
              TM_FAST_LMASS:     0., $
              TM_FAST_LMASS_L68: 0., $
              TM_FAST_LMASS_U68: 0., $
              TM_FAST_LSFR:      0., $
              TM_FAST_LSFR_L68:  0., $
              TM_FAST_LSFR_U68:  0., $
              TM_FAST_LA2T:      0., $
              TM_FAST_LA2T_L68:  0., $
              TM_FAST_LA2T_U68:  0., $
              TM_FAST_CHI2:      0.}
  savedata = replicate(savedata, ngals)

  for ii = 0, ngals - 1 do begin

     savedata[ii].ID_TM             = id[ii]                    
     savedata[ii].TM_F435           = F435[ii]
     savedata[ii].TM_E435           = E435[ii]
     savedata[ii].TM_F606           = F606[ii]
     savedata[ii].TM_E606           = E606[ii]
     savedata[ii].TM_F814           = F814[ii]   
     savedata[ii].TM_E814           = E814[ii]
     savedata[ii].TM_F105           = F105[ii]
     savedata[ii].TM_E105           = E105[ii]
     savedata[ii].TM_F125           = F125[ii]
     savedata[ii].TM_E125           = E125[ii]
     savedata[ii].TM_F140           = F140[ii]
     savedata[ii].TM_E140           = E140[ii]
     savedata[ii].TM_F160           = F160[ii]
     savedata[ii].TM_E160           = E160[ii]
     savedata[ii].TM_Z_USED         = z[ii]                 
     savedata[ii].TM_FAST_LTAU      = ltau[ii]        
     savedata[ii].TM_FAST_LTAU_L68  = l68_ltau[ii]    
     savedata[ii].TM_FAST_LTAU_U68  = u68_ltau[ii]   
     savedata[ii].TM_FAST_METAL     = metal[ii]         
     savedata[ii].TM_FAST_METAL_L68 = l68_metal[ii]
     savedata[ii].TM_FAST_METAL_U68 = u68_metal[ii]  
     savedata[ii].TM_FAST_LAGE      = lage[ii]            
     savedata[ii].TM_FAST_LAGE_L68  = l68_lage[ii]   
     savedata[ii].TM_FAST_LAGE_U68  = u68_lage[ii]   
     savedata[ii].TM_FAST_AV        = Av[ii]               
     savedata[ii].TM_FAST_AV_L68    = l68_Av[ii]       
     savedata[ii].TM_FAST_AV_U68    = u68_Av[ii]       
     savedata[ii].TM_FAST_LMASS     = lmass[ii]         
     savedata[ii].TM_FAST_LMASS_L68 = l68_lmass[ii]  
     savedata[ii].TM_FAST_LMASS_U68 = u68_lmass[ii]  
     savedata[ii].TM_FAST_LSFR      = lsfr[ii]           
     savedata[ii].TM_FAST_LSFR_L68  = l68_lsfr[ii]   
     savedata[ii].TM_FAST_LSFR_U68  = u68_lsfr[ii]   
     savedata[ii].TM_FAST_LA2T      = la2t[ii]        
     savedata[ii].TM_FAST_LA2T_L68  = l68_la2t[ii]
     savedata[ii].TM_FAST_LA2T_U68  = u68_la2t[ii] 
     savedata[ii].TM_FAST_CHI2      = chi2[ii]     
                                       
  endfor                               

  mwrfits, savedata, outname, /create

  print, ''
  print, ' >>> TAKA FAST RESULT file output to : '+outname
  print, ''
  
end                                    
                                       
;;
;;
;;

pro join2onid, masterCat, slaveCat, $
               OUTNAME = outname

  if NOT keyword_set(OUTNAME) then outname = 'tmpJoin.fits'

  master = mrdfits(masterCat, 1) ;; THIS IS THE CATALOG LENGTH TO BE MATCHED
  n1 = n_elements(master)

  slave = mrdfits(slaveCat, 1)   ;; THIS IS THE CATALOG THAT MIGHT NEED TO BE EXPANDED
  n2 = n_elements(slave)
  
  newtags = tag_names(slave)
  ntags = n_elements(newtags)

  minds = lonarr(n2)
  sinds = lonarr(n2)
  for ii = 0, n2 - 1 do begin
     hits = where(master.ID_TM eq slave[ii].ID_TM, nhits)
     if nhits eq 1 then begin
        minds[ii] = hits[0]
        sinds[ii] = ii
     endif
  endfor
  
  striptag = 'ID_TM'
  if n_elements(striptag) gt 0 then begin
     check = []    
     for ii = 0, n_elements(striptag) - 1 do $
        check = [check, where(newtags eq striptag[ii])]
  endif else $
     check = !VALUES.F_NAN
  
  tagsdone = 0
  for ii = 0, ntags - 1 do begin

     if ~total(check eq ii) then begin
     
        tmpdata = slave.(ii)
        type    = size(tmpdata, /tname)
        case type of
           'STRING' : tmpVal = 'string(-99)'
           'FLOAT'  : tmpVal = '-99.'
           'LONG'   : tmpVal = '-99L'
           'INT'    : tmpVal = '-99'
           'DOUBLE' : tmpVal = '-99.d'
        endcase
        
        if tagsdone eq 0 then begin
           com = 'tmpstruct = {'+newtags[ii]+': '+string(tmpVal)+'}'
           foo = execute(com)
           com = 'tmpstruct = replicate(tmpstruct, n1)'
           foo = execute(com)           
        endif else begin
           com = 'tmpstruct1 = {'+newtags[ii]+': '+string(tmpVal)+'}'
           foo = execute(com)
           com = 'tmpstruct1 = replicate(tmpstruct1, n1)'
           foo = execute(com)
           tmpstruct = struct_addtags(tmpstruct, tmpstruct1)
        endelse
        
        tmpstruct[minds].(tagsdone) = tmpdata[sinds]

        tagsdone++
        
     endif

  endfor

  master = struct_addtags(master, tmpstruct)

  mwrfits, master, outname, /create

  print, ''
  print, ' >>> '+slaveCat+' JOINED TO '+masterCat+' AND DUMPED TO '+outname
  print, ''
  
end

;;
;;
;;

pro combine_taka_cats, fieldName, $
                       COORDS_NAME = coords_name, $
                       GALFIT_NAME = galfit_name, $
                       EAZY_NAME   = eazy_name, $
                       FAST_NAME   = fast_name, $
                       OUTFILE = outfile

;  dataDir   = "$DATA_DIR/GLASS_official_products/"+fieldName+"/Catalogs/"

;  if NOT keyword_set(COORDS_NAME) then coords_name = "TAKA_ID_TO_RADEC.fits"
;  if NOT keyword_set(GALFIT_NAME) then galfit_name = "TAKA_GALFIT_RESULTS.fits"
;  if NOT keyword_set(EAZY_NAME)   then eazy_name   = "TAKA_EAZY_RESULTS.fits"
;  if NOT keyword_set(FAST_NAME)   then fast_name   = "TAKA_FAST_RESULTS.fits"

;  coordfile = dataDir+coords_name
;  galFile   = dataDir+galfit_name
;  ezFIle    = dataDir+eazy_name
;  fastFIle  = dataDir+fast_name

  coordfile = coords_name
  galFile   = galfit_name
  ezFile    = eazy_name
  fastFile  = fast_name
  
  join2onid, coordfile, galFile
  join2onid, 'tmpJoin.fits', ezFIle
  join2onid, 'tmpJoin.fits', fastFIle, OUTNAME = outfile
  
end

;;
;;
;;

;; Read in both the roman photmoetry file AND the derived photo-zs,
;; merge them, and spit-out a FITS binary table

pro readRomanCats, photometryCat, redshiftCat, $
                   OUTFITS = outfits

  ;; Read the photometry catalog
  readcol, photometryCat, $
           ID, RA, DEC, $
           B435, V606, I814, $
           Y105, J125, JH140,$
           H160, Ks, CH1,$
           CH2 ,$
           errB435, errV606, errI814,$
           errY105, errJ125, errJH140,$
           errH160, errKs, errCH1, errCH2, $
           f = 'L,D,D,'+$
           'F,F,F,'+$
           'F,F,F,'+$
           'F,F,F,'+$
           'F,'+$
           'F,F,F,'+$
           'F,F,F,'+$
           'F,F,F', comment = '#'
  ngals = n_elements(ID)
  
  ;; Read the photo-z catalog
  readcol, redshiftCat, $
           zID, ZBEST, $
           ZSPEC_FLAG, ZSPEC_ID, $
           PHOTOZ_FLAG, $
           f = 'L,'+$
           'F,F,'+$
           'F,F,'+$
           'F', comment = '#', /quick

  ;; Check that they have the same number of objects.
  ;; If so, dump to output.
  if total(ID - zID) ne 0 then begin
     print, ''
     print, ' !!! CATALOGS DO NOT ALIGN --> ABORTING !!!'
     stop
  endif else begin
     print, ''
     print, ' >>> Catalogs look aligned; dumping to output ... '

     savedata = {$
                ID_ROME:          0L , $
                RA:               0.d, $
                DEC:              0.d, $
                ROME_ZBEST:       0. , $
                ROME_ZSPEC_FLAG:  0  , $
                ROME_ZSPEC_ID:    0L , $
                ROME_PHOTOZ_FLAG: 0B , $
                ROME_F435:        0. , $
                ROME_F606:        0. , $
                ROME_F814:        0. , $
                ROME_F105:        0. , $
                ROME_F125:        0. , $
                ROME_F140:        0. , $
                ROME_F160:        0. , $
                ROME_Ks  :        0. , $
                ROME_F36 :        0. , $
                ROME_F45 :        0. , $
                ROME_E435:        0. , $
                ROME_E606:        0. , $
                ROME_E814:        0. , $
                ROME_E105:        0. , $
                ROME_E125:        0. , $
                ROME_E140:        0. , $
                ROME_E160:        0. , $
                ROME_EKs :        0. , $
                ROME_E36 :        0. , $
                ROME_E45 :        0.   $
                }
     
     savedata = replicate(savedata, ngals)
     for ii = 0, ngals - 1 do begin
        savedata[ii].ID_ROME          = ID[ii]
        savedata[ii].RA               = RA[ii]              
        savedata[ii].DEC              = DEC[ii]                           
        savedata[ii].ROME_ZBEST       = ZBEST[ii]             
        savedata[ii].ROME_ZSPEC_FLAG  = ZSPEC_FLAG[ii]   
        savedata[ii].ROME_ZSPEC_ID    = ZSPEC_ID[ii]       
        savedata[ii].ROME_PHOTOZ_FLAG = PHOTOZ_FLAG[ii] 
        savedata[ii].ROME_F435        = B435[ii]                       
        savedata[ii].ROME_F606        = V606[ii]                       
        savedata[ii].ROME_F814        = I814[ii]                       
        savedata[ii].ROME_F105        = Y105[ii]                       
        savedata[ii].ROME_F125        = J125[ii]                       
        savedata[ii].ROME_F140        = JH140[ii]                     
        savedata[ii].ROME_F160        = H160[ii]                       
        savedata[ii].ROME_Ks          = Ks[ii]                           
        savedata[ii].ROME_F36         = CH1[ii]                         
        savedata[ii].ROME_F45         = CH2[ii]                         
        savedata[ii].ROME_E435        = errB435[ii]                 
        savedata[ii].ROME_E606        = errV606[ii]                 
        savedata[ii].ROME_E814        = errI814[ii]                 
        savedata[ii].ROME_E105        = errY105[ii]                 
        savedata[ii].ROME_E125        = errJ125[ii]                 
        savedata[ii].ROME_E140        = errJH140[ii]               
        savedata[ii].ROME_E160        = errH160[ii]                 
        savedata[ii].ROME_EKs         = errKs[ii]                     
        savedata[ii].ROME_E36         = errCH1[ii]                   
        savedata[ii].ROME_E45         = errCH2[ii]         
     endfor
     
     mwrfits, savedata, outfits, /create
     print, ''
     print, ' >>> Catalog output to : '+outfits

  endelse

end

;;
;;
;;

;; Convert a GLASS Master ASCII Catalog to a FITS table

pro readGlassList, incat, outname

;NUMBER
;X_IMAGE Y_IMAGE
;X_WORLD Y_WORLD
;A_IMAGE B_IMAGE THETA_IMAGE
;A_WORLD B_WORLD THETA_WORLD
;FLUX_APER FLUX_APER2 FLUX_APER3
;FLUXERR_APER FLUXERR_APER2 FLUXERR_APER3
;MAG_APER MAG_APER2 MAG_APER3
;MAGERR_APER MAGERR_APER2 MAGERR_APER3
;FLUX_AUTO FLUXERR_AUTO
;MAG_AUTO MAGERR_AUTO

;KRON_RADIUS PETRO_RADIUS
;BACKGROUND THRESHOLD
;XWIN_IMAGE YWIN_IMAGE
;AWIN_IMAGE BWIN_IMAGE THETAWIN_IMAGE
;MU_THRESHOLD FLAGS FWHM_IMAGE
;FLUX_RADIUS FLUX_RADIUS2 CLASS_STAR X_FLT Y_FLT
  
  ;; Do fluxes (mainly)
  print, ''
  print, 'Reading MASTER GLASS catalog ... '
  print, ''
  readcol, incat, $
           NUMBER, $
           X_IMAGE, Y_IMAGE, $
           X_WORLD, Y_WORLD, $
           A_IMAGE, B_IMAGE, THETA_IMAGE, $
           A_WORLD, B_WORLD, THETA_WORLD, $
           FLUX_APER, FLUX_APER2, FLUX_APER3, $
           FLUXERR_APER, FLUXERR_APER2, FLUXERR_APER3, $
           MAG_APER, MAG_APER2, MAG_APER3, $
           MAGERR_APER, MAGERR_APER2, MAGERR_APER3, $
           FLUX_AUTO, FLUXERR_AUTO, $
           MAG_AUTO, MAGERR_AUTO, $
           f = $
           'F,'+$
           'F,F,'+$
           'D,D,'+$
           'F,F,F,'+$
           'F,F,F,'+$
           'F,F,F,'+$
           'F,F,F,'+$
           'F,F,F,'+$
           'F,F,F,'+$
           'F,F,'+$
           'F,F'
  number = long(number)
  
  ;; Do the rest
  readcol, incat, $
           KRON_RADIUS, PETRO_RADIUS, $
           BACKGROUND, THRESHOLD, $
           XWIN_IMAGE, YWIN_IMAGE, $
           AWIN_IMAGE, BWIN_IMAGE, THETAWIN_IMAGE, $
           MU_THRESHOLD, FLAGS, FWHM_IMAGE, $
           FLUX_RADIUS, FLUX_RADIUS2, CLASS_STAR, X_FLT, Y_FLT, $
           f = 'X,'+$
           'X,X,'+$
           'X,X,'+$
           'X,X,X,'+$
           'X,X,X,'+$
           'X,X,X,'+$
           'X,X,X,'+$
           'X,X,X,'+$
           'X,X,X,'+$
           'X,X,'+$
           'X,X'+$  
           'F,F,'+$
           'F,F,'+$
           'F,F,'+$
           'F,F,F,'+$
           'F,L64,F,'+$
           'F,F,F,F,F'

  savedata = {$
             NUMBER:        0L  ,$
             X_IMAGE:       0.  ,$
             Y_IMAGE:       0.  ,$
             X_WORLD:       0.d ,$
             Y_WORLD:       0.d ,$
             A_IMAGE:       0.  ,$
             B_IMAGE:       0.  ,$
             THETA_IMAGE:   0.  ,$
             A_WORLD:       0.  ,$
             B_WORLD:       0.  ,$
             THETA_WORLD:   0.  ,$
             FLUX_APER:     0.  ,$
             FLUX_APER2:    0.  ,$
             FLUX_APER3:    0.  ,$
             FLUXERR_APER:  0.  ,$
             FLUXERR_APER2: 0.  ,$
             FLUXERR_APER3: 0.  ,$
             MAG_APER:      0.  ,$
             MAG_APER2:     0.  ,$
             MAG_APER3:     0.  ,$
             MAGERR_APER:   0.  ,$
             MAGERR_APER2:  0.  ,$
             MAGERR_APER3:  0.  ,$
             FLUX_AUTO:     0.  ,$
             FLUXERR_AUTO:  0.  ,$
             MAG_AUTO:      0.  ,$
             MAGERR_AUTO:   0.  ,$
             KRON_RADIUS:   0.  ,$
             PETRO_RADIUS:  0.  ,$
             BACKGROUND:    0.  ,$
             THRESHOLD:     0.  ,$
             XWIN_IMAGE:    0.  ,$
             YWIN_IMAGE:    0.  ,$
             AWIN_IMAGE:    0.  ,$
             BWIN_IMAGE:    0.  ,$
             THETAWIN_IMAGE:0.  ,$
             MU_THRESHOLD:  0.  ,$
             FLAGS:  long64(0)  ,$
             FWHM_IMAGE:    0.  ,$
             FLUX_RADIUS:   0.  ,$
             FLUX_RADIUS2:  0.  ,$
             CLASS_STAR:    0.  ,$
             X_FLT:         0.  ,$
             Y_FLT:         0.   $
             }
  savedata = replicate(savedata, n_elements(X_WORLD))

  for ii = 0, n_elements(savedata) - 1 do begin
     savedata[ii].NUMBER        = NUMBER[ii]
     savedata[ii].X_IMAGE       = X_IMAGE[ii]      
     savedata[ii].Y_IMAGE       = Y_IMAGE[ii]  
     savedata[ii].X_WORLD       = X_WORLD[ii]            
     savedata[ii].Y_WORLD       = Y_WORLD[ii]
     savedata[ii].A_IMAGE       = A_IMAGE[ii]           
     savedata[ii].B_IMAGE       = B_IMAGE[ii]           
     savedata[ii].THETA_IMAGE   = THETA_IMAGE[ii]   
     savedata[ii].A_WORLD       = A_WORLD[ii]             
     savedata[ii].B_WORLD       = B_WORLD[ii]             
     savedata[ii].THETA_WORLD   = THETA_WORLD[ii]    
     savedata[ii].FLUX_APER     = FLUX_APER[ii]         
     savedata[ii].FLUX_APER2    = FLUX_APER2[ii]    
     savedata[ii].FLUX_APER3    = FLUX_APER3[ii]       
     savedata[ii].FLUXERR_APER  = FLUXERR_APER[ii]  
     savedata[ii].FLUXERR_APER2 = FLUXERR_APER2[ii]
     savedata[ii].FLUXERR_APER3 = FLUXERR_APER3[ii] 
     savedata[ii].MAG_APER      = MAG_APER[ii]           
     savedata[ii].MAG_APER2     = MAG_APER2[ii]         
     savedata[ii].MAG_APER3     = MAG_APER3[ii]         
     savedata[ii].MAGERR_APER   = MAGERR_APER[ii]     
     savedata[ii].MAGERR_APER2  = MAGERR_APER2[ii]   
     savedata[ii].MAGERR_APER3  = MAGERR_APER3[ii]  
     savedata[ii].FLUX_AUTO     = FLUX_AUTO[ii]         
     savedata[ii].FLUXERR_AUTO  = FLUXERR_AUTO[ii]   
     savedata[ii].MAG_AUTO      = MAG_AUTO[ii]           
     savedata[ii].MAGERR_AUTO   = MAGERR_AUTO[ii]     
     savedata[ii].KRON_RADIUS   = KRON_RADIUS[ii]     
     savedata[ii].PETRO_RADIUS  = PETRO_RADIUS[ii]   
     savedata[ii].BACKGROUND    = BACKGROUND[ii]    
     savedata[ii].THRESHOLD     = THRESHOLD[ii]       
     savedata[ii].XWIN_IMAGE    = XWIN_IMAGE[ii]      
     savedata[ii].YWIN_IMAGE    = YWIN_IMAGE[ii]    
     savedata[ii].AWIN_IMAGE    = AWIN_IMAGE[ii]        
     savedata[ii].BWIN_IMAGE    = BWIN_IMAGE[ii]        
     savedata[ii].THETAWIN_IMAGE= THETAWIN_IMAGE[ii]
     savedata[ii].MU_THRESHOLD  = MU_THRESHOLD[ii] 
     savedata[ii].FLAGS         = FLAGS[ii]  
     savedata[ii].FWHM_IMAGE    = FWHM_IMAGE[ii]        
     savedata[ii].FLUX_RADIUS   = FLUX_RADIUS[ii]
     savedata[ii].FLUX_RADIUS2  = FLUX_RADIUS2[ii]    
     savedata[ii].CLASS_STAR    = CLASS_STAR[ii]        
     savedata[ii].X_FLT         = X_FLT[ii]                  
     savedata[ii].Y_FLT         = Y_FLT[ii] 
  endfor

  print, 'Saving to FITS ... '
  print, ''
  mwrfits, savedata, outname, /create
  print, ' >>> New FITS catalog output to : '+outname
  print, ''
  
end

;;
;;
;;

pro readRedshift, catIn, fitsOut

  ;;ID   RA    DEC
  ;;redshift  redshift_quality   multiple_redshift_solutions   

  readcol, catIn, $
           glassID, ra, dec, $
           z, zq, mzs, $
           f = 'A,D,D,F,F,F'
  ngals = n_elements(glassID)

  savedata = {GLASS_ID         : 0L, $
              RA               : 0.d, $
              DEC              : 0.d, $
              Z_GLASS          : 0., $
              Z_GLASS_QUAL     : 0., $
              Z_GLASS_MULTI_SOL: 0.}
  savedata = replicate(savedata, ngals)

  for ii = 0, ngals - 1 do begin

     savedata[ii].GLASS_ID          = long(glassID[ii])
     savedata[ii].RA                = ra[ii]
     savedata[ii].DEC               = dec[ii]
     savedata[ii].Z_GLASS           = z[ii]
     savedata[ii].Z_GLASS_QUAL      = zq[ii]
     savedata[ii].Z_GLASS_MULTI_SOL = mzs[ii]
     
  endfor

  mwrfits, savedata, fitsOut, /create
  
end

;;
;;
;;

pro readGiG, gigIn, fitsOut

;; ID RA DEC
;; CLASH_ID CLASH_DISTANCE MAG_SELECTION MAGERR_SELECTION
;; CONTAM_G102_20 DEFECT_G102_20 CONTAM_DEFECT_G102_20
;; CONTAM_G102_280 CONTAM_DEFECT_G102_280 DEFECT_G102_280
;; CONTAM_G141_20 CONTAM_DEFECT_G141_20 DEFECT_G141_20
;; CONTAM_G141_280    CONTAM_DEFECT_G141_280 DEFECT_G141_280
;; STAR

  readcol, gigIn, $
           glassID, ra, dec, $
           clashID, clash_dpos, magsel, emagsel, $
           contam_102_p1, defect_102_p1, contam_defect_102_p1, $
           contam_102_p2, defect_102_p2, contam_defect_102_p2, $
           contam_141_p1, defect_141_p1, contam_defect_141_p1, $
           contam_141_p2, defect_141_p2, contam_defect_141_p2, $
           starflag, $
           f = 'A,D,D,'+$
           'A,D,F,F,'+$
           'F,F,F,'+$
           'F,F,F,'+$
           'F,F,F,'+$
           'F,F,F,'+$
           'F', comment = '#'
  ngals = n_elements(glassID)
  
  spawn, 'cat '+gigIn+' | grep "PA1=" > tmp.txt'
  readcol, 'tmp.txt', str, f = 'A', $
           comment = '', delimit = '_'

  x1 = strpos(str, "PA1=")
  x2 = strpos(str, "PA2=")
  pa1 = strmid(str, x1+4, 3)
  pa2 = strmid(str, x2+4, 3)

  savedata = {GLASS_ID: 0L, RA: 0.d, DEC: 0.d, PA1: 0., PA2: 0., $
              CLASH_ID: 0L, CLASH_DPOS: 0.d, $
              MAG_SELECT: 0., E_MAG_SELECT:0., $
              CONTAM_G102_PA1: 0., CONTAM_G102_PA2: 0., $
              DEFECT_G102_PA1: 0., DEFECT_G102_PA2: 0., $
              CONDEF_G102_PA1: 0., CONDEF_G102_PA2: 0., $
              CONTAM_G141_PA1: 0., CONTAM_G141_PA2: 0., $
              DEFECT_G141_PA1: 0., DEFECT_G141_PA2: 0., $
              CONDEF_G141_PA1: 0., CONDEF_G141_PA2: 0.}
  savedata = replicate(savedata, ngals)

  for ii = 0, ngals - 1 do begin

     savedata[ii].GLASS_ID        = long(glassId[ii])
     savedata[ii].RA              = ra[ii]
     savedata[ii].DEC             = dec[ii]
     savedata[ii].PA1             = float(pa1)
     savedata[ii].PA2             = float(pa2)
     savedata[ii].CLASH_ID        = clashID[ii]
     savedata[ii].CLASH_DPOS      = clash_dpos[ii]
     savedata[ii].MAG_SELECT      = magsel[ii]
     savedata[ii].E_MAG_SELECT    = emagsel[ii]
     savedata[ii].CONTAM_G102_PA1 = contam_102_p1[ii]
     savedata[ii].CONTAM_G102_PA2 = contam_102_p2[ii]
     savedata[ii].DEFECT_G102_PA1 = defect_102_p1[ii]
     savedata[ii].DEFECT_G102_PA2 = defect_102_p2[ii]
     savedata[ii].CONDEF_G102_PA1 = contam_defect_102_p1[ii]
     savedata[ii].CONDEF_G102_PA2 = contam_defect_102_p2[ii]
     savedata[ii].CONTAM_G141_PA1 = contam_141_p1[ii]
     savedata[ii].CONTAM_G141_PA2 = contam_141_p2[ii]
     savedata[ii].DEFECT_G141_PA1 = defect_141_p1[ii]
     savedata[ii].DEFECT_G141_PA2 = defect_141_p2[ii]
     savedata[ii].CONDEF_G141_PA1 = contam_defect_141_p1[ii]
     savedata[ii].CONDEF_G141_PA2 = contam_defect_141_p2[ii]
     
  endfor

  mwrfits, savedata, fitsOut, /create
  
end

;;
;;
;;

pro join2cats, masterIn, otherCatIn, output, $
               STRIPTAG = striptag, $
               POINTING = pointing, $
               DPOSNAME = dposname
  
  master = mrdfits(masterIn, 1)
  n1 = n_elements(master)
  mra  = master.X_WORLD
  mdec = master.Y_WORLD

  slave  = mrdfits(otherCatIn, 1)
  n2 = n_elements(slave)

  newtags = tag_names(slave)
  ntags = n_elements(newtags)
  if total(newtags eq 'ALPHA_J2000') eq 0 then begin
     sra  = slave.RA
     sdec = slave.DEC
  endif else begin
     sra  = slave.ALPHA_J2000  ;; Kuang puts his photozs this way
     sdec = slave.DELTA_J2000  
  endelse
     
  spherematch, mra, mdec, sra, sdec, 0.2d/3600, $
               minds, sinds, len
  if n_elements(striptag) gt 0 then begin
     check = []    
     for ii = 0, n_elements(striptag) - 1 do $
        check = [check, where(newtags eq striptag[ii])]
;        slave = mod_struct(slave, striptag[ii], /delete)
  endif else check = !VALUES.F_NAN
;  newtags = tag_names(slave)
;  ntags = n_elements(newtags)

  tagsdone = 0
  for ii = 0, ntags - 1 do begin

     if ~total(check eq ii) then begin
     
        tmpdata = slave.(ii)
        type    = size(tmpdata, /tname)
        case type of
           'STRING' : tmpVal = 'string(-99)'
           'FLOAT'  : tmpVal = '-99.'
           'LONG'   : tmpVal = '-99L'
           'INT'    : tmpVal = '-99'
           'BYTE'   : tmpVal = '-99'
           'DOUBLE' : tmpVal = '-99.d'
        endcase
        
        if tagsdone eq 0 then begin
           com = 'tmpstruct = {'+newtags[ii]+': '+string(tmpVal)+'}'
           foo = execute(com)
           com = 'tmpstruct = replicate(tmpstruct, n1)'
           foo = execute(com)           
        endif else begin
           com = 'tmpstruct1 = {'+newtags[ii]+': '+string(tmpVal)+'}'
           foo = execute(com)
           com = 'tmpstruct1 = replicate(tmpstruct1, n1)'
           foo = execute(com)
           tmpstruct = struct_addtags(tmpstruct, tmpstruct1)
        endelse
        
        tmpstruct[minds].(tagsdone) = tmpdata[sinds]

        tagsdone++
        
     endif

  endfor

  ;; Add a column for source offsets
  if KEYWORD_SET(DPOSNAME) then begin
     com = 'tmpstruct1 = {'+dposname+': -99.d}'
     foo = execute(com)
     com = 'tmpstruct1 = replicate(tmpstruct1, n1)'
     foo = execute(com)
     tmpstruct = struct_addtags(tmpstruct, tmpstruct1)
     tmpstruct[minds].(tagsdone) = len[sinds] * 3600
  endif
  
  master = struct_addtags(master, tmpstruct)

  if keyword_set(POINTING) then begin
     tmp = {POINTING: pointing}
     tmp = replicate(tmp, n_elements(master))
     master = struct_addtags(master, tmp)
  endif

  mwrfits, master, output, /create

end

;;
;;
;;

pro combine_for_release, fieldName, outputName, $
                         POINTING = pointing, $
                         OUTDIR = outdir   
  case fieldName of
     'ABEL0370': begin
        rfname = 'A0370'
        sfname = 'a370'
     end
     'ABEL2744': begin
        rfname = 'A2744'
        sfname = 'a2744'
     end
     'MACS0416': begin
        rfname = 'M0416'
        sfname = 'macs0416'
     end
     'MACS0717': begin
        rfname = 'M0717'
        sfname = 'macs0717'
     end
     'MACS1149': begin
        rfname = 'M1149'
        sfname = 'macs1149'
     end
  endcase

  
  ;; Set up file structure
  ;; ---------------------
  
  dataDir = "$DATA_DIR/GLASS_official_products/"+fieldName+"/Catalogs/"
  if NOT keyword_set(OUTDIR) then outdir = dataDir
  if NOT keyword_set(POINTING) then pointing = fieldName+'_CLS'  ;; Default to cluster pointing per Morishita+16a
  if n_elements(outputName) eq 0 then outputName = outDir+pointing+'_allData.fits'

  ;; GLASS official things plus TM's results
  glassMaster  = dataDir+'hlsp_glass_hst_wfc3_'+sfname+'-fullfov-pa999_ir_v001_glassmaster.cat'
  glassGig     = dataDir+'hlsp_glass_hst_wfc3_'+sfname+'-fullfov-pa999_ir_v001_glassgigcatalog.txt'
  glassReds    = dataDir+'hlsp_glass_hst_wfc3_'+sfname+'-fullfov-pa999_ir_v001_redshiftcatalog.txt'  
  takaIDFile   = dataDir+'TAKA_IDtoRADEC.cat'
  takaGalFile  = dataDir+'TAKA_GALFIT.output'
  takaEZfile   = dataDir+'TAKA_EAZY.output'
  takaFASTfile = dataDir+'TAKA_FAST.output'

  ;; Stuff from gitHub; Roman catalogs + Kuang's versions thereof
  romanZdir    = "/home/labramson/LocalData/Abramson/GLASS/ROMAN_CATALOGS/"+rfname
  romanPhotCat = romanZdir+'/'+rfname+'_CLUSTER.cat'
  romanRedCat  = romanZdir+'/'+rfname+'_PHOTOZ_FULLSAMPLE.cat'
  

  ;; Manipulate files
  ;; ----------------
  
  ;; Make Roman photo-z catalog
  readromancats, romanPhotCat, romanRedCat, $
                 outfits = outdir+'ROMAN_PHOTOZ_CATALOG.fits'

  ;; Combine Takahiro's stuff
  make_taka_coord_fits, takaIDfile, outdir+'TM_ID_TO_RADEC.fits'
  make_taka_gal_fits, takaGalFile, outdir+'TM_GALFIT_RESULTS.fits'
  make_taka_ez_fits, takaEZfile, outdir+'TM_EAZY_RESULTS.fits'
  make_taka_fast_fits, takaFASTfile, outdir+'TM_FAST_RESULTS.fits'
  combine_taka_cats, fieldName, $
                     COORDS_NAME = outdir+'TM_ID_TO_RADEC.fits', $
                     GALFIT_NAME = outdir+'TM_GALFIT_RESULTS.fits', $
                     EAZY_NAME   = outdir+'TM_EAZY_RESULTS.fits', $
                     FAST_NAME   = outdir+'TM_FAST_RESULTS.fits', $
                     OUTFILE     = outdir+'ALL_TM_RESULTS.fits'

  ;; Combine Glass official things
  readglasslist, glassMaster, outdir+'GLASS_SRCLIST.fits'    
  readredshift , glassReds  , outdir+'GLASS_REDSHIFTS.fits'  
  readgig      , glassGig   , outdir+'GLASS_GIG_RESULTS.fits'

  master    = outdir+'GLASS_SRCLIST.fits'
  gigCat    = outdir+'GLASS_GIG_RESULTS.fits'
  redshifts = outdir+'GLASS_REDSHIFTS.fits'  
  print, ''
  print, 'APPENDING GiG RESULTS TO MASTER GLASS SOURCE CATALOG ... '
  print, ''
  join2cats, master, gigCat, 'tmp1.fits', $
             striptag = 'GLASS_ID'
  print, ''
  print, 'APPENDING GiGZ RESULTS TO GLASS SOURCE + GiG CATALOG ... '
  print, ''
  join2cats, 'tmp1.fits', redshifts, 'tmp2.fits', $
             striptag = ['GLASS_ID', 'RA', 'DEC']


  ;; Combine the independent pieces
  ;; ------------------------------
  
  print, ''
  print, 'APPENDING GALFIT+MASS+TM_PHOTO_Z RESULTS TO GLASS SOURCE + GiG + GiGz CATALOG ... '
  print, ''
  join2cats, 'tmp2.fits', outdir+'ALL_TM_RESULTS.fits', 'tmp3.fits', $
             striptag = ['RA', 'DEC'], $
             dposname = 'TM2GLASS_DPOS'

  ;; Add Roman/Kuang Photo-zs
  print, ''
  print, 'APPENDING KUANG/ROMAN PHOTO-ZS TO GLASS SOURCE + GiG + GiGz CATALOG ... '
  print, ''
  join2cats, 'tmp3.fits', outdir+'ROMAN_PHOTOZ_CATALOG.fits', $
             outputName, $
             striptag = ['RA', 'DEC'], $
             dposname = 'ROMAN2GLASS_DPOS', $
             POINTING = pointing

  
  ;; Define a log file
  ;; -----------------
  close, 1
  openw, 1, repstr(outputName, '.fits', '.log'), width = 256
  printf, 1, 'File summary for master catalog : '+outputName
  printf, 1, '----------------------------------------------'
  printf, 1, ' '
  printf, 1, 'Field    : '+fieldName
  printf, 1, 'Pointing : '+pointing
  printf, 1, ''
  printf, 1, ' ** Syntax is TXT_FILE >> CORRESPONDING_FITS_TABLE ** '
  printf, 1, ''
  printf, 1, 'GLASS master catalog   : '+glassMaster+' >> '+outdir+'GLASS_SRCLIST.fits'
  printf, 1, 'GLASS GiG catalog      : '+glassGig+' >> '+outdir+'GLASS_GIG_RESULTS.fits'
  printf, 1, 'GLASS redshift catalog : '+glassReds+' >> '+outdir+'GLASS_REDSHIFTS.fits'
  printf, 1, ''
  printf, 1, 'GALFIT catalog         : '+takaGalFile+' >> '+outdir+'TM_GALFIT_RESULTS.fits'
  printf, 1, 'TM_PHOTZ catalog       : '+takaEZFile+' >> '+outdir+'TM_EAZY_RESULTS.fits'
  printf, 1, 'TM_FAST catalog        : '+takaFASTFile+' >> '+outdir+'TM_FAST_RESULTS.fits'
  printf, 1, 'TM ID to RADEC list    : '+takaIDfile+' >> '+outdir+'TM_ID_TO_RADEC.fits'
  printf, 1, ''
  printf, 1, 'ROMAN PHOTZ catalog    : '+outdir+'ROMAN_PHOTOZ_CATALOG.fits' ;; this may still change as of 20160711
  printf, 1, ''
  printf, 1, 'Catalog constructed    : '+systime()
  close, 1
  
  print, ''
  print, ''
  print, ''
  print, ' * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * '
  print, ' * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * '
  print, ' MASTER CATALOG FOR POINTING '+pointing+' OUTPUT TO : '+outputName
  print, ' * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * '
  print, ' * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * '
  print, ''
  print, ''
  
end
;; combine_for_release, 'ABEL2744'
