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
                ZBEST_ROME:       0. , $
                ZSPEC_FLAG_ROME:  0  , $
                ZSPEC_ID_ROME:    0L , $
                PHOTOZ_FLAG_ROME: 0B , $
                B435R:            0. , $
                V606R:            0. , $
                I814R:            0. , $
                Y105R:            0. , $
                J125R:            0. , $
                JH140R:           0. , $
                H160R:            0. , $
                KsR:              0. , $
                CH1R:             0. , $
                CH2R:             0. , $
                errB435R:         0. , $
                errV606R:         0. , $
                errI814R:         0. , $
                errY105R:         0. , $
                errJ125R:         0. , $
                errJH140R:        0. , $
                errH160R:         0. , $
                errKsR:           0. , $
                errCH1R:          0. , $
                errCH2R:          0.   $
                }
     
     savedata = replicate(savedata, ngals)
     for ii = 0, ngals - 1 do begin
        savedata[ii].ID_ROME          = ID[ii]
        savedata[ii].RA               = RA[ii]              
        savedata[ii].DEC              = DEC[ii]                           
        savedata[ii].ZBEST_ROME       = ZBEST[ii]             
        savedata[ii].ZSPEC_FLAG_ROME  = ZSPEC_FLAG[ii]   
        savedata[ii].ZSPEC_ID_ROME    = ZSPEC_ID[ii]       
        savedata[ii].PHOTOZ_FLAG_ROME = PHOTOZ_FLAG[ii] 
        savedata[ii].B435R            = B435[ii]                       
        savedata[ii].V606R            = V606[ii]                       
        savedata[ii].I814R            = I814[ii]                       
        savedata[ii].Y105R            = Y105[ii]                       
        savedata[ii].J125R            = J125[ii]                       
        savedata[ii].JH140R           = JH140[ii]                     
        savedata[ii].H160R            = H160[ii]                       
        savedata[ii].KsR              = Ks[ii]                           
        savedata[ii].CH1R             = CH1[ii]                         
        savedata[ii].CH2R             = CH2[ii]                         
        savedata[ii].errB435R         = errB435[ii]                 
        savedata[ii].errV606R         = errV606[ii]                 
        savedata[ii].errI814R         = errI814[ii]                 
        savedata[ii].errY105R         = errY105[ii]                 
        savedata[ii].errJ125R         = errJ125[ii]                 
        savedata[ii].errJH140R        = errJH140[ii]               
        savedata[ii].errH160R         = errH160[ii]                 
        savedata[ii].errKsR           = errKs[ii]                     
        savedata[ii].errCH1R          = errCH1[ii]                   
        savedata[ii].errCH2R          = errCH2[ii]         
     endfor
     
     mwrfits, savedata, outfits, /create
     print, ''
     print, ' >>> Catalog output to : '+outfits

  endelse

end

pro doa744

readromancats, $
   '/home/labramson/LocalData/Abramson/GLASS/ROMAN_CATALOGS/A2744/A2744_CLUSTER.cat', $
   '/home/labramson/LocalData/Abramson/GLASS/ROMAN_CATALOGS/A2744/A2744_PHOTOZ_FULLSAMPLE.cat', $
   OUTFITS = '$DATA_DIR/GLASS_official_products/ABEL2744/Catalogs/ROMAN_PHOTOZ.fits'

end
