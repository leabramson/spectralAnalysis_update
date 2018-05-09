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
           'F,F,F', comment = '#', /fast

  ;; Read the photo-z catalog
  readcol, redshiftCat, $
           zID, ZBEST, $
           ZSPEC_FLAG, ZSPEC_ID, $
           PHOTOZ_FLAG, $
           f = 'L,'+$
           'F,F,'+$
           'F,F,'+$
           'F', comment = '#', /fast

  ;; Check that they have the same number of objects
  if total(ID - zID) ne 0 then begin
     print, ''
     print, ' !!! CATALOGS DO NOT ALIGN -- ABORTING !!!'
     stop
  endif else begin
     print, ''
     print, ' >>> Catalogs look aligned; dumping to output ... '

     savedata = {$
                ID:               0L , $
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

     mwrfits, savedata, outfits, /create
     print, ''
     print, ' >>> Catalog output to : '+outfits
  endelse

end
