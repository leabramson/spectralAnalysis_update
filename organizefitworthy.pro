pro organizeFitworthy, field, $
                       DOMASKED = domasked

  if NOT keyword_set(DOMASKED) then $
     domasked = 0 $
  else $
     domasked = 1

  ;; Go to the data
  cd, field
  
  spawn, 'ls *.pyspec > TOFIT.list'
  spawn, 'ls *masked.pyspec > TOFIT.MASKED.list'

  ;; Find the things that have extracted spectra (which are "good")
  if NOT domasked then $
     readcol, 'TOFIT.list', tofit, f = 'A' $
  else $
     readcol, 'TOFIT.MASKED.list', tofit, f = 'A' 

  if n_elements(tofit) gt 0 then begin
  
     idsToFit = strarr(n_elements(tofit))
     for ii = 0, n_elements(tofit) - 1 do $
        idsToFit[ii] = strmid(tofit[ii], 0, 7)
     
     idsToFit = idsToFit(sort(idsToFit))
     idsToFit = idsToFit(sort(idsToFit))
     idsToFit = idsToFit(UNIQ(idsToFit))
     
     ntofit = n_elements(idsToFit)
     
     ;; Split ID numbers from PAs
     ids = strarr(ntofit)
     pas = intarr(ntofit)
     for ii = 0, ntofit - 1 do begin
        ids[ii] = strmid(idsToFit[ii], 0, 5)
        pas[ii] = fix(strmid(idsToFit[ii], 6, 1))
     endfor 
     
     ;; Read the master file info list
     readcol, 'fileInfo.list', $
              FILENAME, $
              GLASS_ID, $
              PA, $
              GRISM, $
              Z_PHOT, $
              Z_SPEC, $
              Z_SPEC_QUAL, $
              f = 'A,A,I,A,F,F,F'
     
     ;; Match and print
     close, 1
     openw, 1, 'pyspecWorthy_sourceInfo.list', width = 256
     printf, 1, '#FILENAME GLASS_ID PA GRISM Z_PHOT Z_SPEC Z_SPEC_QUAL'
     for ii = 0, ntofit - 1 do begin
        hit = where(GLASS_ID eq ids[ii] AND PA eq pas[ii], nhit)
        if nhit gt 2 then stop
        for jj = 0, 1 do $
           printf, 1, $
                   FILENAME[hit[jj]], ' ', $
                   GLASS_ID[hit[jj]], ' ', $
                   PA[hit[jj]], ' ', $
                   GRISM[hit[jj]], ' ', $
                   Z_PHOT[hit[jj]], $
                   Z_SPEC[hit[jj]], $
                   Z_SPEC_QUAL[hit[jj]]
     endfor
     close, 1

  endif else begin

     print, ''
     print, ' >>> NO FIT-WORTHY FILES FOUND !!! <<< '
     print, ''
  endelse
  
  cd, '..'
  
end

;;
;;
;;

pro organizeAllFitworthy

  organizefitworthy, 'MACS0717', /domasked
  organizefitworthy, 'MACS0744', /domasked
  organizefitworthy, 'MACS1149', /domasked
  organizefitworthy, 'MACS1423', /domasked
  organizefitworthy, 'MACS2129', /domasked
  organizefitworthy, 'RXJC1347', /domasked
  organizefitworthy, 'RXJC2248', /domasked
  organizefitworthy, 'ABEL2744', /domasked
  
end
