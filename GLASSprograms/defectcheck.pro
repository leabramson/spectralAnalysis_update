function defectcheck, inlist, $
                      THRESH = thresh

  if NOT keyword_set(THRESH) then thresh = 0.25 ;; Default to 75% good pix

  readcol, inlist, files, f = 'A', /silent
  nfiles = n_elements(files)

  keep = []
  for ii = 0, nfiles - 1 do begin

     data = mrdfits(files[ii], 'sci', head, /silent)
     check = total(data lt 0) / n_elements(data)  ;;  GLASS spectra have negative null pix

     if check le thresh then $
        keep = [keep, ii]
     
  endfor

  RETURN, keep
  
end
