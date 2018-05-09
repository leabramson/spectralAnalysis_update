function getmagnification, inlist

  readcol, inlist, model, mag, $
           f = 'A,F', comment = '#', /silent

  m = sort(mag)

  mag   = mag[m]
  model = model[m]

  ;; Kill top and bottom estimate

  killMag   = [mag[0], mag[-1]]
  killModel = [model[0], model[-1]]
  
  mag   = mag[1:-2]
  model = model[1:-2]

  results = {MEDIAN_MAG: median(mag), $
             MEAN_MAG: mean(mag, /nan), $
             MAG_ERR: stddev(mag, /nan), $
             MINMAG: min(mag), $
             MAXMAG: max(mag), $
             MAGS: mag, $
             MODELS: model, $
             EXCLUDED_MAGS: killMag, $
             EXCLUDED_MODELS: killModel}  

  RETURN, results
  
end
