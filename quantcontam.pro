;; Return the fraction of pixels w/in a given extraction region
;; that is contaminated to a certain percentage of the source flux

function quantContam, infits, $
                      CFRACT = cfract, $
                      DETECTONLY = detectonly

  data  = mrdfits(infits, 1, /silent)
  file  = data.FILENAME

;  sci    = mrdfits(file, 'sci', /silent)
  sci    = mrdfits(file, 'model', /silent)
  err    = mrdfits(file, 'wht', /silent)
  contam = mrdfits(file, 'contam', /silent)
 
  cfrac = contam / sci ;; Contamination fraction map
  zones = data.EXTRACT_ZONES
  
  if keyword_set(DETECTONLY) then begin
     use = where(sci/err gt 2, compl = bad)
     cfrac = cfrac[use]
     zones = zones[use]
  endif
  
  clevels = {TOTAL   : total(cfrac ge cfract) / n_elements(cfrac), $
             OPT_ZONE: total(cfrac[where(zones eq  5)] ge cfract) / total(zones eq  5), $
             OUT_ZONE: total(cfrac[where(zones eq 10)] ge cfract) / total(zones eq 10), $
             INT_ZONE: total(cfrac[where(zones eq 20)] ge cfract) / total(zones eq 20), $
             INN_ZONE: total(cfrac[where(zones eq 40)] ge cfract) / total(zones eq 40), $ 
             INN_NE  : total(data.SPECORIG[where(zones eq 40)] le 0) / total(zones eq 40)} ; total(data.F_INNER le 0) / n_elements(data.F_INNER)
  
  RETURN, clevels
end
