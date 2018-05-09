;; - - - - - - - - - - - - - - - - - - - - - - - - - - - -
;; - - - - - - - - - - - - - - - - - - - - - - - - - - - -
;; Get a region-specific LSF; don't stack the whole object
;; - - - - - - - - - - - - - - - - - - - - - - - - - - - -
;; - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function getLSF, di, $
                 X1 = x1, $
                 X2 = x2, $
                 Y1 = y1, $
                 Y2 = y2

  if NOT keyword_set(X1) then x1 = 20
  if NOT keyword_set(X2) then x2 = 40
  if NOT keyword_set(Y1) then y1 = 20
  if NOT keyword_set(Y2) then y2 = 40
  
  ;; SPECTRAL DIRECTION
  profile = total(di[x1:x2,y1:y2], 2)
  s       = size(di, /dim)
  lsf     = mpfitpeak(findgen(n_elements(profile)), profile, /moffat, out)
  tx      = findgen(s[0]) - (out[1] + 20)
  u       = tx / out[2]
  lsf     = out[0] / (u^2 + 1)^out[3]
  lsf /= total(lsf)

  dat = {LSF: LSF}
  
  RETURN, dat; profile
end
