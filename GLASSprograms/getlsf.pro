;; - - - - - - - - - - - - - - - - - - - - - - - - - - - -
;; - - - - - - - - - - - - - - - - - - - - - - - - - - - -
;; Get a region-specific LSF; don't stack the whole object
;; - - - - - - - - - - - - - - - - - - - - - - - - - - - -
;; - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function getLSF, di, $
                 X1 = x1, $
                 X2 = x2, $
                 Y1 = y1, $
                 Y2 = y2, $
                 ACS= acs, $
                 DEBUG = debug

  on_error, 2

  if NOT keyword_Set(ACS) then begin
     if N_ELEMENTS(X1) eq 0    then x1 = 20
     if N_ELEMENTS(X2) eq 0    then x2 = 40
     if N_ELEMENTS(Y1) eq 0    then y1 = 20
     if N_ELEMENTS(Y2) eq 0    then y2 = 40
  endif else begin
     if N_ELEMENTS(X1) eq 0    then x1 = 10
     if N_ELEMENTS(X2) eq 0    then x2 = 30
     if N_ELEMENTS(Y1) eq 0    then y1 = 10
     if N_ELEMENTS(Y2) eq 0    then y2 = 30
  endelse
  if NOT keyword_set(DEBUG) then debug = 0 else debug = 1
  
  ;; SPECTRAL DIRECTION
  profile = total(di[x1:x2,y1:y2], 2)
  s       = size(di, /dim)
  lsf     = mpfitpeak(findgen(n_elements(profile)), profile, /moffat, out)
  tx      = findgen(s[0]) - (out[1] + 20)
  u       = tx / out[2]
  lsf     = out[0] / (u^2 + 1)^out[3]
  lsf    /= total(lsf)

  if debug then begin
     im = di
     im[*,y1:y2] = 10 * max(di)
     plot, findgen(n_elements(di[*,0])), findgen(n_elements(di[0,*])), /iso
     cgimage, im, stretch = 2, /over
     k = get_kbrd(1)
  endif
  
  dat = {LSF: LSF, LSF_PARAMS: out}
  
  RETURN, dat; profile
end
