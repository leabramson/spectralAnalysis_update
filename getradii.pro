function getradii, infitsn

  ;; Get spatial & wavelegth scales
  infits = mrdfits(infitsn, 1)
  bdi = infits.DIRECT_IM
  s = size(bDI, /dim)
  foo = mrdfits(infits.FILENAME, 'sci', head, /silent)
  pscale = sxpar(head, 'CD2_2') ;; It's downsampled by ~2x from WFC3 plate scale
  lscale = sxpar(head, 'CD1_1')
  pivot  = sxpar(head, 'CRVAL1')
;  tpixz  = lscale / pivot * 10 * (1+z) ;; 10 pix offset in z space

  run = (findgen(s[0]) - (s[0]+1)/2) * pscale
  run += pscale/2

;  kpc       = 1. / zang(1.0, z, /silent)    ;; kpc / asec
;  kpcPix    = kpc * pscale                  ;; kpc / pix
;  hlrKpc    = infits.RE * kpcPix            ;; re in Kpc
  hlrObs    = infits.RE * pscale            ;; re in asec


  regions = ['INNFL', 'INTER', 'OUTER', $
             'INTUP', 'INTDN', 'OUTUP', 'OUTDN']
  
  radii = {FILE: string(0), $
           REGION: string(0), $
           RNORM: 0. ,$
           RE_OBS: 0.}
  radii = replicate(radii, n_elements(regions))
  radii[*].FILE   = infitsn
  radii[*].RE_OBS = hlrObs
  for ii = 0, n_elements(regions) - 1 do begin
     radii[ii].REGION = regions[ii]
     case regions[ii] of
        'INNFL': radii[ii].RNORM = 0
        'INTER': radii[ii].RNORM =  0.5 * (infits.EXTR_INTUP_OFF + infits.EXTR_INTDN_OFF) / infits.RE
        'OUTER': radii[ii].RNORM =  0.5 * (infits.EXTR_OUTUP_OFF + infits.EXTR_OUTDN_OFF) / infits.RE
        'INTUP': radii[ii].RNORM =  infits.EXTR_INTUP_OFF / infits.RE
        'INTDN': radii[ii].RNORM = -infits.EXTR_INTDN_OFF / infits.RE
        'OUTUP': radii[ii].RNORM =  infits.EXTR_OUTUP_OFF / infits.RE
        'OUTDN': radii[ii].RNORM = -infits.EXTR_OUTDN_OFF / infits.RE
     endcase
  endfor

  RETURN, radii
  
end
