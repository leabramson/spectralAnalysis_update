function setzprior, zctr, fwhm, dz, $
                    FLAT = flat

  if N_ELEMENTS(dz) eq 0 then dz = 0.01
  if NOT keyword_set(flat) then flat = 0 else flat = 1
  
  redshift = findgen(10. / dz + 1) * dz

  if NOT flat then begin
     sigma = fwhm / 2.35
     prior = sqrt(2*!pi*sigma^2)^(-1) $
             * exp(-0.5 * (zctr - redshift)^2 / sigma^2)
  endif else begin
     qui = where(abs(redshift - zctr) le FWHM/2., nqui)
     prior = replicate(0.0, 10./dz+1)
     prior[qui] = 1. / FWHM
  endelse
  
  savedata = {Z: redshift, P: prior}
  
  RETURN, savedata
end
