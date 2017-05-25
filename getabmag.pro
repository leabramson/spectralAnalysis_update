;; Compute an F140W AB mag from a unified pyspec input spectrum

function getabmag, spectrum, master, $
                   FILTER = filter

  readcol, spectrum, l, f, e, m

  mstr = mrdfits(master, 1, /silent)
  exptime = mstr.EXPTIME
;  mastermag = mstr.MAG
  
  f /= exptime  ;; 1e-17 erg/s/cm^2/A

  case filter of
     'F814W': file = 'f814w.fil'
     'F125W': file = 'f125w.fil'
     'F140W': file = 'F140W.dat'
     'F160W': file = 'f160w.fil'
  endcase

  readcol, '~/LocalData/Abramson/FILTER_CURVES/'+file, $
           lambda, throughput, f = 'X,F,F'

  resp = interpol(throughput, lambda, l)

  c = 2.99792458d18  ;; \AA/s
  fnu  = l^2 / c * f ;; 1e-17 erg/s/cm^2/Hz; l * f_l = nu * f_nu

  qui = where(m eq 0)
  
;  abmag     = -2.5 * alog10(tsum(l[qui], resp[qui] * fnu[qui] *
;  1e-17) / tsum(l[qui], resp[qui])) - 48.6
  abmag     = -2.5 * alog10(total(resp[qui] * fnu[qui] * 1e-17) / total(resp[qui])) - 48.6
;  offTomstr = abmag - mastermag
  
  RETURN, abmag
end
