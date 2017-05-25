function getmag, p, x = X, y = Y, err = Err

  zp = (x - (y + p[0])) / err
  
  RETURN, zp
end

;;
;;

pro spectrophot, spectrum, photometry, $
                 FILTER = filter

  readcol, spectrum, l, f, e, m, $
           comment = '#', /quick, /silent

  phot = mrdfits(photometry, 1, /silent)

  case filter of
     'F814W': file = 'f814w.fil'
     'F125W': file = 'f125w.fil'
     'F140W': file = 'F140W.dat'
     'F160W': file = 'f160w.fil'
  endcase


  filter = fix(strmid(filter, 1,3))
  if filter ge 105 then filter *= 10
  
  readcol, '~/LocalData/Abramson/FILTER_CURVES/'+file, $
           lambda, throughput, f = 'X,F,F'

  resp = interpol(throughput, lambda, l)

  c   = 2.99792458d18  ;; \AA/s
  fnu = f * l^2 / c    ;; 1e-17 erg/s/cm^2/Hz; l * f_l = nu * f_nu

  qui = where(m eq 0 AND e le f)
  tf = total(resp[qui] * fnu[qui] * 1e-17) / total(resp[qui])
  ef = 2.5/alog(10) * total(resp[qui] * e[qui] / f[qui]) / total(resp[qui])
  
  abmag = -2.5 * alog10(tf)

  knownMag = -2.5 * alog10(phot.SED_INNER[where(phot.FILTER eq filter)]) + 23.9
  knownMag = knownMag[0]
  knownErr = 2.5/alog(10) * phot.SED_INNER_ERR[where(phot.FILTER eq filter)] / phot.SED_INNER[where(phot.FILTER eq filter)]
  knownErr = knownErr[0]
  
  p0 = [0.]
  functargs = {X:KnownMag, Y:abMag, Err:sqrt(knownErr^2 + ef^2)}
  norm = mpfit('getmag', p0, functargs = functargs)

  norm = norm[0]  

  print, norm
  
end

pro test
  spectrophot, 'MACS1149/00753_metalFloatResults/00753_1_R_INNFL.masked.pyspec', $
               'MACS1149/00753_1_B_resolvedSED.fits', $
               filter = 'F140W'

end
