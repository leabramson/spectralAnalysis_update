function getmag, p, x = X, y = Y, err = Err

  zp = (x - (y + p[0])) / err
  
  RETURN, zp
end

;;
;;

function spectrophot, spectrum, photometry, $
                 ZONE = zone, $
                 FILTER = filter

  if NOT keyword_set(ZONE) then zone = 'INNFL'
  
;  readcol, spectrum, l, f, e, m, $
;           comment = '#', /quick, /silent

  data = mrdfits(spectrum, 1)
  l = data.lambda
  case zone of
     'INNFL': begin
          f = data.F_INNER        
          e = sqrt(data.VAR_INNER)
          m = data.MASK_INNER
       end
     'INTUP': begin
          f = data.F_INTER_UP       
          e = sqrt(data.VAR_INTER_UP)
          m = data.MASK_INTER_UP
     end
     'INTDN': begin
          f = data.F_INTER_DN     
          e = sqrt(data.VAR_INTER_DN)
          m = data.MASK_INTER_DN
     end
     'OUTUP': begin
          f = data.F_OUTER_UP    
          e = sqrt(data.VAR_OUTER_UP)
          m = data.MASK_OUTER_UP
     end
     'OUTDN': begin
          f = data.F_OUTER_DN    
          e = sqrt(data.VAR_OUTER_DN)
          m = data.MASK_OUTER_DN
     end
  endcase
  f*= 1d-17 / data.SENSITIVITY / data.EXPTIME
  e*= 1d-17 / data.SENSITIVITY / data.EXPTIME
     
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
  fnu = f * l^2 / c    ;; erg/s/cm^2/Hz; l * f_l = nu * f_nu
  
  qui = where(m eq 0 AND e le f)
  tf = total(resp[qui] * fnu[qui]) / total(resp[qui])
  ef = 2.5/alog(10) * total(resp[qui] * e[qui] / f[qui]) / total(resp[qui])
  
  abmag = -2.5 * alog10(tf) - 48.6

  case zone of
     'INNFL': begin
        pd = phot[where(phot.FILTER eq filter)].SED_INNER
        pe = phot[where(phot.FILTER eq filter)].SED_INNER_ERR / $
             phot[where(phot.FILTER eq filter)].SED_INNER
     end
     'INTUP': begin
        pd = phot[where(phot.FILTER eq filter)].SED_INTER_UP
        pe = phot[where(phot.FILTER eq filter)].SED_INTER_UP_ERR / $
             phot[where(phot.FILTER eq filter)].SED_INTER_UP
     end
     'INTDN': begin
        pd = phot[where(phot.FILTER eq filter)].SED_INTER_DN
        pe = phot[where(phot.FILTER eq filter)].SED_INTER_DN_ERR / $
             phot[where(phot.FILTER eq filter)].SED_INTER_DN
     end
     'OUTUP': begin
        pd = phot[where(phot.FILTER eq filter)].SED_OUTER_UP
        pe = phot[where(phot.FILTER eq filter)].SED_OUTER_UP_ERR / $
             phot[where(phot.FILTER eq filter)].SED_OUTER_UP
     end
     'OUTDN': begin
        pd = phot[where(phot.FILTER eq filter)].SED_OUTER_UP
        pe = phot[where(phot.FILTER eq filter)].SED_OUTER_UP_ERR / $
             phot[where(phot.FILTER eq filter)].SED_OUTER_UP
     end
  endcase
  
  knownMag = -2.5 * alog10(pd) + 23.9
  knownMag = knownMag[0]
  knownErr = 2.5/alog(10) * pe
  knownErr = knownErr[0]

  print, knownMag, knownErr
  print, abmag, ef
  print, 10.^(-(abmag-knownMag)/2.5), 10.^(sqrt(knownErr^2 + ef^2)/2.5) - 1.
  
;  stop
;  
;  p0 = [0.]
;  functargs = {X:KnownMag, Y:abMag, Err:sqrt(knownErr^2 + ef^2)}
;  norm = mpfit('getmag', p0, functargs = functargs)
;
;  norm = norm[0]  
;
;  print, norm

  RETURN, [10.^(-(abmag-knownMag)/2.5), 10.^(sqrt(knownErr^2 + ef^2)/2.5) - 1.]
  
end

pro test
;  spectrophot, 'MACS1149/00753_metalFloatResults/00753_1_R_INNFL.masked.pyspec', $
;               'MACS1149/00753_1_B_resolvedSED.fits', $
;               filter = 'F140W'

  m9inn = spectrophot('MACS1149/00900_1_R.fits', $
                      'MACS1149/00900_1_B_resolvedSED.fits', $
                      zone = 'INNFL', $
                      filter = 'F140W')
  m9intup = spectrophot('MACS1149/00900_1_R.fits', $
                        'MACS1149/00900_1_B_resolvedSED.fits', $
                        zone = 'INTUP', $
                        filter = 'F140W')
  m9intdn = spectrophot('MACS1149/00900_1_R.fits', $
                        'MACS1149/00900_1_B_resolvedSED.fits', $
                        zone = 'INTDN', $
                        filter = 'F140W')
  m9outup = spectrophot('MACS1149/00900_1_R.fits', $
                        'MACS1149/00900_1_B_resolvedSED.fits', $
                        zone = 'OUTUP', $
                        filter = 'F140W')
  m9outdn = spectrophot('MACS1149/00900_1_R.fits', $
                        'MACS1149/00900_1_B_resolvedSED.fits', $
                        zone = 'OUTDN', $
                        filter = 'F140W')
  m9 = [m9inn[0], m9intup[0], m9intdn[0], m9outup[0], m9outdn[0]]
  m9e = [m9inn[1], m9intup[1], m9intdn[1], m9outup[1], m9outdn[1]]
  
  m6inn = spectrophot('MACS0744/00660_2_R.fits', $
                      'MACS0744/00660_2_B_resolvedSED.fits', $
                      zone = 'INNFL', $
                      filter = 'F140W')
  m6intup = spectrophot('MACS0744/00660_2_R.fits', $
                        'MACS0744/00660_2_B_resolvedSED.fits', $
                        zone = 'INTUP', $
                        filter = 'F140W')
  m6intdn = spectrophot('MACS0744/00660_2_R.fits', $
                        'MACS0744/00660_2_B_resolvedSED.fits', $
                        zone = 'INTDN', $
                        filter = 'F140W')
  m6outup = spectrophot('MACS0744/00660_2_R.fits', $
                        'MACS0744/00660_2_B_resolvedSED.fits', $
                        zone = 'OUTUP', $
                        filter = 'F140W')
  m6outdn = spectrophot('MACS0744/00660_2_R.fits', $
                        'MACS0744/00660_2_B_resolvedSED.fits', $
                        zone = 'OUTDN', $
                        filter = 'F140W')
  m6 = [m6inn[0], m6intup[0], m6intdn[0], m6outup[0], m6outdn[0]]
  m6e = [m6inn[1], m6intup[1], m6intdn[1], m6outup[1], m6outdn[1]]

  m19inn = spectrophot('MACS1423/01916_2_R.fits', $
                      'MACS1423/01916_2_B_resolvedSED.fits', $
                      zone = 'INNFL', $
                      filter = 'F140W')
  m19intup = spectrophot('MACS1423/01916_2_R.fits', $
                        'MACS1423/01916_2_B_resolvedSED.fits', $
                        zone = 'INTUP', $
                        filter = 'F140W')
  m19intdn = spectrophot('MACS1423/01916_2_R.fits', $
                        'MACS1423/01916_2_B_resolvedSED.fits', $
                        zone = 'INTDN', $
                        filter = 'F140W')
  m19outup = spectrophot('MACS1423/01916_2_R.fits', $
                        'MACS1423/01916_2_B_resolvedSED.fits', $
                        zone = 'OUTUP', $
                        filter = 'F140W')
  m19outdn = spectrophot('MACS1423/01916_2_R.fits', $
                        'MACS1423/01916_2_B_resolvedSED.fits', $
                        zone = 'OUTDN', $
                        filter = 'F140W')
  m19 = [m19inn[0], m19intup[0], m19intdn[0], m19outup[0], m19outdn[0]]
  m19e = [m19inn[1], m19intup[1], m19intdn[1], m19outup[1], m19outdn[1]]

  m4inn = spectrophot('MACS2129/00451_2_R.fits', $
                      'MACS2129/00451_2_B_resolvedSED.fits', $
                      zone = 'INNFL', $
                      filter = 'F140W')
  m4intup = spectrophot('MACS2129/00451_2_R.fits', $
                        'MACS2129/00451_2_B_resolvedSED.fits', $
                        zone = 'INTUP', $
                        filter = 'F140W')
  m4intdn = spectrophot('MACS2129/00451_2_R.fits', $
                        'MACS2129/00451_2_B_resolvedSED.fits', $
                        zone = 'INTDN', $
                        filter = 'F140W')
  m4outup = spectrophot('MACS2129/00451_2_R.fits', $
                        'MACS2129/00451_2_B_resolvedSED.fits', $
                        zone = 'OUTUP', $
                        filter = 'F140W')
  m4outdn = spectrophot('MACS2129/00451_2_R.fits', $
                        'MACS2129/00451_2_B_resolvedSED.fits', $
                        zone = 'OUTDN', $
                        filter = 'F140W')
  m4 = [m4inn[0], m4intup[0], m4intdn[0], m4outup[0], m4outdn[0]]
  m4e = [m4inn[1], m4intup[1], m4intdn[1], m4outup[1], m4outdn[1]]
  
  print, ''
  print, mean(m9), mean(m6), mean(m19), mean(m4)
  print, mean(abs(1-m9)/m9e), mean(abs(1-m6)/m6e), mean(abs(1-m19)/m19e), mean(abs(1-m4)/m4e)

  h = histogram([m4,m6,m9,m19], min = 0, max = 10, bins= 0.1, loc = bins)

  stop

end
