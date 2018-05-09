function loadsummarydata, dir

  if n_elements(dir) eq 0 then $
     dir = '..'

;  #FIELD ID+PA RA DEC ZFIT MU EMU SURVEY Z_GLASS ZQ KZ_PHOT MAG_SEL  
  
  readcol, dir+'/'+'20170613_finalSampleDetails.dat', $
           field, id, ra, dec, z, mag, emag, survey, $
           z_glass, zq_glass, z_phot, f140_select, $
           f = 'A,A,D,D,F,F,F,A,F,F,F,F', /silent, comment = '#'
  ngals = n_elements(id)
  
  summary = {FIELD:   string(0), $
             ID:      string(0), $
             RA:      0.d, $
             DEC:     0.d, $
             Z:       0., $
             MU:      0., $
             EMU:     0., $
             SURVEY:  string(0), $
             Z_GLASS: 0., $
             Z_QUAL:  0., $
             KZ_PHOT: 0., $
             MAG_SEL: 0.}
  summary = replicate(summary, ngals)

  for ii = 0, ngals - 1 do begin

     summary[ii].FIELD   = field[ii]
     summary[ii].ID      = id[ii]
     summary[ii].RA      = ra[ii]
     summary[ii].DEC     = dec[ii]
     summary[ii].Z       = z[ii]
     summary[ii].MU      = mag[ii]
     summary[ii].EMU     = emag[ii]
     summary[ii].SURVEY  = survey[ii]
     summary[ii].Z_GLASS = z_glass[ii]
     summary[ii].Z_QUAL  = zq_glass[ii]
     summary[ii].KZ_PHOT = z_phot[ii]
     summary[ii].MAG_SEL = f140_select[ii]
     
  endfor

  RETURN, summary
  
end
