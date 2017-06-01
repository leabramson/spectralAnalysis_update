;; find the average n2/ha ratio given an oiii/hb ratio using the fit
;; from Kewley+13.

function deblendhalpha, ha, $
                        oiii = oiii, $
                        hb = hb, $
                        z = z

  logO = alog10(oiii/hb)

  logN = 0.61 / (logO - 1.2 - 0.03 * z) $
         +0.02 + 0.1833 * z

  n2ha = 10.^logN 

  ha_corr = (1 - n2ha) * ha
  
  RETURN, ha_corr
end
