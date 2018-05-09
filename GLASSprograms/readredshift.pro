pro readRedshift, catIn, fitsOut

  ;;ID   RA    DEC
  ;;redshift  redshift_quality   multiple_redshift_solutions   

  readcol, catIn, $
           glassID, ra, dec, $
           z, zq, mzs, $
           f = 'A,D,D,F,F,F'
  ngals = n_elements(glassID)

  savedata = {GLASS_ID         : 0L, $
              RA               : 0.d, $
              DEC              : 0.d, $
              Z_GLASS          : 0., $
              Z_GLASS_QUAL     : 0., $
              Z_GLASS_MULTI_SOL: 0.}
  savedata = replicate(savedata, ngals)

  for ii = 0, ngals - 1 do begin

     savedata[ii].GLASS_ID          = long(glassID[ii])
     savedata[ii].RA                = ra[ii]
     savedata[ii].DEC               = dec[ii]
     savedata[ii].Z_GLASS           = z[ii]
     savedata[ii].Z_GLASS_QUAL      = zq[ii]
     savedata[ii].Z_GLASS_MULTI_SOL = mzs[ii]
     
  endfor

  mwrfits, savedata, fitsOut, /create
  
end
