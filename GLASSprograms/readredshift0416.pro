pro readredshift0416

  file = '~/LocalData/Abramson/GLASS_official_products/'+$
         'MACS0416/Catalogs/hlsp_glass_hst_wfc3_macs0416-'+$
         'fullfov-pa999_ir_v001_redshiftcatalog-rmatch0p1.txt'

;ID_zcat
;ID_GLASSv001
;RA
;DEC
;redshift
;redshift_quality
;ID_CLASH_zcat
;CLASH_z
;CLASH_Q
;z_diff

  readcol, file, foo1, glassID, ra, dec, z, zq, $
           f = 'L,L,D,D,F,F'
  ngals = n_elements(glassID)

  savedata = {GLASS_ID         : 0L, $
              RA               : 0.d, $
              DEC              : 0.d, $
              Z_GLASS          : 0., $
              Z_GLASS_QUAL     : 0., $
              Z_GLASS_MULTI_SOL: 0}
  savedata = replicate(savedata, ngals)

  for ii = 0L, ngals - 1 do begin

     savedata[ii].GLASS_ID          = long(glassID[ii])
     savedata[ii].RA                = ra[ii]
     savedata[ii].DEC               = dec[ii]
     savedata[ii].Z_GLASS           = z[ii]
     savedata[ii].Z_GLASS_QUAL      = zq[ii]
     savedata[ii].Z_GLASS_MULTI_SOL = -99
     
  endfor

  mwrfits, savedata, '../MasterCats_V1/MACS0416_glass_redshifts.fits', /create

end
