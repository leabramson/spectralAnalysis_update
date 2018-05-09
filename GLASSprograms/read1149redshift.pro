;; reads the redshift list (many sources) for MACS1149 supplied by
;; the GLASS MAST HST archive cite. File is the "master" redshift file
;; "hlsp_glass_hst-vlt-keck_wfc3-muse-deimos_macs1149_optical-ir_v001_redshiftcatalog.txt" not the one pre-crossmatched by Kasper.

pro read1149redshift, outfits

  ;;ID   RA    DEC
  ;;redshift  redshift_quality redshift_src (1=HST,2=MUSE,
  ;;3=HST+MUSE,4=MUSE+DEIMOS)
  
  readcol, '$GLASS_DIR/MACS1149/Catalogs/'+$
           'hlsp_glass_hst-vlt-keck_wfc3-'+$
           'muse-deimos_macs1149_optical-ir'+$
           '_v001_redshiftcatalog.txt', $
           zID, ra, dec, $
           z, zq, src, $
           f = 'L,D,D,'+$
           'F,F,I'
  ngals = n_elements(zID)
  
  output = {SPEC_ID:      0L , $
            RA:           0.d, $
            DEC:          0.d, $
            Z_GLASS:      0. , $
            Z_GLASS_QUAL: 0. , $
            Z_SRC:        0}
  output = replicate(output, ngals)

  for ii = 0, ngals - 1 do begin
     output[ii].SPEC_ID      = zID[ii]
     output[ii].RA           = ra[ii]
     output[ii].DEC          = dec[ii]
     output[ii].Z_GLASS      = z[ii]
     output[ii].Z_GLASS_QUAL = zq[ii]
     output[ii].Z_SRC        = src[ii]
  endfor

  mwrfits, output, outfits, /create
  
end
