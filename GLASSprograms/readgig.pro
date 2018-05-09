pro readGiG, gigIn, fitsOut

;; ID RA DEC
;; CLASH_ID CLASH_DISTANCE MAG_SELECTION MAGERR_SELECTION
;; CONTAM_G102_20 DEFECT_G102_20 CONTAM_DEFECT_G102_20
;; CONTAM_G102_280 CONTAM_DEFECT_G102_280 DEFECT_G102_280
;; CONTAM_G141_20 CONTAM_DEFECT_G141_20 DEFECT_G141_20
;; CONTAM_G141_280    CONTAM_DEFECT_G141_280 DEFECT_G141_280
;; STAR

  readcol, gigIn, $
           glassID, ra, dec, $
           clashID, clash_dpos, magsel, emagsel, $
           contam_102_p1, defect_102_p1, contam_defect_102_p1, $
           contam_102_p2, defect_102_p2, contam_defect_102_p2, $
           contam_141_p1, defect_141_p1, contam_defect_141_p1, $
           contam_141_p2, defect_141_p2, contam_defect_141_p2, $
           starflag, $
           f = 'A,D,D,'+$
           'A,D,F,F,'+$
           'F,F,F,'+$
           'F,F,F,'+$
           'F,F,F,'+$
           'F,F,F,'+$
           'F', comment = '#'
  ngals = n_elements(glassID)
  
  spawn, 'cat '+gigIn+' | grep "PA1=" > tmp.txt'
  readcol, 'tmp.txt', str, f = 'A', $
           comment = '', delimit = '_'

  x1 = strpos(str, "PA1=")
  x2 = strpos(str, "PA2=")
  pa1 = strmid(str, x1+4, 3)
  pa2 = strmid(str, x2+4, 3)

  savedata = {GLASS_ID: 0L, RA: 0.d, DEC: 0.d, PA1: 0., PA2: 0., $
              CLASH_ID: 0L, CLASH_DPOS: 0.d, $
              MAG_SELECT: 0., E_MAG_SELECT:0., $
              CONTAM_G102_PA1: 0., CONTAM_G102_PA2: 0., $
              DEFECT_G102_PA1: 0., DEFECT_G102_PA2: 0., $
              CONDEF_G102_PA1: 0., CONDEF_G102_PA2: 0., $
              CONTAM_G141_PA1: 0., CONTAM_G141_PA2: 0., $
              DEFECT_G141_PA1: 0., DEFECT_G141_PA2: 0., $
              CONDEF_G141_PA1: 0., CONDEF_G141_PA2: 0.}
  savedata = replicate(savedata, ngals)

  for ii = 0, ngals - 1 do begin

     savedata[ii].GLASS_ID        = long(glassId[ii])
     savedata[ii].RA              = ra[ii]
     savedata[ii].DEC             = dec[ii]
     savedata[ii].PA1             = float(pa1)
     savedata[ii].PA2             = float(pa2)
     savedata[ii].CLASH_ID        = clashID[ii]
     savedata[ii].CLASH_DPOS      = clash_dpos[ii]
     savedata[ii].MAG_SELECT      = magsel[ii]
     savedata[ii].E_MAG_SELECT    = emagsel[ii]
     savedata[ii].CONTAM_G102_PA1 = contam_102_p1[ii]
     savedata[ii].CONTAM_G102_PA2 = contam_102_p2[ii]
     savedata[ii].DEFECT_G102_PA1 = defect_102_p1[ii]
     savedata[ii].DEFECT_G102_PA2 = defect_102_p2[ii]
     savedata[ii].CONDEF_G102_PA1 = contam_defect_102_p1[ii]
     savedata[ii].CONDEF_G102_PA2 = contam_defect_102_p2[ii]
     savedata[ii].CONTAM_G141_PA1 = contam_141_p1[ii]
     savedata[ii].CONTAM_G141_PA2 = contam_141_p2[ii]
     savedata[ii].DEFECT_G141_PA1 = defect_141_p1[ii]
     savedata[ii].DEFECT_G141_PA2 = defect_141_p2[ii]
     savedata[ii].CONDEF_G141_PA1 = contam_defect_141_p1[ii]
     savedata[ii].CONDEF_G141_PA2 = contam_defect_141_p2[ii]
     
  endfor

  mwrfits, savedata, fitsOut, /create
  
end

pro do0717

  readGiG, 'GLASS_OFFICIAL/GIG_RESULTS/hlsp_glass_hst_wfc3_macs0717-fullfov-pa999_ir_v001_glassgigcatalog.txt', $
           'M0717_gig_results.fits'
end
