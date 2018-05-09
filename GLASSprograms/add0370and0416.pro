;; For some reason, MACS0416 and ABEL0370 lack all of the standard
;; data comapred to the other fields. This program will pull out the
;; necessary "cutOnZ" and "cutOnMag"-formatted files for whatever you
;; can scrounge together for those fields

pro add0370and0416

  add0370
  add0416
  
end

pro add0370

  ;; Do 0370
  master    = '../MasterCats_V1/ABEL0370_glass_srclist.fits'
  gigCat    = '../MasterCats_V1/ABEL0370_glass_GiG_results.fits'
  redshifts = '../MasterCats_V1/ABEL0370_glass_redshifts.fits'
  tids     = '/data/code/mtaka/HFF/Public/SEXT/f160w_05.cat'
  photo_z  = '/data/code/mtaka/HFF/Public/EAZY/photz_05.zout'
  
  print, f = $
         '(%"\n\nAPPENDING ABEL0370 GiG RESULTS TO MASTER GLASS SOURCE CATALOG ...\n\n")'
  join2cats, master, gigCat, 'tmp1.fits', $
             striptag = 'GLASS_ID'

  print, f = $
         '(%"\n\nAPPENDING ABEL0370 GiGZ RESULTS TO GLASS SOURCE + GiG CATALOG ...\n\n")'
  join2cats, 'tmp1.fits', redshifts, 'tmp2.fits', $
             striptag = ['GLASS_ID', 'RA', 'DEC']

  print, f = $
         '(%"\n\nTRANSLATING ABEL0370 TAKAHIRO CATALOG ...\n\n")'
  readtakamaster, tids, 'ABEL0370_takaID.fits'

  print, f = $
         '(%"\n\nTRANSLATING ABEL0370 TAKAHIRO PHOTO-Z CATALOG ...\n\n")'
  readtakaphotoz, photo_z, 'ABEL0370_takaPhotoz.fits'

  print, f = $
         '(%"\n\nCOMBINING ABEL0370 TAKAHIRO CATALOGS ...\n\n")'
  
  tid = mrdfits('ABEL0370_takaID.fits', 1)
  tpz = mrdfits('ABEL0370_takaPhotoz.fits', 1)

  ttot = joinOnIDs(tid, tpz, striptag = 'ID_TAKA')
  mwrfits, ttot, 'ABEL0370_taka_IDplusPZ.fits', /create
  
  print, f = $
         '(%"\n\nCOMBINING ABEL0370 TAKAHIRO CATALOGS W/ GLASS ...\n\n")'
  join2cats, 'tmp2.fits', 'ABEL0370_taka_IDplusPZ.fits', '../MasterCats_V1/ABEL0370_allAvailableDataTaka.fits', $
             striptag = ['MAG_AUTO', 'MAGERR_AUTO', 'KRON_RADIUS', 'BACKGROUND', 'RA', 'DEC',$
                         'THETA_IMAGE', 'A_IMAGE', 'B_IMAGE', 'CLASS_STAR', 'FLUX_RADIUS', 'FLAGS', $
                         'FLUX_AUTO', 'FLUXERR_AUTO', 'FLUX_APER2', 'FLUXERR_APER2', 'FLUX_APER3', 'FLUXERR_APER3'], $
             POINTING = 'ABEL0370_CTR'
  
end

pro add0416

  master    = '../MasterCats_V1/MACS0416_glass_srclist.fits'
  redshifts = '../MasterCats_V1/MACS0416_glass_redshifts.fits'
  tids     = '/data/code/mtaka/HFF/Public/SEXT/f160w_02.cat'
  photo_z  = '/data/code/mtaka/HFF/Public/EAZY/photz_02.zout'
  
  print, f = $
         '(%"\n\nAPPENDING MACS0416 GiGZ RESULTS TO GLASS SOURCE CATALOG ...\n\n")'
  join2cats, master, redshifts, 'tmp2.fits', $
             striptag = 'GLASS_ID'

  print, f = $
         '(%"\n\nTRANSLATING MACS0416 TAKAHIRO CATALOG ...\n\n")'
  readtakamaster, tids, 'MACS0416_takaID.fits'

  print, f = $
         '(%"\n\nTRANSLATING MACS0416 TAKAHIRO PHOTO-Z CATALOG ...\n\n")'
  readtakaphotoz, photo_z, 'MACS0416_takaPhotoz.fits'

  print, f = $
         '(%"\n\nCOMBINING MACS0416 TAKAHIRO CATALOGS ...\n\n")'
  
  tid = mrdfits('MACS0416_takaID.fits', 1)
  tpz = mrdfits('MACS0416_takaPhotoz.fits', 1)

  ttot = joinOnIDs(tid, tpz, striptag = 'ID_TAKA')
  mwrfits, ttot, 'MACS0416_taka_IDplusPZ.fits', /create
  
  print, f = $
         '(%"\n\nCOMBINING MACS0416 TAKAHIRO CATALOGS W/ GLASS ...\n\n")'
  join2cats, 'tmp2.fits', 'MACS0416_taka_IDplusPZ.fits', '../MasterCats_V1/MACS0416_allAvailableDataTaka.fits', $
             striptag = ['MAG_AUTO', 'MAGERR_AUTO', 'KRON_RADIUS', 'BACKGROUND', 'RA', 'DEC',$
                         'THETA_IMAGE', 'A_IMAGE', 'B_IMAGE', 'CLASS_STAR', 'FLUX_RADIUS', 'FLAGS', $
                         'FLUX_AUTO', 'FLUXERR_AUTO', 'FLUX_APER2', 'FLUXERR_APER2', 'FLUX_APER3', 'FLUXERR_APER3'], $
             POINTING = 'MACS0416_CTR'
end
