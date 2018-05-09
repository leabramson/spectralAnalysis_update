;; Make a "master" GLASS catalog -- WON'T WORK ON MACS1149.
;; For that, use "make1149.pro" contained below.

pro makemaster, FIELD = field, $
                POINTING = pointing, $
                OUTDIR = outdir

  catDir = '$DATA_DIR/GLASS_official_products/'+field+'/Catalogs'

  master    = outdir+'/'+field+'_glass_srclist.fits'
  redshifts = outdir+'/'+field+'_glass_redshifts.fits'
  gigCat    = outdir+'/'+field+'_glass_GiG_results.fits'
  photozOut = outdir+'/'+field+'_Kuang_photozs.fits' 
  
  readglasslist, catDir+'/*glassmaster.cat'    , master
  readredshift , catDir+'/*redshiftcatalog.txt', redshifts
  readgig      , catDir+'/*glassgigcatalog.txt', gigCat

  spawn, 'ls '+catDir+'/*_phot*z*.fits > tmp.list'
  readcol, 'tmp.list', photozs, f = 'A'
  photozs   = photozs[0]
  translatekuang, photozs, photozOut

  print, ''
  print, 'APPENDING GiG RESULTS TO MASTER GLASS SOURCE CATALOG ... '
  print, ''
  join2cats, master, gigCat, 'tmp1.fits', $
             striptag = 'GLASS_ID'
  print, ''
  print, 'APPENDING GiGZ RESULTS TO GLASS SOURCE + GiG CATALOG ... '
  print, ''
  join2cats, 'tmp1.fits', redshifts, 'tmp2.fits', $
             striptag = ['GLASS_ID', 'RA', 'DEC']
  print, ''
  print, 'APPENDING KUANG PHOTO-ZS TO GLASS SOURCE + GiG + GiGz CATALOG ... '
  print, ''
  join2cats, 'tmp2.fits', photozOut, outdir+'/'+field+'_allAvailableData.fits', $
             striptag = ['RA', 'DEC'], $
             dposname = 'KUANG2GLASS_DPOS', $
             POINTING = pointing
  
  print, ''
  print, ' >>> ALL DONE! <<< '
  
end

pro doMost

  makemaster, field = 'ABEL0370', pointing = 'ABEL0370_CTR', outdir = '../MasterCats_V1'
  makemaster, field = 'ABEL2744', pointing = 'ABEL2744_CTR', outdir = '~/Desktop/MasterCats_V1'
;  makemaster, field = 'MACS0416', pointing = 'MACS0416_CTR', outdir =
;  '~/Desktop/MasterCats_V1' ;; MISSING GIG RESULTS
  makemaster, field = 'MACS0717', pointing = 'MACS0717_CTR', outdir = '~/Desktop/MasterCats_V1'
  makemaster, field = 'MACS0744', pointing = 'MACS0744_CTR', outdir = '~/Desktop/MasterCats_V1'
  makemaster, field = 'MACS1423', pointing = 'MACS1423_CTR', outdir = '~/Desktop/MasterCats_V1'
  makemaster, field = 'MACS2129', pointing = 'MACS2129_CTR', outdir = '~/Desktop/MasterCats_V1'
  makemaster, field = 'RXJC1347', pointing = 'RXJC1347_CTR', outdir = '~/Desktop/MasterCats_V1'
  makemaster, field = 'RXJC2248', pointing = 'RXJC2248_CTR', outdir = '../MasterCats_V1'
  
end

pro do1149, outdir = outdir, $
            pointing = pointing

  field = 'MACS1149'
  catDir = '$DATA_DIR/GLASS_official_products/'+field+'/Catalogs'

  master    = outdir+'/'+field+'_glass_srclist.fits'
  redshifts = outdir+'/'+field+'_glass_redshifts.fits'
  gigCat    = outdir+'/'+field+'_glass_GiG_results.fits'
  photozOut = outdir+'/'+field+'_Kuang_photozs.fits' 
  
  readglasslist   , catDir+'/*glassmaster.cat'    , master
  read1149redshift, redshifts  ;;  1149 is a special case
  readgig         , catDir+'/*glassgigcatalog.txt', gigCat

  photozs   = catDir+'/*_photoz.fits'
  translatekuang, photozs, photozOut

  print, ''
  print, 'APPENDING GiG RESULTS TO MASTER GLASS SOURCE CATALOG ... '
  print, ''
  join2cats, master, gigCat, 'tmp1.fits', $
             striptag = 'GLASS_ID'
    print, ''
  print, 'APPENDING REDSHIFTS TO GLASS SOURCE + GiG CATALOG ... '
  print, ''
  join2cats, 'tmp1.fits', redshifts, 'tmp2.fits', $
             striptag = ['RA', 'DEC'], dposname = 'Z2GLASS_DOPS'
  print, ''
  print, 'APPENDING KUANG PHOTO-ZS TO GLASS SOURCE + GiG + GiGz CATALOG ... '
  print, ''
  join2cats, 'tmp2.fits', photozOut, outdir+'/'+field+'_allAvailableData.fits', $
             striptag = ['RA', 'DEC'], $
             dposname = 'KUANG2GLASS_DPOS', $
             POINTING = pointing
  print, ''
  print, ' >>> ALL DONE! <<< '
  
end
;do1149, outdir = '~/Projects/GLASS/MasterCats_V1', pointing = 'MACS1149_CTR'
