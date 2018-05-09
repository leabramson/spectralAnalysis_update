;; Extract sources in a given redshfit range from the Glass "master"
;; FITS tables

function doContam, infits, $
                   BOTHPA = bothPA, $
                   CLEVEL = clevel
  
  if NOT bothPA then $
     keep = where((    infits.CONTAM_G102_PA1 le clevel    $
                   AND infits.CONTAM_G141_PA1 le clevel    $
                   AND infits.DEFECT_G102_PA1 eq 0         $
                   AND infits.DEFECT_G141_PA1 eq 0)        $
                  OR $
                  (    infits.CONTAM_G102_PA2 le clevel    $
                   AND infits.CONTAM_G141_PA2 le clevel    $
                   AND infits.DEFECT_G102_PA2 eq 0         $
                   AND infits.DEFECT_G141_PA2 eq 0), nuse) $ ;; Take it as long as 1 good PA
  else $
     keep = where((    infits.CONTAM_G102_PA1 le clevel    $
                   AND infits.CONTAM_G141_PA1 le clevel    $
                   AND infits.DEFECT_G102_PA1 eq 0         $
                   AND infits.DEFECT_G141_PA1 eq 0)        $
                  AND $
                  (    infits.CONTAM_G102_PA2 le clevel    $
                   AND infits.CONTAM_G141_PA2 le clevel    $
                   AND infits.DEFECT_G102_PA2 eq 0         $
                   AND infits.DEFECT_G141_PA2 eq 0), nuse) ;; Take it only if 2 good PAs
  
  RETURN, infits[keep]
end

;;
;;
;;

pro cutOnZ, infits, $
            ZMIN     = zmin, $
            ZMAX     = zmax, $
            PHOTOZ   = photoz, $
            DOCONTAM = docontam, $
            BOTHPA   = bothPa, $
            CLEVEL   = clevel, $
            OUTPUT   = output
  
  if NOT keyword_set(ZMIN) then zmin = 0.3
  if NOT keyword_set(ZMIN) then zmax = 0.7
  if NOT keyword_set(PHOTOZ) then photoz = 0 else photoz = 1
  if NOT keyword_set(DOCONTAM) then docontam = 0 else docontam = 1
  if docontam then begin
     if NOT keyword_set(CLEVEL) then clevel = 1.5
     if NOT keyword_set(BOTHPA) then bothPA = 0 else bothPA = 1
  endif
  
;  readcol, fileList, files, f = 'A'
;  nfiles = n_elements(files)

  fields = []
  ids    = []
  ras    = []
  decs   = []
  zs     = []
  zqs    = []
  pzs    = []
  mags   = []
  
  d = mrdfits(infits, 1, /silent)
  use = where(d.Z_GLASS ge zmin AND d.Z_GLASS lt zmax, nuse)

  if photoz then $
     pzadd = where(d.KUANG_PZ_A ge zmin AND d.KUANG_PZ_A lt zmax, nadd)

  use = [use,pzadd]
  tu  = use[sort(use)]
  use = tu[UNIQ(tu)]

  d2 = d[use]

  if docontam then $
     d2 = docontam(d2, BOTHPA = bothPA, CLEVEL = clevel)
  
  fields = [fields, d2.POINTING]
  ids    = [ids,    d2.NUMBER]
  ras    = [ras,    d2.X_WORLD]
  decs   = [decs,   d2.Y_WORLD]
  zs     = [zs,     d2.Z_GLASS]
  zqs    = [zqs,    d2.Z_GLASS_QUAL]
  pzs    = [pzs,    d2.KUANG_PZ_A]
  mags   = [mags,   d2.MAG_AUTO]
  
  nToPrint = n_elements(ids)
  
  close, 1
  openw, 1, output, width = 128
  printf, 1, '#FIELD GLASS_ID RA DEC Z ZQ KZ_PHOT MAG'
  for ii = 0, nToPrint - 1 do $
     printf, 1, fields[ii], $
             ids[ii], $
             ras[ii], $
             decs[ii], $
             zs[ii], $
             zqs[ii], $
             pzs[ii], $
             mags[ii]
  close, 1

  print, ''
  print, string(nToPrint, f = '(I)')+' galaxies at '+$
         string(zmin, f = '(F4.2)')+' < z < '+$
         string(zmax, f = '(F4.2)')+' dumped to '+output
  
end
