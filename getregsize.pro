pro getregsize

;  readcol, 'FITS_FILES_ONLY_1.0_1.8_14.0_21.8_c99.9_pacut0_ALLSRCS_QUANTITATIVE_3031.list', files, f= 'A'
  readcol, 'FITS_FILES_ONLY_1.0_1.8_14.0_21.8_c99.9_pacut0_ALLSRCS_QUANTITATIVE.list', files, f= 'A'
  nfiles = n_elements(files)

  names = strarr(nfiles)
  widths = fltarr(nfiles, 2)
  for ii = 0, nfiles - 1 do begin
     d = mrdfits(files[ii],1, /silent)
     names = files[ii]
     widths[ii,*] = [median(d.INTERUP - d.INNERUP), median(d.OUTERUP - d.INTERUP)]
     print, files[ii], widths[ii,0], widths[ii,1]
  endfor
  
;  print, transpose([[names],[reform(widths[*,0])],[reform(widths[*,1])]])
end

pro nada3031
;  ABEL0370/00675_2_B.fits      1.00000      1.00000
;  ABEL0370/01865_2_B.fits      1.00000      1.00000
;  MACS0717/01236_2_B.fits      1.00000      1.00000
  ABEL0370/03762_1_B.fits      5.00000      8.00000 ** bad contam
  ABEL2744/00793_2_B.fits      2.00000      3.00000 *merger/superpos
  MACS0717/01266_1_B.fits      2.00000      3.00000 ** low-z interloper
  MACS0744/00660_2_B.fits      3.00000      4.00000 in
  MACS0744/01475_1_B.fits      3.00000      4.00000 ** low-z interloper
  MACS1149/00900_1_B.fits      3.00000      5.00000 in
  MACS1149/01931_1_B.fits      2.00000      3.00000 ** merger
  MACS1423/01306_2_B.fits      3.00000      4.00000 ** no lgn converge
  MACS2129/00451_2_B.fits      4.00000      6.00000 in
  MACS1423/01916_2_B.fits      2.00000      3.00000 in
  MACS2129/01050_1_B.fits      2.00000      3.00000 ** red/red merger
;  MACS0744/00181_2_B.fits      1.00000      1.00000
;  MACS0744/00660_1_B.fits      1.00000      1.00000
;  MACS0744/01475_2_B.fits      1.00000      1.00000
;  MACS1423/01365_1_B.fits      1.00000      1.00000
;  MACS1423/01365_2_B.fits      1.00000      1.00000
;  MACS2129/00346_1_B.fits      1.00000      1.00000
;  MACS2129/00346_2_B.fits      1.00000      1.00000
;  MACS2129/00431_1_B.fits      1.00000      1.00000
;  MACS2129/00845_2_B.fits      1.00000      1.00000
;  MACS2129/00965_1_B.fits      1.00000      1.00000
;  MACS2129/01126_2_B.fits      1.00000      1.00000
;  MACS2129/01547_1_B.fits      1.00000      1.00000
;  RXJC1347/00752_1_B.fits      1.00000      1.00000
;  RXJC1347/00752_2_B.fits      1.00000      1.00000
;  RXJC2248/00538_2_B.fits      1.00000      1.00000
end

pro nada3030

;  ABEL0370/00675_2_B.fits      1.00000      1.00000
;  ABEL0370/01865_2_B.fits      1.00000      1.00000
;  MACS0717/01236_2_B.fits      1.00000      1.00000
;  MACS0744/00181_2_B.fits      1.00000      1.00000
;  MACS0744/01475_2_B.fits      1.00000      1.00000
  MACS0744/01475_1_B.fits      3.00000      4.00000 *low-z interloper
  ABEL0370/03762_1_B.fits      5.00000      8.00000 *bad contam
  ABEL2744/00793_2_B.fits      2.00000      3.00000 *merger/superpos
  MACS0717/01266_1_B.fits      2.00000      3.00000 *low-z interloper
  MACS1149/00900_1_B.fits      3.00000      5.00000 in
  MACS1149/01931_1_B.fits      2.00000      3.00000 *merger
  MACS1423/01306_2_B.fits      3.00000      4.00000 *no lgn converge
  MACS1423/01916_2_B.fits      2.00000      3.00000 in
  MACS2129/00451_2_B.fits      4.00000      6.00000 in
  MACS2129/01050_1_B.fits      2.00000      3.00000 *red/red merger
;  MACS1423/01365_1_B.fits      1.00000      1.00000
;  MACS1423/01365_2_B.fits      1.00000      1.00000
;  MACS2129/00431_1_B.fits      1.00000      1.00000
;  MACS2129/01547_1_B.fits      1.00000      1.00000
;  RXJC1347/00752_1_B.fits      1.00000      1.00000
;  RXJC2248/00111_2_B.fits      1.00000      1.00000

end
