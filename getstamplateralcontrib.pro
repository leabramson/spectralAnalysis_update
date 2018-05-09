;; Actually extract photometry from the stamps made by getsrcimage.pro

function getstamplateralcontrib, master = master, $
                                 multiband = multiband, $
                                 zeropoint = zeropoint, $
                                 cluster = cluster, $
                                 dustfile = dustfile, $
                                 OUTPUT = output, $
                                 SHIFT = shift

  if NOT keyword_set(zeropoint) then zeropoint = 23.9 ;; AB mags
  if NOT keyword_set(dustfile) then dustfile = '../'+cluster+'/GALACTIC_EXTINCTION.dat'
  if NOT keyword_set(OUTPUT) then output = repstr(master, '.fits', '_resolvedSED.fits')
  if NOT keyword_Set(shift) then shift = 0 else shift = 1
  
  mastern = master
  master = mrdfits(master, 1)

  scale = 3.
  
  data   = mrdfits(multiband, 1)
  nfilts = n_elements(data)
  td     = data[0].IMAGE
  s      = size(td, /dim)
  ims    = dblarr(scale*s[0], scale*s[1], nfilts) ;; mag 2x
  ers    = ims
  
  for ii = 0, nfilts - 1 do begin
     ims[*,*,ii] = congrid(data[ii].IMAGE * data[ii].EXPTIME, scale * s[0], scale * s[1], /center)
     ers[*,*,ii] = congrid(sqrt(data[ii].ERR^2 + data[ii].IMAGE) * sqrt(data[ii].EXPTIME), $
                           scale * s[0], scale * s[1], /center)
  endfor
  
  ;; Find the center
  irbs = where(data.FILTER le 160, nirbs)
  sncheck = fltarr(nirbs)
  for ii = 0, nirbs - 1 do $
     sncheck[ii] = mean(ims[15:45,15:45,irbs[ii]]/ers[15:45,15:45,irbs[ii]],/nan)
  irbs = irbs[where(sncheck gt 5. / scale)]
  
  irs  = mean(ims[*,*,irbs], dim = 3, /nan)
  eirs = sqrt(mean(ers[*,*,irbs]^2, dim = 3, /nan))

  s = size(irs, /dim)
  
;  if ~FINITE(total(eirs)) then stop
  
  mwrfits, irs, 'tmp.fits', /create
  spawn, 'sex tmp.fits -c locateObjectOnStamp.sex'
  cat = mrdfits('tmp.cat', 2)
  foo = min(cat.MAG_AUTO, obj)
  cat = cat[obj]
  xx = cat.X_IMAGE - 1 ;; to idl coords
  yy = cat.Y_IMAGE - 1  
  ctr = yy
  y = findgen(s[1])
  tmp = findgen(180) * 2*!pi/179

  mwrfits, congrid(master.DIRECT_IM, s[0], s[1], /center),  'tmp2.fits', /create
;  mwrfits, congrid(master.DIRECT_IM, s[0], s[1]), 'tmp2.fits', /create
  spawn, 'sex tmp2.fits -c locateObjectOnStamp.sex'
  cat2 = mrdfits('tmp.cat', 2)
  foo = min(cat2.MAG_AUTO, obj)
  cat2 = cat2[obj]
  xx2 = cat2.X_IMAGE - 1 ;; to idl coords
  yy2 = cat2.Y_IMAGE - 1  
  ctr2 = yy2
  
  ;; Make the apertures
  innerUpOff = scale * mean(master.INNERUP - master.TRACE) + ctr
  innerDnOff = scale * mean(master.INNERDN - master.TRACE) + ctr
  interUpOff = scale * mean(master.INTERUP - master.TRACE) + ctr
  interDnOff = scale * mean(master.INTERDN - master.TRACE) + ctr
  outerUpOff = scale * mean(master.OUTERUP - master.TRACE) + ctr
  outerDnOff = scale * mean(master.OUTERDN - master.TRACE) + ctr

  innfl = where(y ge INNERDNOFF and y le INNERUPOFF)
  intup = where(y gt INNERUPOFF and y le INTERUPOFF)
  intdn = where(y ge INTERDNOFF and y lt INNERDNOFF)
  outup = where(y gt INTERUPOFF and y le OUTERUPOFF)
  outdn = where(y ge OUTERDNOFF and y lt INTERDNOFF)

  if 0 then begin
     !p.multi = [0,2,0]
     plot, findgen(s[0]), findgen(s[1]), /iso
     cgimage, congrid(master.DIRECT_IM, s[0], s[1], /center), /over, stretch = 2
;     cgimage, congrid(master.DIRECT_IM, s[0], s[1]), /over, stretch = 2
     oplot, xx2 + cos(tmp) * scale * master.RE, yy2 + sin(tmp) * scale * master.RE, thick = 1, col = 'ffa500'x
     oplot, xx2 + cos(tmp) * cat2.FLUX_RADIUS, yy2 + sin(tmp) * cat2.FLUX_RADIUS, thick = 1, col = 'ffa500'x, linesty = 2
     oplot, !X.CRANGE, replicate(innerUpOff, 2) - ctr + ctr2, col = 255
     oplot, !X.CRANGE, replicate(innerDnOff, 2) - ctr + ctr2, col = 255
     oplot, !X.CRANGE, replicate(interUpOff, 2) - ctr + ctr2, col = '00ff00'x
     oplot, !X.CRANGE, replicate(interDnOff, 2) - ctr + ctr2, col = '00ff00'x
     oplot, !X.CRANGE, replicate(outerUpOff, 2) - ctr + ctr2, col = 'ff0000'x
     oplot, !X.CRANGE, replicate(outerDnOff, 2) - ctr + ctr2, col = 'ff0000'x
     
     plot, findgen(s[0]), findgen(s[1]), /iso
     cgimage, irs, /over, stretch = 2
     oplot, xx + cos(tmp) * scale * master.RE, yy + sin(tmp) * scale * master.RE, thick = 1, col = 'ffa500'x
     oplot, xx + cos(tmp) * cat.FLUX_RADIUS, yy + sin(tmp) * cat.FLUX_RADIUS, thick = 1, col = 'ffa500'x, linesty = 2
     oplot, !X.CRANGE, replicate(innerUpOff, 2), col = 255
     oplot, !X.CRANGE, replicate(innerDnOff, 2), col = 255
     oplot, !X.CRANGE, replicate(interUpOff, 2), col = '00ff00'x
     oplot, !X.CRANGE, replicate(interDnOff, 2), col = '00ff00'x
     oplot, !X.CRANGE, replicate(outerUpOff, 2), col = 'ff0000'x
     oplot, !X.CRANGE, replicate(outerDnOff, 2), col = 'ff0000'x
     !p.multi = 0
     k = get_kbrd(1)
  endif
;  stop
  
  ty = replicate(0, s[1])
  ty[innfl] = 100
  ty[intUP] = 51
  ty[intDN] = 50
  ty[outUP] = 26
  ty[outDN] = 25     

  extractIm = reform(ims[*,*,0])
  extractIm[*,*] = 0
  for ii = 0, s[0] - 1 do $
     extractIm[ii,*] = ty
  innfl = where(extractIm eq 100,   nin)
  intup = where(extractIm eq 51, nintup)
  intdn = where(extractIm eq 50, nintdn)
  outUp = where(extractIm eq 26, noutup)
  outDn = where(extractIm eq 25, noutdn)
  inter = where(extractIm eq 50 OR extractIm eq 51, ninter) ;; LEA 20170729
  outer = where(extractIm eq 25 OR extractIm eq 26, nouter)

  
  ;; LEA 20180424 -- calculate rotated INNER, INTER, and
  ;; OUTER contribution to INNER exraction to check LSF calc
  ;; for ref.
  x = findgen(s[0])
  tx = replicate(0, s[1])

  innw = scale * master.NPIX_INNER[0]/2.
  intw = scale * master.NPIX_INTER_UP[0]
  outw = scale * master.NPIX_OUTER_UP[0]
  
  foo = array_indices(irs, innfl)
  tprof = mean(irs[*,reform(foo[1,*])], dim = 2, /nan)
  cprof = total(tprof, /cum)

;  xxx = cat.XPEAK_IMAGE
;  oplot, replicate(xxx, 2), !Y.CRANGE

  ;; 2d centroids don't actually capture the INNER profile
  ;; peak in the case of SSF. Hence, center on peak of the
  ;; collapsed 1D light profile.

  txxx = max(tprof, xxx)

  innInn  = where(x ge (xxx - innW - shift)        AND x le (xxx + innW - shift))
  intInn  = where(x gt (xxx + innW - shift)        AND x le (xxx + innW - shift + intw))
  intInn2 = where(x lt (xxx - innW - shift)        AND x ge (xxx - innw - shift - intw))
  outInn  = where(x gt (xxx + innW - shift + intW) AND x le (xxx - shift + innW + intW + outW))
  outInn2 = where(x lt (xxx - innW - shift - intW) AND x ge (xxx - shift - innW - intW - outW))
 
  tx[innInn]  = 100
  tx[intInn]  = 50
  tx[intInn2] = 50
  tx[outInn]  = 25
  tx[outInn2] = 25
  
;  if (2 * innW MOD 0) then stop
  
  if 0 then begin
     !p.multi = [0,2,0]
     plot, x, tprof / total(tprof)
     oplot, x[innInn], (tprof / total(tprof))[innInn], col = 200
     oplot, x[intInn], (tprof / total(tprof))[intInn], col = '00ff00'x
     oplot, x[intInn2], (tprof / total(tprof))[intInn2], col = '00ff00'x
     oplot, x[outInn], (tprof / total(tprof))[outInn], col = 'ff0000'x
     oplot, x[outInn2], (tprof / total(tprof))[outInn2], col = 'ff0000'x
;     oplot, x, congrid(master.LSF_INNER/scale, scale * 60., /center), col = 255

     plot, x, tx
     k = get_kbrd(1)
     !p.multi = 0
  endif
  
  c1 = total(tprof[innInn])
  c2 = total(tprof[[intInn,intInn2]])
  c3 = total(tprof[[outInn,outInn2]])
  c4 = total(tprof)
  
  extractIm2 = fltarr(s[0],s[1])
  qui = where(y ge INNERDNOFF and y le INNERUPOFF)
  for ii = 0, n_elements(qui) - 1 do $
     extractIm2[*, qui[ii]] = tx
  extractIm2 *= extractIm
  
  if 0 then begin
     !p.multi = [0,2.0]
     plot, findgen(s[0]), findgen(s[1]), /iso
     cgimage, extractIm2 * irs, stretch = 5, /over
     oplot, [xxx], [yy], psym = 1, col = 255

     plot, findgen(s[0]), findgen(s[1]), /iso
     cgimage, extractIm2, stretch = 2, /over
     oplot, [xxx], [yy], psym = 1, col = 255
     !p.multi = 0
     k = get_kbrd(1)
  endif

;  stop
  
  ;; Now get the fluxes/fractions and errors

  tf  = total(irs[where(extractIm eq 100)])
  inninnfl = total(irs[where(extractIm2 eq 100.^2, nInnInn)])
  intinnfl = total(irs[where(extractIm2 eq 50. * 100., nIntInn)])
  outinnfl = total(irs[where(extractIm2 eq 25. * 100., nOutInn)])

  tfe = sqrt(total(eirs[where(extractIm eq 100)]^2))
  einninnfl = sqrt(total(eirs[where(extractIm2 eq 100.^2)]^2))
  eintinnfl = sqrt(total(eirs[where(extractIm2 eq 50. * 100.)]^2))
  eoutinnfl = sqrt(total(eirs[where(extractIm2 eq 25. * 100.)]^2))
  
  savedata = {APERS: [0,0.75,1.5], $
              INN_INN: inninnfl, E_INN_INN: einninnfl, $
              INT_INN: intinnfl, E_INT_INN: eintinnfl, $
              OUT_INN: outinnfl, E_OUT_INN: eoutinnfl, $
              INN_FLUX: tf, E_INN_FLUX: tfe, $
              PROF_INN_INN: c1, PROF_INT_INN: c2, PROF_OUT_INN: c3, $
              PROF_TOT: c4}

  return, savedata
  
;  print, [inninnfl, intinnfl, outinnfl]
;  print, [inninnfl, intinnfl, outinnfl] / total([inninnfl, intinnfl, outinnfl])
;  print, [nInnInn, nIntInn, nOutInn]

end

pro dogood, shift = shift

  if not keyword_set(shift) then shift = 0

  mfiles = ['../MACS0744/00660_2_B.fits',$
            '../MACS1149/00900_1_B.fits',$
            '../MACS1423/01916_2_B.fits',$
            '../MACS2129/00451_2_B.fits']
  files = ['../MACS0744/00660_2_B_allbandImages.fits',$
           '../MACS1149/00900_1_B_allbandImages.fits',$
           '../MACS1423/01916_2_B_allbandImages.fits',$
           '../MACS2129/00451_2_B_allbandImages.fits']
  clusters = 'MACS'+['0744','1149','1423','2129']
  ngals = 4
  names = ['PAS','SSF','PSB','CSF']

  apers  = fltarr(3, ngals)
  fluxes = fltarr(4, ngals)
  cfluxes = fltarr(4, ngals)
  errors = fltarr(4, ngals)

  for ii = 0, ngals - 1 do begin
     foo = getStampLateralContrib(master = mfiles[ii], $
                                  multiband = files[ii], $
                                  cluster = clusters[ii], shift = shift)
     apers[*,ii]  = foo.APERS
     fluxes[*,ii] = [foo.INN_INN, foo.INT_INN, foo.OUT_INN, $
                     foo.INN_FLUX]
     errors[*,ii] = [foo.E_INN_INN, foo.E_INT_INN, foo.E_OUT_INN, $
                     foo.E_INN_FLUX]
     cfluxes[*,ii] = [foo.PROF_INN_INN, foo.PROF_INT_INN, foo.PROF_OUT_INN, $
                      foo.PROF_TOT]
  endfor

  cgloadct, 33, ncolors = ngals, clip = [10,240]

  plotsym, 0, /fill
  plot, findgen(15), findgen(15), /nodat, $
        yran = [0,1], xran = [-0.1,2]
  oplot, !X.CRANGE, [0.5,0.5], linesty = 5
  for ii = 0, ngals - 1 do begin
     oplot, apers[*,ii], fluxes[0:2,ii]/fluxes[-1,ii], $
            col = cgcolor(string(fix(ii)))
     oplot, apers[*,ii], cfluxes[0:2,ii]/cfluxes[-1,ii], $
            col = cgcolor(string(fix(ii))), linesty = 2
     oploterror, apers[*,ii], fluxes[0:2,ii]/fluxes[-1,ii], $
                 sqrt(4 * errors[0:2,ii]^2)/fluxes[0:2,ii], $
                 psym = 8, errcolor = cgcolor(string(fix(ii))), $
                 color = cgcolor(string(fix(ii))), /nohat
  endfor
  legend, /top, /right, box = 0, $
          names, linesty = 0, $
          col = cgcolor(['0','1','2','3']), $
          pspacing = 1

  stop

end
