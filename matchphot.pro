;; Actually extract photometry from the stamps made by getsrcimage.pro

pro matchphot, master = master, $
               multiband = multiband, $
               zeropoint = zeropoint, $
               cluster = cluster, $
               dustfile = dustfile, $
               OUTPUT = output

  if NOT keyword_set(zeropoint) then zeropoint = 23.9 ;; AB mags
  if NOT keyword_set(dustfile) then dustfile = cluster+'/GALACTIC_EXTINCTION.dat'
  if NOT keyword_set(OUTPUT) then output = repstr(master, '.fits', '_resolvedSED.fits')
  
  mastern = master
  master = mrdfits(master, 1)

  data   = mrdfits(multiband, 1)
  nfilts = n_elements(data)
  td     = data[0].IMAGE
  s      = size(td, /dim)
  ims    = dblarr(s[0], s[1], nfilts)

  for ii = 0, nfilts - 1 do $
     ims[*,*,ii] = data[ii].IMAGE
  
  ;; Find the center
  irbs = where(data.FILTER le 160)
  irs  = mean(ims[*,*,irbs], dim = 3)
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

  mwrfits, master.DIRECT_IM, 'tmp2.fits', /create
  spawn, 'sex tmp2.fits -c locateObjectOnStamp.sex'
  cat2 = mrdfits('tmp.cat', 2)
  foo = min(cat2.MAG_AUTO, obj)
  cat2 = cat2[obj]
  xx2 = cat2.X_IMAGE - 1 ;; to idl coords
  yy2 = cat2.Y_IMAGE - 1  
  ctr2 = yy2
  
  ;; Make the apertures
  innerUpOff = mean(master.INNERUP - master.TRACE) + ctr
  innerDnOff = mean(master.INNERDN - master.TRACE) + ctr
  interUpOff = mean(master.INTERUP - master.TRACE) + ctr
  interDnOff = mean(master.INTERDN - master.TRACE) + ctr
  outerUpOff = mean(master.OUTERUP - master.TRACE) + ctr
  outerDnOff = mean(master.OUTERDN - master.TRACE) + ctr

  fluxes = fltarr(5, n_elements(data))
  errs   = fltarr(5, n_elements(data))

  innfl = where(y ge INNERDNOFF and y le INNERUPOFF)
  intup = where(y gt INNERUPOFF and y le INTERUPOFF)
  intdn = where(y ge INTERDNOFF and y lt INNERDNOFF)
  outup = where(y gt INTERUPOFF and y le OUTERUPOFF)
  outdn = where(y ge OUTERDNOFF and y lt INTERDNOFF)

  !p.multi = [0,2,0]
  plot, findgen(s[0]), findgen(s[1]), /iso
  cgimage, master.DIRECT_IM, /over, stretch = 2
  oplot, xx2 + cos(tmp) * master.RE, yy2 + sin(tmp) * master.RE, thick = 1, col = 'ffa500'x
  oplot, xx2 + cos(tmp) * cat2.FLUX_RADIUS, yy2 + sin(tmp) * cat2.FLUX_RADIUS, thick = 1, col = 'ffa500'x, linesty = 2
  oplot, !X.CRANGE, replicate(innerUpOff, 2) - ctr + ctr2, col = 255
  oplot, !X.CRANGE, replicate(innerDnOff, 2) - ctr + ctr2, col = 255
  oplot, !X.CRANGE, replicate(interUpOff, 2) - ctr + ctr2, col = '00ff00'x
  oplot, !X.CRANGE, replicate(interDnOff, 2) - ctr + ctr2, col = '00ff00'x
  oplot, !X.CRANGE, replicate(outerUpOff, 2) - ctr + ctr2, col = 'ff0000'x
  oplot, !X.CRANGE, replicate(outerDnOff, 2) - ctr + ctr2, col = 'ff0000'x

  plot, findgen(s[0]), findgen(s[1]), /iso
  cgimage, irs, /over, stretch = 2
  oplot, xx + cos(tmp) * master.RE, yy + sin(tmp) * master.RE, thick = 1, col = 'ffa500'x
  oplot, xx + cos(tmp) * cat.FLUX_RADIUS, yy + sin(tmp) * cat.FLUX_RADIUS, thick = 1, col = 'ffa500'x, linesty = 2
  oplot, !X.CRANGE, replicate(innerUpOff, 2), col = 255
  oplot, !X.CRANGE, replicate(innerDnOff, 2), col = 255
  oplot, !X.CRANGE, replicate(interUpOff, 2), col = '00ff00'x
  oplot, !X.CRANGE, replicate(interDnOff, 2), col = '00ff00'x
  oplot, !X.CRANGE, replicate(outerUpOff, 2), col = 'ff0000'x
  oplot, !X.CRANGE, replicate(outerDnOff, 2), col = 'ff0000'x
  !p.multi = 0
  
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

  ; data.SPEC2D - fan(smooth(data.F_OPTI, 100, /edge_truncate), 60) * data.MODEL_TRACE

  prof = master.MODEL_TRACE[n_elements(master.LAMBDA)/2-1, *] ;; Spatial profile
  foo  = max(prof, hit)
  tx   = findgen(n_elements(prof)) - hit + ctr  ;; For weighting the phototmetry in the optimal extraction
  foo = interpol(prof, findgen(n_elements(prof)), tx)  

  twim = fltarr(s[0],s[1])
  for jj = 0, s[0] - 1 do $
     twim[jj,*] = prof

  owts = [1., $
          mean(twim[intDN]) / mean(twim[innfl]), $
          mean(twim[intUP]) / mean(twim[innfl]), $
          mean(twim[outDN]) / mean(twim[innfl]), $
          mean(twim[outUP]) / mean(twim[innfl])]

  optiflux = fltarr(nfilts)
  optifler = fltarr(nfilts)
  
  for ii = 0, nfilts - 1 do begin
     tim = data[ii].IMAGE
     ter = data[ii].ERR

     fluxes[*,ii] = [total(tim[innfl], /nan), $
                     total(tim[intDN], /nan), $
                     total(tim[intUp], /nan), $
                     total(tim[outDn], /nan), $
                     total(tim[outUp], /nan)]
     errs[*,ii] = sqrt([total(ter[innfl]^2, /nan), $
                        total(ter[intDN]^2, /nan), $
                        total(ter[intUP]^2, /nan), $
                        total(ter[outDN]^2, /nan), $
                        total(ter[outUP]^2, /nan)])          

     optiflux[ii] = total(fluxes[*,ii] * owts)
     optifler[ii] = sqrt(total(errs[*,ii] ^2 * owts))
     
;     print, fluxes[*,ii]/errs[*,ii]
;     crap = where(fluxes[*,ii] gt 0 AND fluxes[*,ii] lt 1d5 AND errs[*,ii] gt 1d5)
;     errs[crap,ii] = sqrt(fluxes[crap,ii])

     lims = where(fluxes[*,ii] le 2 * errs[*,ii], compl = hits, nlims) ;; 2 sigma non-detection floor
     if nlims gt 0 then $
        fluxes[lims,ii] = 0
    
;     errs[hits,ii] = errs[hits,ii] > (0.03 * fluxes[hits,ii]) ;;
;     floor at 3% -- MULTINEST is telling me this was too inflated of
;                    an error bar...
;     errs[*,ii] /= fluxes[*,ii]
;     errs[*,ii] = abs([mean(ter[innfl]/tim[innfl]) / sqrt(nin), $
;                       mean(ter[intDN]/tim[intDN]) / sqrt(nintdn), $
;                       mean(ter[intUp]/tim[intUP]) / sqrt(nintup), $
;                       mean(ter[outDn]/tim[outDN]) / sqrt(noutdn), $
;                       mean(ter[outUp]/tim[outUP]) / sqrt(noutup)])
  endfor

  lims = where(optiflux le 2 * optifler, compl = hits, nlims)
  if nlims gt 0 then $
     optiflux[lims] = 0
  
  filters = data.FILTER
  filters[irbs] *= 10
  fs = sort(filters)
  
  mags = -2.5 * alog10(fluxes) + 25
  newfluxes = 10.^((zeropoint - mags)/2.5)

  omag = -2.5 * alog10(optiflux) + 25
  newoflux = 10.^((zeropoint - omag)/2.5)
  
  cgloadct, 33
  window, 1
  plot, filters[fs], newfluxes[0,fs], /nodat, $
        xran = [300,1800], yran = [0,1.1*max(newfluxes)]
  oploterror, filters[fs], newfluxes[0,fs], errs[0,fs], $; * newfluxes[0,fs], $
              col = '0000ff'x, errcol = '0000ff'x, linesty = 0
  oploterror, filters[fs], newfluxes[1,fs], errs[1,fs], $;  * newfluxes[1,fs], $
              col = '00aa00'x, errcol = '00aa00'x, linesty = 0
  oploterror, filters[fs], newfluxes[2,fs], errs[2,fs], $;  * newfluxes[2,fs], $
              col = '00ff00'x, errcol = '00ff00'x, linesty = 0
  oploterror, filters[fs], newfluxes[3,fs], errs[3,fs], $;  * newfluxes[3,fs], $
              col = 'ff0000'x, errcol = 'ff0000'x, linesty = 0
  oploterror, filters[fs], newfluxes[4,fs], errs[4,fs], $;  * newfluxes[4,fs], $
              col = 'ffa500'x, errcol = 'ffa500'x, linesty = 0
  oploterror, filters[fs], newoflux[fs], optifler[fs], $
              col = '777777'x, errcol = '777777'x
  
  ;; Correct for dust extinction
  readcol, dustfile, bpasses, av, f = 'I,F', $
           /silent, comment = '#', /quick

  for ii = 0, nfilts - 1 do begin
     hit = where(bpasses eq filters[ii], nhit)
     if nhit ne 1 then stop
     tav = av[hit[0]]
     newfluxes[*,ii] *= 10.^(tav/2.5)
     newoflux[ii] *= 10.^(tav/2.5)
  endfor

  oploterror, filters[fs], newfluxes[0,fs], errs[0,fs], $ ;*newfluxes[0,fs], $
              col = '0000ff'x, errcol = '0000ff'x, linesty = 2
  oploterror, filters[fs], newfluxes[1,fs], errs[1,fs], $ ;*newfluxes[1,fs], $
              col = '00aa00'x, errcol = '00aa00'x, linesty = 2
  oploterror, filters[fs], newfluxes[2,fs], errs[2,fs], $ ;*newfluxes[2,fs], $
              col = '00ff00'x, errcol = '00ff00'x, linesty = 2
  oploterror, filters[fs], newfluxes[3,fs], errs[3,fs], $ ;*newfluxes[3,fs], $
              col = 'ff0000'x, errcol = 'ff0000'x, linesty = 2
  oploterror, filters[fs], newfluxes[4,fs], errs[4,fs], $ ;*newfluxes[4,fs], $
              col = 'ffa500'x, errcol = 'ffa500'x, linesty = 2
  oploterror, filters[fs], newoflux[fs], optifler[fs], $
              col = '777777'x, errcol = '777777'x, linesty = 2
  
  savedata = {FILTER:           0, $
              SED_OPTI:         0., $
              SED_OPTI_ERR:     0., $
              SED_INNER:        0., $
              SED_INNER_ERR:    0., $
              SED_INTER_DN:     0., $       
              SED_INTER_DN_ERR: 0., $
              SED_INTER_UP:     0., $       
              SED_INTER_UP_ERR: 0., $
              SED_OUTER_DN:     0., $       
              SED_OUTER_DN_ERR: 0., $
              SED_OUTER_UP:     0., $       
              SED_OUTER_UP_ERR: 0.}
  savedata = replicate(savedata, nfilts)
  
  for ii = 0, nfilts - 1 do begin
     jj = fs[ii]
     savedata[ii].FILTER           = filters[jj]
     savedata[ii].SED_OPTI         = newoflux[jj]
     savedata[ii].SED_OPTI_ERR     = optifler[jj]     
     savedata[ii].SED_INNER        = newfluxes[0,jj]
     savedata[ii].SED_INNER_ERR    = errs[0,jj]     
     savedata[ii].SED_INTER_DN     = newfluxes[1,jj]
     savedata[ii].SED_INTER_DN_ERR = errs[1,jj]
     savedata[ii].SED_INTER_UP     = newfluxes[2,jj]
     savedata[ii].SED_INTER_UP_ERR = errs[2,jj]
     savedata[ii].SED_OUTER_DN     = newfluxes[3,jj]
     savedata[ii].SED_OUTER_DN_ERR = errs[3,jj]
     savedata[ii].SED_OUTER_UP     = newfluxes[4,jj]
     savedata[ii].SED_OUTER_UP_ERR = errs[4,jj]     
  endfor
  
  mwrfits, savedata, output, /create
  
end
;matchPhot, master = 'ABEL2744/00793_1_B.fits', multiband = 'ABEL2744/00793_1_B_allbands.fits', cluster = 'ABEL2744'
;matchPhot, master = 'MACS1149/00753_1_B.fits', multiband = 'MACS1149/00753_1_B_allbands.fits', cluster = 'MACS1149'
;matchPhot, master = 'MACS0744/00660_2_B.fits', multiband = 'MACS0744/00660_2_B_allbands.fits', cluster = 'MACS0744'
;matchPhot, master = 'MACS1149/00900_1_B.fits', multiband = 'MACS1149/00900_1_B_allbands.fits', cluster = 'MACS1149'
