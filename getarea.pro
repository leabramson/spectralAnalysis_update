function getarea, fitsfile, z

  objval = -1
  
  z = z[0]
  master = mrdfits(fitsfile, 1)
  image = master.DIRECT_IM
  s = size(image, /dim)

  foo = mrdfits(master.FILENAME, 'sci', head)
  pscale = sxpar(head, 'CD2_2')
  kpcPix = pscale / zang(1.0, z) ;; asec / (asec/kpc)
  apix   = kpcPix^2
  
  mwrfits, image, 'tmp.fits', /create
  spawn, 'sex tmp.fits -c locateObjectOnStamp.sex'+$
         ' -CHECKIMAGE_TYPE segmentation -CHECKIMAGE_NAME tmpSeg.fits'
  cat2 = mrdfits('tmp.cat', 2)
  seg  = readfits('tmpSeg.fits')
  foo = min(cat2.MAG_AUTO, obj)
  seg[where(seg eq obj+1)] = objval
  cat2 = cat2[obj]
  xx2 = cat2.X_IMAGE - 1 ;; to idl coords
  yy2 = cat2.Y_IMAGE - 1  
  ctr = yy2
  
  ;; Make the apertures
  innerUpOff = mean(master.INNERUP - master.TRACE) + ctr
  innerDnOff = mean(master.INNERDN - master.TRACE) + ctr
  interUpOff = mean(master.INTERUP - master.TRACE) + ctr
  interDnOff = mean(master.INTERDN - master.TRACE) + ctr
  outerUpOff = mean(master.OUTERUP - master.TRACE) + ctr
  outerDnOff = mean(master.OUTERDN - master.TRACE) + ctr

  y = findgen(s[1])
  innfl = where(y ge INNERDNOFF and y le INNERUPOFF)
  intup = where(y gt INNERUPOFF and y le INTERUPOFF)
  intdn = where(y ge INTERDNOFF and y lt INNERDNOFF)
  outup = where(y gt INTERUPOFF and y le OUTERUPOFF)
  outdn = where(y ge OUTERDNOFF and y lt INTERDNOFF)

  ty = replicate(0, s[1])
  ty[innfl] = 100
  ty[intUP] = 51
  ty[intDN] = 50
  ty[outUP] = 26
  ty[outDN] = 25     

  extractIm = seg
  for ii = 0, s[0] - 1 do $
     extractIm[ii,*] += ty

  innfl = where(extractIm eq 100 + objval,   nin)
  intup = where(extractIm eq  51 + objval, nintup)
  intdn = where(extractIm eq  50 + objval, nintdn)
  outUp = where(extractIm eq  26 + objval, noutup)
  outDn = where(extractIm eq  25 + objval, noutdn)

  inter = where(extractIm eq 51 + objval OR $    
                extractIm eq 50 + objval, ninter)
  outer = where(extractIm eq 26 + objval OR $    
                extractIm eq 25 + objval, nouter)
  
  npix = [nin, ninter, nouter, $
          nintup, nintdn, noutup, noutdn]
  npix = [npix, total(npix)]
  npix = [npix, total(seg eq objval)]

  regions = ['INNFL', 'INTER', 'OUTER', $
             'INTUP', 'INTDN', 'OUTUP', 'OUTDN', $
             'TOTAL_COVERED', 'TOTAL_DETECTED']
  areas = {REGION:    string(0), $
           AREA_PHYS: 0., $
           AREA_OBS:  0.}
  areas = replicate(areas, n_elements(regions))

  for ii = 0, n_elements(regions) - 1 do begin
     areas[ii].REGION    = regions[ii]
     areas[ii].AREA_PHYS = npix[ii] * apix
     areas[ii].AREA_OBS  = npix[ii] * pscale^2
  endfor
  
  RETURN, areas
  
end
