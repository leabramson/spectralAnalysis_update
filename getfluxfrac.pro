function getfluxfrac, fitsFile, z, mag

  pixscale = 0.065
  kpcpix = pixscale / zang(1.0, z)
  aper = 2.5 * 2 / kpcpix * sqrt(mag) ;; image plane

  apers = floor(aper + [-2,-1,0,1,2])
  aps = []
  for ii = 0, n_elements(apers) - 1 do $
     aps = [aps, string(apers[ii], f = '(I)')]

  aps = strcompress(aps[0]+','+aps[1]+','+aps[2]+','+aps[3]+','+aps[4], /rem)
  
  spawn, 'sex '+fitsFile+' -c locateObjectOnStamp.sex -PARAMETERS_NAME cutoutff.param'+$
         ' -CATALOG_NAME tmp.fits -PHOT_APERTURES '+aps

  data = mrdfits('tmp.fits', 2)

  foo = max(data.FLUX_AUTO, hit)
  
;  print, fitsFile

  bts = data[hit].FLUX_APER / data[hit].FLUX_AUTO
  
  RETURN, bts

end

;pro doem
;
;  spawn, 'ls *f160w.fits > imfiles.list'
;  readcol, imfiles, files, f = 'A'
;
;  savedata = {ID: string(0), $
;              SBTS: fltarr(5)}
;  savedata = replicate(savedata, n_elements(files))
;
;  for ii = 0, n_elements(files) - 1 do begin
;     savedata[ii].ID = strmid(files[ii], 0, 7)
;     savedata[ii].SBTS[*] = getfluxfrac(files[ii], 1.0)
;  endfor
;  
;end
