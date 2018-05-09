function getfluxzonefrac, fitsFile, masterfile

  tdata = mrdfits(masterFile, 1, /silent)
  innApp = tdata.NPIX_INNER[0]
  intApp = innApp + tdata.NPIX_INTER_UP[0] * 2
  outApp = intApp + tdata.NPIX_OUTER_UP[0] * 2
  
  apers = [innApp, intApp, outApp]
  aps = []
  for ii = 0, n_elements(apers) - 1 do $
     aps = [aps, string(apers[ii], f = '(I)')]

  aps = strcompress(aps[0]+','+aps[1]+','+aps[2], /rem)
  
  spawn, 'sex '+fitsFile+' -c locateObjectOnStamp.sex -PARAMETERS_NAME dolateral.param'+$
         ' -CATALOG_NAME tmp.fits -PHOT_APERTURES '+aps+' -BACK_VALUE 0'

  data = mrdfits('tmp.fits', 2)

  foo = max(data.FLUX_AUTO, hit)
  data = data[hit]

  innA = !pi * apers[0]^2/4.
  intA = !pi * (apers[1] - apers[0])^2 / 4.
  outA = !pi * (apers[2] - apers[1])^2 / 4.
  
  innSB = data.FLUX_APER[0] / innA
  intSB = (data.FLUX_APER[1] - data.FLUX_APER[0]) / intA
  outSB = (data.FLUX_APER[2] - data.FLUX_APER[1]) / outA

  exra = innApp * [innApp, $
                   tdata.NPIX_INTER_UP[0] * 2, $
                   tdata.NPIX_OUTER_UP[0] * 2]
  
  fs = [innSB, intSB, outSB] * exra

  aper, (mrdfits(fitsFile, 0, /silent)), 30., 30., $
        flux, eflux, sky, skyerr, 1, apers/2., /flux, /nan, setskyval = 0
  cfs = [flux[0], flux[1] - flux[0], flux[2] - flux[1]] $
        / [innA, intA, outA] * exra
  
  cgloadct, 33
  plot, apers / 2., fs / total(fs), $
        yran = [0,1.1]
  oplot, apers / 2., cfs / total(cfs), col = 'ffa500'x
  oplot, !X.CRANGE, [1,1]

;  stop
  
  ans = {APRAD: apers / 2., FLUX_INNER: fs}
  
  RETURN, ans

end

pro doem

  files = ['00451_2_f160w.fits', $
           '00660_2_f160w.fits', $
           '00900_1_f160w.fits', $
           '01916_2_f160w.fits']
  mfiles = ['../MACS2129/00451_2_B.fits', $
            '../MACS0744/00660_2_B.fits', $
            '../MACS1149/00900_1_B.fits', $
            '../MACS1423/01916_2_B.fits']
  names = ['CSF', 'PAS', 'SSF', 'PSB']

  rads = fltarr(3, n_elements(names))
  res = fltarr(3,n_elements(names))
  for ii = 0, n_elements(files) - 1 do begin
     foo = getfluxzonefrac(files[ii], mfiles[ii])
     print, names[ii], foo.FLUX_INNER / total(foo.FLUX_INNER)
     rads[*,ii] = foo.APRAD
     res[*,ii] = foo.FLUX_INNER / total(foo.FLUX_INNER)
  endfor
  plot, rads[*,0], res[*,0], xran = [0,15], yran = [0,1]
  for ii = 0, 3 do $
     oplot, rads[*,ii], res[*,ii]

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
