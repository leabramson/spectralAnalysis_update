function getuvj, file, z

  ;/home/labramson/Projects/GLASS/spectralAnalysis/Bessel_U.dat
  ;/home/labramson/Projects/GLASS/spectralAnalysis/Bessel_V.dat 
  ;/home/labramson/Projects/GLASS/spectralAnalysis/2MASS_J.dat
  
  c = 2.998e+18
 
  readcol, '~/LocalData/Abramson/FILTER_CURVES/u.fil', lu, ru, f = 'X,F,F', /silent
  readcol, '~/LocalData/Abramson/FILTER_CURVES/v.fil', lv, rv, f = 'X,F,F', /silent
  readcol, '~/LocalData/Abramson/FILTER_CURVES/j.fil', lj, rj, f = 'X,F,F', /silent

  readcol, file, l, f, /silent

  l /= 1+z
  f  = f * l^2 / c

  ru = interpol(ru, lu, l) 
  rv = interpol(rv, lv, l) 
  rj = interpol(rj, lj, l)

  u = total((ru > 0) * f) / total((ru > 0))
  v = total((rv > 0) * f) / total((rv > 0))
  j = total((rj > 0) * f) / total((rj > 0))

  uvj = {UV: -2.5*alog10(u/v), VJ: -2.5*alog10(v/j)}
  
  RETURN, uvj
end

pro test
  
  file1 = 'MACS1423/01916_2_pyspecfitLogNormalResults/OUTER.masked.bestmodel'
  file2 = 'MACS1423/01916_2_pyspecfitPhotResults/OUTER.masked.bestmodel'
  a1 = mrdfits('resultsSummaries/01916_2_lgn_bestfit.fits', 1)
  a2 = mrdfits('resultsSummaries/01916_2_exp_bestfit.fits', 1)

  a1 = a1[where(a1.REGION eq 'OUTER')]
  a2 = a2[where(a2.REGION eq 'OUTER')]
  
  z1 = a1.Z[0]
  z1e = 0.5 * (a1.Z[2] - a1.Z[1])
  z2 = a2.Z[0]
  z2e = 0.5 * (a2.Z[2] - a2.Z[1])

  uv1 = fltarr(20)
  vj1 = fltarr(20)
  uv2 = fltarr(20)
  vj2 = fltarr(20)
  for ii = 0, 19 do begin
     c1 = getuvj(file1, (z1+randomn(seed, 1) * 2*z1e)[0])
     uv1[ii] = c1.UV
     vj1[ii] = c1.VJ
     c2 = getuvj(file2, (z2+randomn(seed, 1) * 2*z2e)[0])
     uv2[ii] = c2.UV
     vj2[ii] = c2.VJ
  endfor

  plot, [0.5,2.5], [0.5,2.5], /iso, /ynoz
  oplot, vj1, uv1, psym = 3
  oplot, vj2, uv2, psym = 1, col = 255
  
end
