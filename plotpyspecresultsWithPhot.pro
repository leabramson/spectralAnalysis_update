function grabPyspecResults, filename

  spawn, 'cat '+filename+' | grep "z " > tmpz.dat'
  spawn, 'cat '+filename+' | grep "age" > tmpAge.dat'
  spawn, 'cat '+filename+' | grep "metal" > tmpMetal.dat'
;  spawn, 'cat '+filename+' | grep "H_ew" > tmpEW.dat'
;  if oiii then $
  spawn, 'cat '+filename+' | grep "OIII_ew" > tmpEW_OIII.dat'
  spawn, 'cat '+filename+' | grep "H_ew" > tmpEW.dat'
;  spawn, 'cat '+filename+' | grep "OIIIpsEW" > tmpEW_OIII.dat'
  
  readcol, 'tmpz.dat'    , zFit    , zMin    , zMax    , f = 'X,F,F,F', $
           /silent, /quick, comment = '#'
  readcol, 'tmpAge.dat'  , ageFit  , ageMin  , ageMax  , f = 'X,F,F,F', $
           /silent, /quick, comment = '#'
  readcol, 'tmpMetal.dat', metalFit, metalMin, metalMax, f = 'X,F,F,F', $
           /silent, /quick, comment = '#'
  readcol, 'tmpEW.dat'   , EWFit   , EWMin   , EWMax   , f = 'X,F,F,F', $
           /silent, /quick, comment = '#'

  if file_lines('tmpEW_OIII.dat') gt 0 then $
     readcol, 'tmpEW_OIII.dat', EWFit_OIII, EWMin_OIII, EWMax_OIII, f = 'X,F,F,F', $
              /silent, /quick, comment = '#' $
  else begin
     EWFit_OIII = !VALUES.F_NaN
     EWMin_OIII = !VALUES.F_NaN
     EWMax_OIII = !VALUES.F_NaN
  endelse
     
  savedata = {ZFIT: zFit[0], ZLO: zMin[0], ZHI: zMax[0], $
              AGEFIT: ageFit[0], AGELO: ageMin[0], AGEHI: ageMax[0], $
              METALFIT: metalFit[0], METALLO: metalMin[0], METALHI: metalMax[0], $
              EWFIT: EWFit[0], EWLO: EWmin[0], EWHI: EWmax[0], $
              EWFIT_OIII: EWFit_OIII[0], EWLO_OIII: EWmin_OIII[0], EWHI_OIII: EWmax_OIII[0]}
  
  RETURN, savedata
end

;;
;;
;;

function getlinelist
	 
  
  waves = [3727.092, $
           3729.875, $
           4102.89, $
           4862.68, $ 	
           4960.295, $ 	
           5008.240, $ 	
           6564.61, $ 	
           6718.29, $ 	
           6732.67, $ 	
           3934.777, $ 		 
           3969.588, $ 		 
           4305.61, $ 		 
           5176.7, $
           5895.6];, $
;           6192, $
;           7095]
  names = ['[OII]', $
           '[OII]', $
           'H'+greek('delta'), $   
           'H'+greek('beta'), $   
           '[OIII]', $
           '[OIII]', $
           'H'+greek('alpha'), $   
           '[SII]', $ 
           '[SII]', $
           'CaK', $  
           'CaH', $  
           'G+H'+greek('gamma'), $  
           'Mg', $
           'Na'];, $
;           'TiO', $
;           'TiO'] 
  absorpFlag = [0, $
                0, $
                1, $
                1, $
                0, $
                0, $
                1, $
                0, $
                0, $
                1, $
                1, $
                1, $
                1, $
                1];, $
;                1, $
;                1]
;  replicate(0, n_elements(waves) - 4), $
;                replicate(1, 4)]

  lineList = {WAVES: waves, $
              NAMES: names, $
              ABSORPTION: absorpFlag}
  
  RETURN, lineList
end

;;
;;
;;

pro plotPyspecResults, resultsDir, fitsDir, $
                       DOMASKED = domasked, $
                       MASS = mass, $
                       MAG = mag, $
                       STRETCH = stretch, $
                       CROP = crop, $
                       NEG = neg, $
                       COMBINE = combine, $
                       COLOR = color, $
                       TIFF_NAME = tiff_name, $
                       STAMPSTACK = stampstack, $
                       TIFF_ONFLY = tiff_onfly
;                       HEIGHT = height

  if NOT keyword_set(DOMASKED) then $
     msksfx = '' $
  else $
     msksfx = '.masked'
  if NOT keyword_set(MASS) then mass = -99
  if NOT keyword_set(MAG) then mag = 1.
  if NOT keyword_set(STRETCH) then stretch = 7
;  if NOT keyword_set(HEIGHT) then height = 2. ;; in asec
  if NOT keyword_set(CROP) then crop = 15 ;; in pix
  if NOT keyword_set(NEG) then neg = 0 else neg = 1
  if NOT keyword_set(COMBINE) then combine = 0 else combine = 1
  if NOT keyword_set(COLOR) then color = 0 else color = 1
     
  mass -= alog10(mag)

  objno = strmid(resultsDir, 0, 5)
  pa    = strmid(resultsDir, 6, 1)

  if keyword_set(stampstack) and keyword_set(tiff_onfly) then $
     maketifffromstack, stampstack, objNo+'_'+pa+'_rgb.tiff', /ds9
  
  if color AND NOT keyword_set(TIFF_NAME) then $
     tiff_name = objNo+'_'+pa+'_rgb.tiff'
  
  spawn, 'pwd > foo.dat'
  readcol, 'foo.dat', workingDir, f = 'A', /silent, /quick
  field = strmid(workingDir, 7, /rev)

  obj = field+' '+objno+' PA'+pa
  
  ;; Read results
  regions  = ['OPTFL', $
              'INNFL', $
              'INTUP', $
              'INTDN', $
              'OUTUP', $
              'OUTDN']
  regions += msksfx
  nregions = n_elements(regions)

  ;; Read raw data
  bluDat = mrdfits(fitsDir+'/'+objNo+'_'+pa+'_B.fits', 1, /silent)
  redDat = mrdfits(fitsDir+'/'+objNo+'_'+pa+'_R.fits', 1, /silent)

  zs = fltarr(3, nregions)
  for ii = 0, nregions - 1 do begin

     prefix = resultsDir+'/'+regions[ii]
     res = grabPyspecResults(prefix+'.analyze')

     zs[*,ii] = [res.ZFIT, res.ZLO, res.ZHI]
     
     ;; Read G102
     readcol, prefix+'.bestspec0', $
              lambdaB, innTraceB, innErrB, innPolyB, innModelB, innUsedB, $
              /silent, /quick, comment = '#'
     ;; Read G141
     readcol, prefix+'.bestspec1', $
              lambdaR, innTraceR, innErrR, innPolyR, innModelR, innUsedR, $
              /silent, /quick, comment = '#'

     spec_B  = innTraceB * innPolyB
     err_B   = innErrB   * innPolyB
     model_B = innModelB * innPolyB

     spec_R  = innTraceR * innPolyR
     err_R   = innErrR   * innPolyR
     model_R = innModelR * innPolyR

     spec_B[where(~innUsedB)]  = !values.F_NAN
     err_B[where(~innUsedB)]   = !values.F_NAN
     model_B[where(~innUsedB)] = !values.F_NAN

     spec_R[where(~innUsedR)]  = !values.F_NAN
     err_R[where(~innUsedR)]   = !values.F_NAN
     model_R[where(~innUsedR)] = !values.F_NAN

     lambda = [lambdaB, lambdaR]; / (1 + innerRes.ZFIT)
     flux   = [spec_B, spec_R]
     err    = [err_B, err_R]
     model  = [model_B, model_R]
     
     savedata = {LAMBDA: lambda, $
                 SPEC  : flux  , $
                 ERR   : err   , $
                 MODEL : model }
     results = struct_addtags(savedata, res)          
    
     case ii of
        0: optFlRes = results
        1: innFlRes = results
        2: intUpRes = results
        3: intDnRes = results
        4: outUpRes = results
        5: outDnRes = results
     endcase

  endfor  

  dz = 0.001
  zran = 1.5 * max(zs[2,*]) - 0.5 * min(zs[1,*])
  zran = findgen(zran / dz + 1) * dz
  zpdf = fltarr(n_elements(zran))
  for ii = 1, nregions - 1 do begin
     sig = mean(abs(zs[1:2,ii] - zs[0,ii]))
     pdf = 1. / sqrt(2 * !pi * sig^2) $
           * exp(-0.5 * (zs[0,ii] - zran)^2 / sig^2)
     zpdf += pdf
  endfor
  zcum = total(zpdf * dz, /cum) / total(zpdf * dz)
  endpts = value_locate(zcum, [0.05,.95])
  foo = max(zpdf, zpick)
  zpick = zran[zpick]
    
  bDI = bluDat.DIRECT_IM
  rDI = redDat.DIRECT_IM
  s = size(bDI, /dim)

  ;; Get spatial & wavelegth scales
  foo = mrdfits(bluDat.FILENAME, 'sci', head, /silent)
  pscale = sxpar(head, 'CD2_2') ;; It's downsampled by ~2x from WFC3 plate scale
  lscale = sxpar(head, 'CD1_1')
  pivot  = sxpar(head, 'CRVAL1')
  tpixz  = lscale / pivot * 10 * (1+zpick) ;; 10 pix offset in z space
;  pscale = 0.13 ;; WFC3IR "/pix    
  run = (findgen(s[0]) - (s[0]+1)/2) * pscale
  run += pscale/2

  ;; Redshift sol'n ordering:
  ;  0: optFlRes
  ;  1: innFlRes
  ;  2: intUpRes
  ;  3: intDnRes
  ;  4: outUpRes
  ;  5: outDnRes
  zplt = [zs[0,5],zs[0,3],zs[0,1],zs[0,2],zs[0,4]] 
  crummy = where(abs(zplt - zpick) ge tpixz) 
  
  b2D = bluDat.SPEC2D
  r2D = redDat.SPEC2D
  s2D = size(b2D, /dim)

  blam = bluDat.lambda
  rlam = redDat.lambda

  innOff = median(bluDat.INNERUP - bluDat.TRACE)
  intOff = median(bluDat.INTERUP - bluDat.TRACE)
  outOff = median(bluDat.OUTERUP - bluDat.TRACE)

  kpc       = 1. / zang(1.0, zpick, /silent);optFlRes.ZFIT) ;; kpc / asec
  kpcPix    = kpc * pscale                  ;; kpc / pix
  hlrKpc    = bluDat.RE * kpcPix            ;; re in Kpc
  hlrObs    = bluDat.RE * pscale            ;; re in asec
;  normRadii = [0, mean([innOff,intOff]),mean([intOff,outOff])] / bluDat.RE
  normRadii = [0, mean([bluDat.EXTR_INTDN_OFF, bludat.EXTR_INTUP_OFF]), $
               mean([bluDat.EXTR_OUTDN_OFF, bludat.EXTR_OUTUP_OFF])] / bluDat.RE
  
  bdots    = [-reverse(normRadii[-2:-1]), 0, normRadii[1:2]] * bluDat.RE
  bctroids = [bluDat.LSF_OUTER_DN_PARS[1], $
              bluDat.LSF_INTER_DN_PARS[1], $
              bluDat.LSF_INNER_PARS[1], $
              bluDat.LSF_INTER_UP_PARS[1], $
              bluDat.LSF_OUTER_UP_PARS[1]]
  geoZoff = 10 * (max(bctroids) - min(bctroids)) * lscale / pivot * (1 + zpick) ;; 10x maximum z-error from morpho effects
  
  ;  mass = 10.8 ;; From takahiro + lensing magnification estimated from HFF site @ 2.5x
;  re   = 4.0 * 0.13 / zang(1, 1.086)  ;; In Kpc; from 0753_1_B.fits
;  read_jpeg, 'MACS1149CLS_1096_rgb_rot.jpeg', bdi, true = 1
;  s = size(bdi)
;  w = s[2]
;  tp1 = 20. / 60. * w - 1
;  tp2 = 40. / 60. * w - 1
  
  set_plot, 'PS'
  outputName = objno+'_'+pa+'_withBestfit.eps'
  device, filename = outputName, $
          /col, /encap, /decomp, bits_per_pixel = 8, $
          xsize = 9, ysize = 5./7 * 9., /in

  xs1 = crop
  xs2 = n_elements(run) - crop
  
  scols = [255, '00aa00'x, '005500'x, 'ffa500'x, 'ff5500'x]
  plot, run[xs1:xs2], run[xs1:xs2], /xsty, /ysty, /nodat, $
        xtitle = greek('delta')+'!18x!X ["]', $
        ytitle = greek('delta')+'!18y!X ["]', /iso, $
        pos = [0.1, 0.75, 0.25, 0.90], $
        xthick = 4, ythick = 4, charthick = 3, charsize = 1, $
        xtickint = 0.5, ytickint = 0.5, title = 'direct image', $
        xminor = 1, yminor = 1
  if NOT color then begin
     cgloadct, 0
     cgimage, bdi[xs1:xs2,xs1:xs2], stretch = stretch, /over ;, /neg
  endif else begin     
     cim = read_tiff(tiff_name)
;     mwrfits, bdi, 'tb.fits', /create
;     mwrfits, bdi, 'tg.fits', /create
;     mwrfits, rdi, 'tr.fits', /create
;     spawn, 'stiff -d > stiff.conf'
;     spawn, 'stiff tr.fits tg.fits tb.fits -c stiff.conf -OUTFILE_NAME tmp.tiff'
;     im = read_tiff('tmp.tiff', r, g, b)
;     tvlct, r, g, b
;     tv, im
;     cgloadct, 3
     cgimage, cim[*,xs1:xs2,xs1:xs2], stretch = 2, /over
  endelse
  x = findgen(101) * (2*!pi) / 100.
  oplot, hlrobs * cos(x), hlrobs * sin(x), col = '00a5ff'x
  oplot, bludat.RE_SEX * cos(x) * pscale, bluDat.RE_SEX * sin(x) * pscale, col = '00a5ff'x, linesty = 1
  oplot, run, replicate(innOff , n_elements(run))  * pscale, col = scols[0], linesty = 0, thick = 1.5
  oplot, run, replicate(-innOff, n_elements(run)) * pscale, col = scols[0], linesty = 0, thick = 1.5
  oplot, run, replicate(intOff , n_elements(run))  * pscale, col = scols[1], linesty = 0, thick = 1.5; '00bb00'x
  oplot, run, replicate(-intOff, n_elements(run)) * pscale, col = scols[2], linesty = 0, thick = 1.5; '00bb00'x
  oplot, run, replicate(outOff , n_elements(run))  * pscale, col = scols[3], linesty = 0, thick = 1.5; 'ff5500'x
  oplot, run, replicate(-outOff, n_elements(run)) * pscale, col = scols[4], linesty = 0, thick = 1.5; 'ff5500'x
  plotsym, 0, /fill
  oplot, (bctroids - bctroids[2]) * pscale, bdots * pscale, psym = 1, col = '0055ff'x, $
         symsize = 0.5
  plot, run[xs1:xs2], run[xs1:xs2], /xsty, /ysty, /nodat, $
        pos = [0.1, 0.75, 0.25, 0.90], $
        xtickname = replicate(' ', 60), ytickname = replicate(' ',60), $
        xthick = 4, ythick = 4, charthick = 3, charsize = 1, $
        xtickint = 0.5, ytickint = 0.5, /noer, /iso, $
        yminor = 1, xminor = 1, yticklen = 0.075, xticklen = 0.075
  
  st = !X.WINDOW[1]+0.025

  fullrun = 0.9 - st - 0.05
  bluestretch = max(bluDat.LAMBDA) - min(bluDat.LAMBDA)
  redstgretch = max(redDat.LAMBDA) - min(redDat.LAMBDA)
  fullstretch = max(redDat.LAMBDA) - min(bluDat.LAMBDA)
  g102End = bluestretch / fullstretch
  g141End = 1 - g102End
  g102End *= fullrun
  g141End *= fullrun
  
  plot, blam, findgen(n_elements(b2D[0,xs1:xs2])), /nodat, $ ; / (1 + innerRes.ZFIT)
        xtitle = greek('lambda')+'!Dobs!N ['+texToIDL('\AA')+']', $
        xthick = 4, ythick = 4, charthick = 3, charsize = 1, $
;        pos = [st,!Y.WINDOW[0],0.55,!Y.WINDOW[1]], /noer, $
        pos = [st,!Y.WINDOW[0],st + g102End,!Y.WINDOW[1]], /noer, $
        ytickname = replicate(' ',60), yran = [xs1,xs2], title = obj, /xsty, /ysty
  wid = !X.WINDOW[1] - st
  if keyword_set(DOMASKED) then $
     spec2d = b2D[*,xs1:xs2] * (1 - bluDat.MASK[*,xs1:xs2]) $
  else $
     spec2d = b2D[*,xs1:xs2]
;  ez     = bludat.EXTRACT_ZONES[*,xs1:xs2] * (1 - bluDat.MASK[*,xs1:xs2])
  if neg then begin
     cgimage, spec2d, stretch = stretch, /over, /neg
     cgtext, blam[0] + 200., 0.9 * xs2, /data, align = 0, $
             'G102', charthick = 4, charsize = 1.2
  endif else begin
     cgimage, spec2d, stretch = stretch, /over
     cgtext, blam[0] + 200., 0.9 * xs2, /data, align = 0, $
             'G102', charthick = 4, col = 'ffffff'x, charsize = 1.2
  endelse
  oplot, blam, bludat.innerUP, col = scols[0], linesty = 0, thick = 1.5 ; / (1 + innerRes.ZFIT) 255      
  oplot, blam, bludat.innerDN, col = scols[0], linesty = 0, thick = 1.5; / (1 + innerRes.ZFIT) 255      
  oplot, blam, bludat.interUp, col = scols[1], linesty = 0, thick = 1.5; / (1 + innerRes.ZFIT) '00bb00'x
  oplot, blam, bludat.interDN, col = scols[2], linesty = 0, thick = 1.5; / (1 + innerRes.ZFIT) '00bb00'x
  oplot, blam, bludat.outerUP, col = scols[3], linesty = 0, thick = 1.5; / (1 + innerRes.ZFIT) 'ff5500'x
  oplot, blam, bludat.outerDN, col = scols[4], linesty = 0, thick = 1.5; / (1 + innerRes.ZFIT) 'ff5500'x
  plot, rlam, findgen(n_elements(b2D[0,xs1-1:xs2-1])), /nodat, $ ;  / (1 + innerRes.ZFIT)
        xthick = 4, ythick = 4, charthick = 3, charsize = 1, $
        pos = [st,!Y.WINDOW[0],st + g102End,!Y.WINDOW[1]], /noer, $
        ytickname = replicate(' ',60), xtickname = replicate(' ',60), $
        yran = [xs1,xs2], /xsty, /ysty, xticklen = 0.075, xminor = 2, yminor = 1, ytickint = 0.5/pscale
  
  st = !X.WINDOW[1]+0.025
  plot, rlam, findgen(n_elements(r2D[0,xs1:xs2])), /nodat, $ ;  / (1 + innerRes.ZFIT)
        xtitle = greek('lambda')+'!Dobs!N ['+texToIDL('\AA')+']', $
        xthick = 5, ythick = 5, charthick = 3, charsize = 1, $
        pos = [st,!Y.WINDOW[0], st + g141end,!Y.WINDOW[1]], /noer, $
        ytickname = replicate(' ',60), yran = [xs1,xs2], /xsty, /ysty, $
        title = string(bluDat.RA, f = '(F8.4)')+', '+string(bluDat.DEC, f = '(F8.4)')
  if keyword_set(DOMASKED) then $
     spec2d = r2D[*,xs1:xs2] * (1 - redDat.MASK[*,xs1:xs2]) $
  else $
     spec2d = r2D[*,xs1:xs2]
  if neg then begin
     cgimage, spec2d, stretch = stretch, /over, /neg
     cgtext, rlam[0] + 200., 0.9 * xs2, /data, align = 0, $
             'G141', charthick = 4, charsize = 1.2
  endif else begin
     cgimage, spec2d, stretch = stretch, /over
     cgtext, rlam[0] + 200., 0.9 * xs2, /data, align = 0, $
             'G141', charthick = 4, col = 'ffffff'x, charsize = 1.2
  endelse
  oplot, rlam, reddat.innerUP, col = scols[0], linesty = 0, thick = 1.5; / (1 + innerRes.ZFIT) 255      
  oplot, rlam, reddat.innerDN, col = scols[0], linesty = 0, thick = 1.5; / (1 + innerRes.ZFIT) 255      
  oplot, rlam, reddat.interUp, col = scols[1], linesty = 0, thick = 1.5; / (1 + innerRes.ZFIT) '00bb00'x
  oplot, rlam, reddat.interDN, col = scols[2], linesty = 0, thick = 1.5; / (1 + innerRes.ZFIT) '00bb00'x
  oplot, rlam, reddat.outerUP, col = scols[3], linesty = 0, thick = 1.5; / (1 + innerRes.ZFIT) 'ff5500'x
  oplot, rlam, reddat.outerDN, col = scols[4], linesty = 0, thick = 1.5; / (1 + innerRes.ZFIT) 'ff5500'x
  plot, rlam, findgen(n_elements(r2D[0,xs1:xs2])), /nodat, $ ;  / (1 + innerRes.ZFIT)
        xthick = 4, ythick = 4, charthick = 3, charsize = 1, $
        pos = [st,!Y.WINDOW[0], st + g141end,!Y.WINDOW[1]], /noer, $
        ytickname = replicate(' ',60), xtickname = replicate(' ',60), $
        yran = [xs1,xs2], /xsty, /ysty, xticklen = 0.075, xminor = 2, yminor = 1, ytickint = 0.5/pscale

  if NOT combine then begin
     yr = [0,12]
     yttl = '!18f!X!D'+greek('lambda')+'!N + offset'
  endif else begin
     qui = where(innFlRes.LAMBDA ge 13000 AND innFlRes.LAMBDA le 15000)
     yr = [0, ceil(1+mean(innFLRes.SPEC[qui] / (0.5 * (outUPRes.SPEC[qui] + outDNRes.SPEC[qui])), /nan))]
     yttl = '!18f!X!D'+greek('lambda')+'!N [arb.]' 
  endelse
  plot, lambda, optFlRes.SPEC, $
        xran = [8000,18000] / (1+zpick), yran = yr, xsty = 8+1, $ ;;1.1 * max(optFlRes.SPEC * norms[0])
        xtitle = greek('lambda')+'!Drest!N ['+texToIDL('\AA')+']', $
        ytitle = yttl, /nodat, pos = [0.1,0.1,0.5,0.6], /noer, $
        xthick = 4, ythick = 4, charthick = 3, charsize = 1, yminor = 2, /ysty
  axis, xaxis = 1, xran = [8000,18000], xtitle = greek('lambda')+'!Dobs!N', $
        xthick = 4, ythick = 4, charthick = 3, charsize = 1
  lines = getLinelist()
  for ii = 0, n_elements(lines.WAVES) - 1 do begin
     if optFlRes.SPEC[value_locate(lambda / (1+zpick), lines.WAVES[ii])] gt 0 then begin
        if NOT lines.ABSORPTION[ii] then $
           oplot, replicate(lines.WAVES[ii], 2), !Y.CRANGE, col = '555555'x, linesty = 1 $
        else begin
           oplot, replicate(lines.WAVES[ii], 2), !Y.CRANGE, col = '555555'x
           cgtext, lines.WAVES[ii] - 25, !Y.CRANGE[0] + 0.1, align = 0, /data, $
                   lines.NAMES[ii], charsize = 1, charthick = 3, col = '555555'x, orien = 90
        endelse
     endif
  endfor
  scols = ['777777'x, 255, '00aa00'x, '005500'x, 'ffa500'x, 'ff5500'x]
  mcols = [0,         120, '005500'x, '002200'x, 'ff5500'x, 'ff0000'x]
  offs  = [6, 5, 4, 3, 2, 1]
  if NOT combine then begin
     for ii = 0, nregions - 1 do begin
        case ii of
           0: results = optFlRes 
           1: results = innFlRes 
           2: results = intUpRes 
           3: results = intDnRes
           4: results = outUpRes
           5: results = outDnRes
        endcase
        lambda = results.LAMBDA
        norm = mean(results.SPEC[where(lambda ge 13000 AND lambda le 15000)], /nan)
        oplot, lambda / (1+zpick), smooth(results.SPEC / norm * offs[ii], 2, /nan), col = scols[ii], thick = 3
        oplot, lambda / (1+zpick), smooth(results.MODEL / norm * offs[ii], 2, /nan), col = mcols[ii], thick = 4
;     endif else begin
;        oplot, lambda, smooth(results.SPEC, 2, /nan), col = scols[ii], thick = 3
;        oplot, lambda, smooth(results.MODEL, 2, /nan), col = mcols[ii], thick = 4
                                ;    endelse
     endfor
  endif else begin
;     offs = [2,5,8]
     offs = [1,1,1]; * 2
     for ii = 0, 2 do begin
        case ii of
           2: begin
              results = innFlRes
              spec = results.SPEC
              model = results.MODEL
              tscol = scols[1]
              tmcol = mcols[1]
           end
           1: begin
              r1 = intUPRes
              r2 = intDNRes
              spec = (r1.SPEC/r1.ERR^2 + r2.SPEC/r2.ERR^2) / (1./r1.ERR^2+1./r2.ERR^2)
              model = (r1.MODEL/r1.ERR^2 + r2.MODEL/r2.ERR^2) / (1./r1.ERR^2+1./r2.ERR^2)
              tscol = scols[2]
              tmcol = mcols[2]
           end
           0: begin
              r1 = outUPRes
              r2 = outDNRes
              spec = (r1.SPEC/r1.ERR^2 + r2.SPEC/r2.ERR^2) / (1./r1.ERR^2+1./r2.ERR^2)
              model = (r1.MODEL/r1.ERR^2 + r2.MODEL/r2.ERR^2) / (1./r1.ERR^2+1./r2.ERR^2)
              tscol = scols[4]
              tmcol = mcols[4]
              norm = mean(spec[where(r1.lambda ge 13000 AND r1.lambda le 15000)], /nan)
           end
        endcase
        lambda = results.LAMBDA
        oplot, lambda / (1+zpick), smooth(spec / norm, 2, /nan), col = tscol, thick = 3
        oplot, lambda / (1+zpick), smooth(model / norm, 2, /nan), col = tmcol, thick = 4
     endfor
  endelse
  io = innOff / bluDat.RE
  to = intOff / bluDat.RE
  oo = outOff / bluDat.RE
  legend, /top, /right, box = 0, $
          ['!18R/R!X!De!N!X<'+string(io,f='(F4.2)'), 'fit', $
           string(io,f='(F4.2)')+'<!18R/R!X!De!N!X<'+string(to,f='(F4.2)'), 'fit', $
           string(to,f='(F4.2)')+'<!18R/R!X!De!N<'+string(oo,f='(F4.2)'), 'fit'], $
          col = [255,100,'00a500'x,'004400'x,'ffa500'x,'ff0000'x], $
          pspacing = 1, thick = [2,4,2,4,2,4], $
          charsize = 0.7, charthick = 2, linesty = replicate(0,6)
  
  rgrid = bdots / bluDat.RE
  ageTrend = [outDnRes.AGEFIT, $
              intDnRes.AGEFIT, $
              innFlRes.AGEFIT, $
              intUpRes.AGEFIT, $
              outUPRes.AGEFIT]
  ageHi    = 0.434 * ([outDnRes.AGEHI, $
                       intDnRes.AGEHI, $
                       innFlRes.AGEHI, $
                       intUpRes.AGEHI, $
                       outUPRes.AGEHI] $
                      - ageTrend) / ageTrend
  ageLo    = 0.434 * (ageTrend - $
                      [outDnRes.AGELO, $
                       intDnRes.AGELO, $
                       innFlRes.AGELO, $
                       intUpRes.AGELO, $
                       outUPRes.AGELO]) / ageTrend

  ageTrend = alog10(ageTrend)
  
  if combine then begin

     rgrid = bdots[2:4] / bluDat.RE
     minn = ageTrend[2]
     e1 = 0.5 * (ageHi[1] + ageLo[1])
     e3 = 0.5 * (ageHi[3] + ageLo[3])
     e0 = 0.5 * (ageHi[0] + ageLo[0])
     e4 = 0.5 * (ageHi[4] + ageLo[4])

     mint = (ageTrend[1] / e1^2  + ageTrend[3] / e3^2) / $
            (1. / e1^2 + 1. / e3^2)
     mout = (ageTrend[0] / e0^2  + ageTrend[4] / e4^2) / $
            (1. / e0^2 + 1. / e4^2)

;     (ageTrend[0] / mean(abs([ageLo[0], ageHi[0]]-ageTrend[0]))^2  + $
;             ageTrend[4] / mean(abs([ageLo[4], ageHi[4]]-ageTrend[4]))^2) / $
;            (1. / mean(abs([ageLo[0], ageHi[0]]-ageTrend[0]))^2 + $
;             1. / mean(abs([ageLo[4], ageHi[4]]-ageTrend[4]))^2)

     ageTrend = [minn, mint, mout]
     ageHi = [minn + ageHi[2], mint+mean([e1,e3]), mout+mean([e0,e4])]
     ageLo = [minn - ageLo[2], mint-mean([e1,e3]), mout-mean([e0,e4])]
;     ageHi = [ageHi[2], sqrt(1./(1./ageHi[1]^2 + 1./ageHi[3]^2)), sqrt(1./(1./ageHi[0]^2 + 1./ageHi[4]^2))]
;     ageLo = [ageLo[2], sqrt(1./(1./ageLo[1]^2 + 1./ageLo[3]^2)), sqrt(1./(1./ageLo[0]^2 + 1./ageLo[4]^2))]
     
  endif
  agehiBar    = agehi - ageTrend
  ageloBar    = ageTrend - ageLo

  metalTrend = [outDnRes.METALFIT, $
                intDnRes.METALFIT, $
                innFlRes.METALFIT, $
                intUpRes.METALFIT, $
                outUPRes.METALFIT]
  metalHi    = 0.434 * ([outDnRes.METALHI, $
                         intDnRes.METALHI, $
                         innFlRes.METALHI, $
                         intUpRes.METALHI, $
                         outUPRes.METALHI] $
                        - metalTrend) / metalTrend
  metalLo    = 0.434 * (metalTrend - $
                        [outDnRes.METALLO, $
                         intDnRes.METALLO, $
                         innFlRes.METALLO, $
                         intUpRes.METALLO, $
                         outUPRes.METALLO]) / metalTrend
  
  metalTrend = alog10(metalTrend)

  zsun = 0.02
  if combine then begin

     minn = metalTrend[2]

     e1 = 0.5 * (metalLo[1] + metalHi[1])
     e3 = 0.5 * (metalLo[3] + metalHi[3])
     e0 = 0.5 * (metalLo[0] + metalHi[0])
     e4 = 0.5 * (metalLo[4] + metalHi[4])
     mint = (metalTrend[1] / e1^2  + metalTrend[3] / e3^2) / $
            (1. / e1^2 + 1. / e3^2)
     mout = (metalTrend[0] / e0^2 + metalTrend[4] / e4^2)  / $
            (1. / e0^2 + 1. / e4^2)

     metalTrend = [minn, mint, mout] - alog10(zsun)
     metalHi = [minn + metalHi[2], mint+mean([e1,e3]), mout+mean([e0,e4])] - alog10(zsun)
     metalLo = [minn - metalLo[2], mint-mean([e1,e3]), mout-mean([e0,e4])] - alog10(zsun)

;     metalHi = [metalHi[2], sqrt(1./(1./metalHi[1]^2 + 1./metalHi[3]^2)), sqrt(1./(1./metalHi[0]^2 + 1./metalHi[4]^2))]
;     metalLo = [metalLo[2], sqrt(1./(1./metalLo[1]^2 + 1./metalLo[3]^2)), sqrt(1./(1./metalLo[0]^2 + 1./metalLo[4]^2))]

     ctrdex = 0
     xr = [-0.25, ceil(max(rgrid))]
     scols = [255, '005500'x, 'ff5500'x]     
     
  endif else begin

     ctrdex = 2
     xr = [-1*ceil(max(rgrid)), ceil(max(rgrid))]
     scols = ['ff5500'x, '005500'x, 255, '00aa00'x, 'ffa500'x]

  endelse
  metalhiBar    = metalhi - metalTrend
  metalloBar    = metalTrend - metalLo
  
  plotsym, 0, /fill
  plot, rgrid, ageTrend, /nodat, $
        yran = [min(ageTrend)-0.5, max(ageTrend)+0.5], $
;        ytitle = 'Age [Gyr]', $
        ytickname = replicate(' ',60), $
        xthick = 4, ythick = 4, charthick = 3, $
        pos = [0.55,0.35, 0.875,0.6], $
        /noer, xtickname = replicate(' ', 60), charsize = 1, ysty = 8+1, $
        xran = xr, /xsty
  axis, yaxis = 1, yran = !Y.CRANGE, ythick = 4, charsize = 1, ytitle = 'log Age [Gyr]', $
        charthick = 3, /ysty
  polyfill, [!X.CRANGE[0], !X.CRANGE[0], !X.CRANGE[1], !X.CRANGE[1]], $
            [agelo[ctrdex], agehi[ctrdex], agehi[ctrdex], agelo[ctrdex]], $
            col = 'aaaaaa'x
  oplot, !X.CRANGE, replicate(ageTrend[ctrdex], 2), linesty = 1, col = '777777'x
  polyfill, [rgrid, reverse(rgrid)], $
            [ageHi, reverse(ageLo)], $
            spacing = 0.05, thick = 1, orien = 45
;  off = [-0.1, -0.05, 0, 0.05, 0.1]
  off = replicate(0, nregions - 1)
  for ii = 0, n_elements(ageTrend) - 1 do begin
     if abs(zplt[ii] - zpick) le (geoZoff < tpixz) then $
        plotsym, 0, /fill $
     else $
        plotsym, 0, thick = 3     
     oploterror, [rgrid[ii]] + off[ii], [ageTrend[ii]], [ageHibar[ii]], $
                 /hibar, psym = 8, symsize = 1, $
                 col = scols[ii], errcol = scols[ii], errthick = 3
     oploterror, [rgrid[ii]] + off[ii], [ageTrend[ii]], [ageLobar[ii]], $
                 /lobar, psym = 8, symsize = 1, $
                 col = scols[ii], errcol = scols[ii], errthick = 3
     plotsym, 0, thick = 4
     oplot, [rgrid[ii]] + off[ii], [ageTrend[ii]], symsize = 1.2, psym = 8 
  endfor
  legend, /bottom, /left, box = 0, $ ; pos = [0.17,8.4], /data
          ['!18z!X!DGLASS!N='+string(bluDat.Z, f = '(F5.3)'), $
           '!18z!X!Dfit!N='+string(zpick, f = '(F5.3)'), $
;           '!18z!X!Dfit!N='+string(optFlRes.ZFIT, f = '(F5.3)'), $
           '!18R!X!De!N='+string(hlrObs, f = '(F4.2)')+'"='+string(hlrKpc, f = '(F3.1)')+' kpc'], $
          charsize = 1, charthick = 3

  plot, rgrid, metalTrend, /nodat, $
        yran = [min(metalTrend)-0.3,max(metalTrend)+0.3], $
;        ytitle = 'Age [Gyr]', $
        ytickname = replicate(' ',60), $
        xthick = 4, ythick = 4, charthick = 3, $
        pos = [0.55,0.1,0.875,0.35], $
        /noer, $
        xtitle = '!18R/R!X!De!N', $
        charsize = 1, ysty = 8+1, $
        xran = xr, /xsty
  axis, yaxis = 1, yran = !Y.CRANGE, ythick = 4, charsize = 1, $
        ytitle = 'log !18Z!X [!18Z!X!D'+sunsymbol()+']', $
        charthick = 3, /ysty
  polyfill, [!X.CRANGE[0], !X.CRANGE[0], !X.CRANGE[1], !X.CRANGE[1]], $
            [metallo[ctrdex], metalhi[ctrdex], metalhi[ctrdex], metallo[ctrdex]], $
            col = 'aaaaaa'x
  oplot, !X.CRANGE, replicate(metalTrend[ctrdex], 2), linesty = 1, col = '777777'x
  polyfill, [rgrid, reverse(rgrid)], $
            [metalHi, reverse(metalLo)], $
            spacing = 0.05, thick = 1, orien = 45
  for ii = 0, n_elements(metalTrend) - 1 do begin
     if abs(zplt[ii] - zpick) le (geoZoff < tpixz) then $
        plotsym, 0, /fill $
     else $
        plotsym, 0, thick = 3
     oploterror, [rgrid[ii]] + off[ii], [metalTrend[ii]], [metalHibar[ii]], $
                 /hibar, psym = 8, symsize = 1, $
                 col = scols[ii], errcol = scols[ii], errthick = 3
     oploterror, [rgrid[ii]] + off[ii], [metalTrend[ii]], [metalLobar[ii]], $
                 /lobar, psym = 8, symsize = 1, $
                 col = scols[ii], errcol = scols[ii], errthick = 3
     plotsym, 0, thick = 4
     oplot, [rgrid[ii]] + off[ii], [metalTrend[ii]], symsize = 1.2, psym = 8
;     oplot, replicate(rgrid[ii], 2), metalTrend[ii] + [0,metalHibar[ii]], col = 'ff00ff'x, thick = 5
;     oplot, replicate(rgrid[ii], 2), metalTrend[ii] - [0,metalLobar[ii]], col = 'ff00ff'x, thick = 5
  endfor
  key =   ['m!DF140W!N='+string(bluDat.MAG, f = '(F5.2)')]
  if keyword_set(MASS) then $
     key[0]+='; log!18M!X!D*!N='+string(mass, f = '(F4.1)') $
  else key = key[0]
  legend, /bottom, /left, box = 0, $ ; pos = [0.17,8.4], /data
          key, $
;           '!18R!X!De!N='+string(hlrObs, f = '(F4.2)')+'"='+string(hlrKpc, f = '(F3.1)')+' kpc'], $
          charsize = 1, charthick = 3

  device, /close
  spawn, 'gv '+outputName+' &'

  ;;
  ;;
  ;;

  rgrid = bdots / bluDat.RE
  scols = ['ff5500'x, '005500'x, 255, '00aa00'x, 'ffa500'x]
  
  set_plot, 'PS'
  outputName = objno+'_'+pa+'_emissionAndRedshift.eps'
  device, filename = outputName, $
          /col, /encap, /decomp, bits_per_pixel = 8, $
          xsize = 9, ysize = 9 * 5./7, /in
  xs1 = crop
  xs2 = n_elements(run) - crop
  
  scols = [255, '00aa00'x, '005500'x, 'ffa500'x, 'ff5500'x]
  plot, run[xs1:xs2], run[xs1:xs2], /xsty, /ysty, /nodat, $
        xtitle = greek('delta')+'!18x!X ["]', $
        ytitle = greek('delta')+'!18y!X ["]', /iso, $
        pos = [0.1, 0.75, 0.25, 0.90], $
        xthick = 4, ythick = 4, charthick = 3, charsize = 1, $
        xtickint = 0.5, ytickint = 0.5, title = 'direct image', $
        xminor = 1, yminor = 1
  cgimage, bdi[xs1:xs2,xs1:xs2], stretch = stretch, /over, /neg
  x = findgen(101) * (2*!pi) / 100.
  oplot, hlrobs * cos(x), hlrobs * sin(x), col = '00a5ff'x
  oplot, bludat.RE_SEX * cos(x) * pscale, bluDat.RE_SEX * sin(x) * pscale, col = '00a5ff'x, linesty = 1
  oplot, run, replicate(innOff , n_elements(run))  * pscale, col = scols[0], linesty = 0, thick = 1.5
  oplot, run, replicate(-innOff, n_elements(run)) * pscale, col = scols[0], linesty = 0, thick = 1.5
  oplot, run, replicate(intOff , n_elements(run))  * pscale, col = scols[1], linesty = 0, thick = 1.5; '00bb00'x
  oplot, run, replicate(-intOff, n_elements(run)) * pscale, col = scols[2], linesty = 0, thick = 1.5; '00bb00'x
  oplot, run, replicate(outOff , n_elements(run))  * pscale, col = scols[3], linesty = 0, thick = 1.5; 'ff5500'x
  oplot, run, replicate(-outOff, n_elements(run)) * pscale, col = scols[4], linesty = 0, thick = 1.5; 'ff5500'x
  plotsym, 0, /fill
  oplot, (bctroids - bctroids[2]) * pscale, bdots * pscale, psym = 1, col = '0055ff'x, $
         symsize = 0.5
  plot, run[xs1:xs2], run[xs1:xs2], /xsty, /ysty, /nodat, $
        pos = [0.1, 0.75, 0.25, 0.90], $
        xtickname = replicate(' ', 60), ytickname = replicate(' ',60), $
        xthick = 4, ythick = 4, charthick = 3, charsize = 1, $
        xtickint = 0.5, ytickint = 0.5, /noer, /iso, $
        yminor = 1, xminor = 1, yticklen = 0.075, xticklen = 0.075
  
  st = !X.WINDOW[1]+0.025

  fullrun = 0.9 - st - 0.05
  bluestretch = max(bluDat.LAMBDA) - min(bluDat.LAMBDA)
  redstgretch = max(redDat.LAMBDA) - min(redDat.LAMBDA)
  fullstretch = max(redDat.LAMBDA) - min(bluDat.LAMBDA)
  g102End = bluestretch / fullstretch
  g141End = 1 - g102End
  g102End *= fullrun
  g141End *= fullrun
  
  plot, blam, findgen(n_elements(b2D[0,xs1:xs2])), /nodat, $ ; / (1 + innerRes.ZFIT)
        xtitle = greek('lambda')+'!Dobs!N ['+texToIDL('\AA')+']', $
        xthick = 4, ythick = 4, charthick = 3, charsize = 1, $
;        pos = [st,!Y.WINDOW[0],0.55,!Y.WINDOW[1]], /noer, $
        pos = [st,!Y.WINDOW[0],st + g102End,!Y.WINDOW[1]], /noer, $
        ytickname = replicate(' ',60), yran = [xs1,xs2], title = obj, /xsty, /ysty
  wid = !X.WINDOW[1] - st
  spec2d = bluDat.SPECORIG[*,xs1:xs2]
;  ez     = bludat.EXTRACT_ZONES[*,xs1:xs2] * (1 - bluDat.MASK[*,xs1:xs2])
  if neg then begin
     cgimage, spec2d, stretch = stretch, /over, /neg
     cgtext, blam[0] + 200., 0.9 * xs2, /data, align = 0, $
             'G102', charthick = 4, charsize = 1.2
  endif else begin
     cgimage, spec2d, stretch = stretch, /over
     cgtext, blam[0] + 200., 0.9 * xs2, /data, align = 0, $
             'G102', charthick = 4, col = 'ffffff'x, charsize = 1.2
  endelse
  oplot, blam, bludat.innerUP, col = scols[0], linesty = 0, thick = 1.5 ; / (1 + innerRes.ZFIT) 255      
  oplot, blam, bludat.innerDN, col = scols[0], linesty = 0, thick = 1.5; / (1 + innerRes.ZFIT) 255      
  oplot, blam, bludat.interUp, col = scols[1], linesty = 0, thick = 1.5; / (1 + innerRes.ZFIT) '00bb00'x
  oplot, blam, bludat.interDN, col = scols[2], linesty = 0, thick = 1.5; / (1 + innerRes.ZFIT) '00bb00'x
  oplot, blam, bludat.outerUP, col = scols[3], linesty = 0, thick = 1.5; / (1 + innerRes.ZFIT) 'ff5500'x
  oplot, blam, bludat.outerDN, col = scols[4], linesty = 0, thick = 1.5; / (1 + innerRes.ZFIT) 'ff5500'x
  plot, rlam, findgen(n_elements(b2D[0,xs1-1:xs2-1])), /nodat, $ ;  / (1 + innerRes.ZFIT)
        xthick = 4, ythick = 4, charthick = 3, charsize = 1, $
        pos = [st,!Y.WINDOW[0],st + g102End,!Y.WINDOW[1]], /noer, $
        ytickname = replicate(' ',60), xtickname = replicate(' ',60), $
        yran = [xs1,xs2], /xsty, /ysty, xticklen = 0.075, xminor = 2, yminor = 1, ytickint = 0.5/pscale
  
  st = !X.WINDOW[1]+0.025
  plot, rlam, findgen(n_elements(r2D[0,xs1:xs2])), /nodat, $ ;  / (1 + innerRes.ZFIT)
        xtitle = greek('lambda')+'!Dobs!N ['+texToIDL('\AA')+']', $
        xthick = 5, ythick = 5, charthick = 3, charsize = 1, $
        pos = [st,!Y.WINDOW[0], st + g141end,!Y.WINDOW[1]], /noer, $
        ytickname = replicate(' ',60), yran = [xs1,xs2], /xsty, /ysty
  spec2d = redDat.SPECORIG[*,xs1:xs2]
  if neg then begin
     cgimage, spec2d, stretch = stretch, /over, /neg
     cgtext, rlam[0] + 200., 0.9 * xs2, /data, align = 0, $
             'G141', charthick = 4, charsize = 1.2
  endif else begin
     cgimage, spec2d, stretch = stretch, /over
     cgtext, rlam[0] + 200., 0.9 * xs2, /data, align = 0, $
             'G141', charthick = 4, col = 'ffffff'x, charsize = 1.2
  endelse
  oplot, rlam, reddat.innerUP, col = scols[0], linesty = 0, thick = 1.5; / (1 + innerRes.ZFIT) 255      
  oplot, rlam, reddat.innerDN, col = scols[0], linesty = 0, thick = 1.5; / (1 + innerRes.ZFIT) 255      
  oplot, rlam, reddat.interUp, col = scols[1], linesty = 0, thick = 1.5; / (1 + innerRes.ZFIT) '00bb00'x
  oplot, rlam, reddat.interDN, col = scols[2], linesty = 0, thick = 1.5; / (1 + innerRes.ZFIT) '00bb00'x
  oplot, rlam, reddat.outerUP, col = scols[3], linesty = 0, thick = 1.5; / (1 + innerRes.ZFIT) 'ff5500'x
  oplot, rlam, reddat.outerDN, col = scols[4], linesty = 0, thick = 1.5; / (1 + innerRes.ZFIT) 'ff5500'x
  plot, rlam, findgen(n_elements(r2D[0,xs1:xs2])), /nodat, $ ;  / (1 + innerRes.ZFIT)
        xthick = 4, ythick = 4, charthick = 3, charsize = 1, $
        pos = [st,!Y.WINDOW[0], st + g141end,!Y.WINDOW[1]], /noer, $
        ytickname = replicate(' ',60), xtickname = replicate(' ',60), $
        yran = [xs1,xs2], /xsty, /ysty, xticklen = 0.075, xminor = 2, yminor = 1, ytickint = 0.5/pscale
 
  plot, zran, zpdf * dz / (nregions - 1), $
;        xran = [min(zs[1,*]) - 0.1, max(zs[2,*]) + 0.1], $
        xran = [min(zs[1,*]) - tpixZ, max(zs[2,*]) + tpixZ], $
        xtitle = '!18z!X', ytitle = '!18P!X(!18z!X)', $
        /noer,  pos = [0.1,0.15,0.5,0.65], $
        charsize = 1.2, charthick = 3, xminor = 2, thick = 6, $
        xthick = 4, ythick = 4              ;, xtickint = 0.5
  oplot, replicate(zs[1,0], 2), !Y.CRANGE, thick = 5, linesty = 1, col = '777777'x
  oplot, replicate(zs[2,0], 2), !Y.CRANGE, thick = 5, linesty = 1, col = '777777'x
  scols = ['777777'x, 255, '00aa00'x, '005500'x, 'ffa500'x, 'ff5500'x]
  for ii = 0, nregions - 2 do $
     oplot, replicate(zs[0,ii], 2), !Y.CRANGE, thick = 4, col = scols[ii]
  oplot, replicate(bluDat.Z, 2), !Y.CRANGE, linesty = 5, thick = 3
  oplot, replicate(bluDat.Z_PHOT, 2), !Y.CRANGE, linesty = 2, thick = 3
  oplot, replicate(zpick, 2), !Y.CRANGE, col = '00a5ff'x, thick = 4
  
  EWTrend = [outDNRes.EWFIT, $
             intDnRes.EWFIT, $
             innFLRes.EWFIT, $
             intUPres.EWFIT, $
             outUPRes.EWFIT]
  EWHi    = [outDNRes.EWHI, $
             intDnRes.EWHI, $
             innFLRes.EWHI, $
             intUPres.EWHI, $
             outUPRes.EWHI]
  EWLo    = [outDNRes.EWLO, $
             intDnRes.EWLO, $
             innFLRes.EWLO, $
             intUPres.EWLO, $
             outUPRes.EWLO]
  hiBar    = EWhi - EWTrend
  loBar    = EWTrend - EWLo

  scols = ['ff5500'x, '005500'x, 255, '00aa00'x, 'ffa500'x]
  plot, rgrid, EWTrend, /nodat, $
;        ytitle = 'Age [Gyr]', $
        yran = [0.75 * min(EWtrend), 1.1 * max(EWtrend)], $
        ytickname = replicate(' ',60), $
        xthick = 4, ythick = 4, charthick = 3, $
        pos = [0.55,0.15,0.875,0.40], $
        /noer, charsize = 1, ysty = 8+1, $
        xtitle = '!18R/R!X!De!N', yminor = 2
  axis, yaxis = 1, ythick = 4, charsize = 1, ytitle = 'EW(H'+greek('alpha')+') ['+texToIDL('\AA')+']', $
        charthick = 3, yran = !y.crange, yminor = 2
  polyfill, [!X.CRANGE[0], !X.CRANGE[0], !X.CRANGE[1], !X.CRANGE[1]], $
            ([EWlo[2], EWhi[2], EWhi[2], EWlo[2]] > !Y.CRANGE[0]) < !Y.CRANGE[1], $
            col = 'aaaaaa'x
  oplot, !X.CRANGE, replicate(EWTrend[2], 2), linesty = 1, col = '555555'x
  for ii = 0, n_elements(EWtrend) - 1 do begin
     if abs(zplt[ii] - zpick) le (geoZoff < tpixz) then $
        plotsym, 0, /fill $
     else $
        plotsym, 0, thick = 3
     oploterror, [rgrid[ii]] + off[ii], [EWTrend[ii]], [HiBar[ii]], $
                 /hibar, psym = 8, symsize = 1, $
                 col = scols[ii], errcol = scols[ii], errthick = 3
     oploterror, [rgrid[ii]] + off[ii], [EWTrend[ii]], [LoBar[ii]], $
                 /lobar, psym = 8, symsize = 1, $
                 col = scols[ii], errcol = scols[ii], errthick = 3
     plotsym, 0, thick = 4
     oplot, [rgrid[ii]] + off[ii], [EWTrend[ii]], psym = 8, symsize = 1.2
  endfor
;  legend, /bottom, /left, box = 0, $ ; pos = [0.17,max(EWtrend) + 5], /data
;          'log!18M!X!D*!N='+string(mass, f = '(F4.1)'), $
;          charsize = 1, charthick = 3

  oiiiTrend = [outDNRes.EWFIT_OIII, $
               intDnRes.EWFIT_OIII, $
               innFLRes.EWFIT_OIII, $
               intUPres.EWFIT_OIII, $
               outUPRes.EWFIT_OIII]
  oiiiHi    = [outDNRes.EWHI_OIII, $
               intDnRes.EWHI_OIII, $
               innFLRes.EWHI_OIII, $
               intUPres.EWHI_OIII, $
               outUPRes.EWHI_OIII]
  oiiiLo    = [outDNRes.EWLO_OIII, $
               intDnRes.EWLO_OIII, $
               innFLRes.EWLO_OIII, $
               intUPres.EWLO_OIII, $
               outUPRes.EWLO_OIII]
  oiiihiBar = oiiihi - oiiiTrend
  oiiiloBar = oiiiTrend - oiiiLo
  
  plot, rgrid, oiiiTrend, /nodat, $
        yran = [0.75 * min(oiiiTrend), 1.1 * max(oiiiTrend)], $
;        ytitle = 'Age [Gyr]', $
        ytickname = replicate(' ',60), $
        xthick = 4, ythick = 4, charthick = 3, $
        pos = [0.55,0.40, 0.875,0.65], $
        /noer, xtickname = replicate(' ', 60), charsize = 1, ysty = 8+1, $
        xran = [-1*ceil(max(rgrid)), ceil(max(rgrid))], /xsty
  axis, yaxis = 1, yran = !Y.CRANGE, ythick = 4, charsize = 1, ytitle = 'EW([OIII]) ['+textoidl('\AA')+']', $
        charthick = 3, /ysty
  polyfill, [!X.CRANGE[0], !X.CRANGE[0], !X.CRANGE[1], !X.CRANGE[1]], $
            ([oiiilo[2], oiiihi[2], oiiihi[2], oiiilo[2]] > !Y.CRANGE[0]) < !Y.CRANGE[1], $
            col = 'aaaaaa'x
  oplot, !X.CRANGE, replicate(oiiiTrend[2], 2), linesty = 1, col = '777777'x
;  off = [-0.1, -0.05, 0, 0.05, 0.1]
  off = replicate(0, nregions - 1)
  for ii = 0, n_elements(oiiiTrend) - 1 do begin
     if abs(zplt[ii] - zpick) le (geoZoff < tpixz) then $
        plotsym, 0, /fill $
     else $
        plotsym, 0, thick = 3
     oploterror, [rgrid[ii]] + off[ii], [oiiiTrend[ii]], [oiiiHibar[ii]], $
                 /hibar, psym = 8, symsize = 1, $
                 col = scols[ii], errcol = scols[ii], errthick = 3
     oploterror, [rgrid[ii]] + off[ii], [oiiiTrend[ii]], [oiiiLobar[ii]], $
                 /lobar, psym = 8, symsize = 1, $
                 col = scols[ii], errcol = scols[ii], errthick = 3
     plotsym, 0, thick = 4
     oplot, [rgrid[ii]] + off[ii], [oiiiTrend[ii]], psym = 8, symsize = 1.2
  endfor

  
  device, /close
  spawn, 'gv '+outputName+' &'
  set_plot, 'X'

  
end

pro doAll

  plotPyspecResults, 'M0717',   '399', '1', mass = 11.34, mag = 4.29 ;       0.99   FIT 20160504 
  plotPyspecResults, 'M0717',  '1236', '1', mass = 10.73, mag = 7.28 ;      2.08   FIT 20160504
  plotPyspecResults, 'M1149',   '753', '1', mass = 11.23, mag = 2.48 ;      0.94   FIT 20160504
  plotPyspecResults, 'M2129',   '839', '1', mass = 11.23, mag = 3.81 ;     -1.03   FIT 20160504
  plotPyspecResults, 'M2129',   '841', '1', mass = 10.98, mag = 3.81 ;     -2.64   FIT 20160504
  plotPyspecResults, 'M2129',   '843', '1', mass = 11.09, mag = 3.81 ;      0.19   FIT 20160504
  plotPyspecResults, 'M2129',   '845', '1', mass = 11.11, mag = 3.81 ;     -6.32   FIT 20160504
  plotPyspecResults, 'M2129',  '1681', '1', mass = 11.14, mag = 3.81 ;      1.39   FIT 20160504

  plotPyspecResults, 'M0717',   '399', '2', mass = 11.34, mag = 4.29 ;       0.99   FIT 20160504 
  plotPyspecResults, 'M0717',  '1236', '2', mass = 10.73, mag = 7.28 ;      2.08   FIT 20160504
  plotPyspecResults, 'M1149',   '753', '2', mass = 11.23, mag = 2.48 ;      0.94   FIT 20160504
  plotPyspecResults, 'M2129',   '839', '2', mass = 11.23, mag = 3.81 ;     -1.03   FIT 20160504
  plotPyspecResults, 'M2129',   '841', '2', mass = 10.98, mag = 3.81 ;     -2.64   FIT 20160504
  plotPyspecResults, 'M2129',   '843', '2', mass = 11.09, mag = 3.81 ;      0.19   FIT 20160504
  plotPyspecResults, 'M2129',   '845', '2', mass = 11.11, mag = 3.81 ;     -6.32   FIT 20160504
  plotPyspecResults, 'M2129',  '1681', '2', mass = 11.14, mag = 3.81 ;      1.39   FIT 20160504

  plotPyspecResults, 'M0717_faint',  '454', '1'
  plotPyspecResults, 'M0717_faint',  '1564', '1'
  plotPyspecResults, 'M1149_faint',  '520', '1'
  plotPyspecResults, 'M1149_faint',  '900', '1'
  plotPyspecResults, 'M1149_faint',  '1931', '1';  1.40710  4.00000  21.1797 ** PA 1   fit 20160507
  plotPyspecResults, 'M1423_faint',  '1090', '1'
  plotPyspecResults, 'M2129_faint',  '1050', '1'  
  plotPyspecResults, 'M2129_faint',  '1126', '1' 
  plotPyspecResults, 'R1347_faint',  '1419', '1'

  plotPyspecResults, 'M0717_faint',  '454', '2'
  plotPyspecResults, 'M0717_faint',  '1564', '2'
  plotPyspecResults, 'M1149_faint',  '520', '2'
  plotPyspecResults, 'M1149_faint',  '900', '2'
  plotPyspecResults, 'M1423_faint',  '1090', '2'
  plotPyspecResults, 'M2129_faint',  '1050', '2'  
  plotPyspecResults, 'M2129_faint',  '1126', '2' 
  plotPyspecResults, 'R1347_faint',  '1419', '2'
  
end
