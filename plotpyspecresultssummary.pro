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
pro plotPyspecResultsSummary, resultsDir, fitsDir, $
                              DOMASKED = domasked, $
                              MASS = mass, $
                              MAG = mag, $
                              STRETCH = stretch, $
                              CROP = crop, $
                              NEG = neg, $
                              COMBINE = combine, $
                              COLOR = color, $
                              TIFF_NAME = tiff_name

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
  if keyword_set(DOMASKED) then $
     analyzetofits, resultsDir, objNo+'_'+pa+'_pyspecSummary.fits', $
                    fitsfile = fitsDir+'/'+objNo+'_'+pa+'_B.fits', sfh = 'LogNormal' $
  else $
     analyzetofits, resultsDir, objNo+'_'+pa+'_pyspecSummary.fits', /notmasked, $
                    fitsfile = fitsDir+'/'+objNo+'_'+pa+'_B.fits', sfh = 'LogNormal'
  tres = mrdfits(objNo+'_'+pa+'_pyspecSummary.fits', 1)
  for ii = 0, nregions - 1 do begin

     res = tres[where(tres.REGION eq repstr(regions[ii], '.masked', ''))]

     zs[*,ii] = res.Z

     prefix = resultsDir+'/'+regions[ii]
     
     ;; Read G102
     readcol, prefix+'.bestspec0', $
              lambdaB, innTraceB, innErrB, innPolyB, innModelB, innUsedB, $
              /silent, /quick, comment = '#'
     ;; Read G141
     readcol, prefix+'.bestspec1', $
              lambdaR, innTraceR, innErrR, innPolyR, innModelR, innUsedR, $
              /silent, /quick, comment = '#'

     ;; Read best-fit model
     readcol, prefix+'.bestmodel', l, m, $
              /silent, /quick, comment = '#'
     
     ;; Read photometry
     readcol, prefix+'.bestphot', $
              lambdap, fp, errp, modelp, $
              /silent, /quick, comment = '#'
     
     spec_B  = innTraceB
     err_B   = innErrB  
     model_B = innModelB

     spec_R  = innTraceR
     err_R   = innErrR  
     model_R = innModelR

     spec_B[where(~innUsedB)]  = !values.F_NAN
     err_B[where(~innUsedB)]   = !values.F_NAN

     spec_R[where(~innUsedR)]  = !values.F_NAN
     err_R[where(~innUsedR)]   = !values.F_NAN

     lambda = [lambdaB, lambdaR]
     flux   = [spec_B, spec_R]
     err    = [err_B, err_R]
     model  = [model_B, model_R]

     ;; Interpolate the best model & convolve
     readcol, fitsdir+'/'+ objNo+'_'+pa+'_R_'+regions[ii]+'.lsf', $
              tx, k, /silent, /quick, comment = '#'
     dl = lambda[1] - lambda[0]
     newl = findgen((max(l) - min(l))/dl+1) * dl + min(l)
     newm = interpol(m, l, newl)
     newm = convol(newm, k, /edge_truncate, /center)
     
     savedata = {LAMBDA     : lambda, $
                 SPEC       : flux  , $
                 ERR        : err   , $
                 MODEL      : model , $
                 FULL_LAMBDA: newl, $
                 FULL_MODEL : newm, $
                 PHOT_LAMBDA: lambdap, $
                 PHOTOMETRY : fp, $
                 EPHOT      : errp, $
                 MODELPHOT  : modelp}
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

  ;; Actually plot
  
  set_plot, 'PS'
  outputName = objno+'_'+pa+'_fixedSolar_withBestfit.eps'
  xs = 9
  device, filename = outputName, $
          /col, /encap, /decomp, bits_per_pixel = 8, $
          xsize = xs, ysize = 5./7 * xs, /in

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
     offs = [3,2,1]
     qui = where(innFlRes.LAMBDA ge 13000 AND innFlRes.LAMBDA le 15000)
;     yr = [0, ceil(1+mean(innFLRes.SPEC[qui] / (0.5 *
;     (outUPRes.SPEC[qui] + outDNRes.SPEC[qui])), /nan))]
     tp = outDnRes.PHOTOMETRY * 1d18
     tp = tp[where(tp gt 0)]
     yr = [0.5*min(tp), 10] ; 1.1 * max(innFlres.FULL_MODEL * offs[0])]
;     yttl = '!18f!X!D'+greek('lambda')+'!N [10!E-18!N erg s!E-1!N
;     cm!E-2!N '+textoidl('\AA')+'!E-1!N]'
     yttl = '!18f!X!D'+greek('lambda')+'!N !N [arb.]' 
  endelse
  plot, lambda, optFlRes.SPEC, $
        xran = [2000,18000] / (1+zpick), yran = yr, xsty = 8+1, $ ;;1.1 * max(optFlRes.SPEC * norms[0])
        xtitle = greek('lambda')+'!Drest!N ['+texToIDL('\AA')+']', $
;        xtickname = replicate(' ', 60), $
        pos = [0.1,0.1,!X.WINDOW[1],0.6], $
        ytitle = yttl, /nodat, /noer, $
        xthick = 4, ythick = 4, charthick = 3, charsize = 1, yminor = 2, /ysty, /ylog
  axis, xaxis = 1, xran = [2000.,18000.], xtitle = greek('lambda')+'!Dobs!N ['+textoidl('\AA')+']', $
        xthick = 4, ythick = 4, charthick = 3, charsize = 1, /xsty
  lines = getLinelist()
;  for ii = 0, n_elements(lines.WAVES) - 1 do begin
;     if optFlRes.SPEC[value_locate(lambda / (1+zpick), lines.WAVES[ii])] gt 0 then begin
;        if NOT lines.ABSORPTION[ii] then $
;           oplot, replicate(lines.WAVES[ii], 2), 10.^!Y.CRANGE, col = '555555'x, linesty = 1 $
;        else begin
;           oplot, replicate(lines.WAVES[ii], 2), 10.^!Y.CRANGE, col = '555555'x
;           cgtext, lines.WAVES[ii] - 25, 10.^!Y.CRANGE[0]*1.1, align = 0, /data, $
;                   lines.NAMES[ii], charsize = 1, charthick = 3, col = '555555'x, orien = 90
;        endelse
;     endif
;  endfor
  scols = ['777777'x, 255, '00aa00'x, '005500'x, 'ffa500'x, 'ff5500'x]
  mcols = [0,         120, '005500'x, '002200'x, 'ff5500'x, 'ff0000'x]
;  offs  = [6, 5, 4, 3, 2, 1]
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
        oploterror, lambda / (1+zpick), results.SPEC, results.ERR, $
                    psym = 3, errthick = 4, errcol = scols[ii]
        oploterror, results.LAMBDAP / (1+zpick), results.PHOTOMETRY, results.EPHOT, col = scols[ii], psym = 8, $
                    errcol = scols[ii], errthick = 4
        oplot, results.FULL_LAMBDA / (1+zpick), results.FULL_MODEL, col = mcols[ii], thick = 4
;        oplot, lambda / (1+zpick), smooth(results.SPEC / norm * offs[ii], 2, /nan), col = scols[ii], thick = 3
;        oplot, lambda / (1+zpick), smooth(results.MODEL / norm * offs[ii], 2, /nan), col = mcols[ii], thick = 4
     endfor
  endif else begin
;     offs = [2,5,8]
;     offs = [1,1,1]; * 2
     for ii = 0, 2 do begin
        case ii of
           2: begin
              results = innFlRes
              spec  = results.SPEC
              model = results.FULL_MODEL
              err   = results.ERR
              phot  = results.PHOTOMETRY
              ephot = results.EPHOT
              mphot = results.MODELPHOT
              tscol = scols[1]
              tmcol = mcols[1]
              off = offs[0]
           end
           1: begin
              r1 = intUPRes
              r2 = intDNRes
;              spec = (r1.SPEC/r1.ERR^2 + r2.SPEC/r2.ERR^2) / (1./r1.ERR^2+1./r2.ERR^2)
;              model = (r1.MODEL/r1.ERR^2 + r2.MODEL/r2.ERR^2) / (1./r1.ERR^2+1./r2.ERR^2)
;              err = 1./(1./r1.ERR^2+1./r2.ERR^2)
              spec = 0.5 * (r1.SPEC + r2.SPEC)
              model = 0.5 * (r1.FULL_MODEL + r2.FULL_MODEL)
              err = sqrt(r1.ERR^2 + r2.ERR^2) / sqrt(2.)              
              phot  = 0.5 * (r1.PHOTOMETRY + r2.PHOTOMETRY)
              ephot = sqrt(r1.EPHOT^2 + r2.EPHOT^2) / sqrt(2.)
              mphot = 0.5 * (r1.MODELPHOT + r2.MODELPHOT)
              tscol = scols[2]
              tmcol = mcols[2]
              off = offs[1]
           end
           0: begin
              r1 = outUPRes
              r2 = outDNRes
;              spec = (r1.SPEC/r1.ERR^2 + r2.SPEC/r2.ERR^2) / (1./r1.ERR^2+1./r2.ERR^2)
;              model = (r1.MODEL/r1.ERR^2 + r2.MODEL/r2.ERR^2) / (1./r1.ERR^2+1./r2.ERR^2)
;              err   = 1./(1./r1.ERR^2+1./r2.ERR^2)
              spec = 0.5 * (r1.SPEC + r2.SPEC)
              model = 0.5 * (r1.FULL_MODEL + r2.FULL_MODEL)
              err = sqrt(r1.ERR^2 + r2.ERR^2) / sqrt(2.)
              phot  = 0.5 * (r1.PHOTOMETRY + r2.PHOTOMETRY)
              ephot = sqrt(r1.EPHOT^2 + r2.EPHOT^2) / sqrt(2.)
              mphot = 0.5 * (r1.MODELPHOT + r2.MODELPHOT)
              tscol = scols[4]
              tmcol = mcols[4]
              off = offs[2]
           end
        endcase
        sn = mean(spec/err, /nan)
        lambda = results.LAMBDA
;        oplot, lambda / (1+zpick), spec, err, psym = 3, col = tscol,
;        errcol=tscol, errthick = 4, /nohat
        oplot, lambda / (1+zpick), spec * off, col = tscol, thick = 3
        oploterror, results.PHOT_LAMBDA / (1+zpick), phot * 1d18 * off, ephot * 1d18 * off, $
                    errthick = 2, psym = 8, symsize = 1.35 ;, $
;                    col = tscol, errcol = tscol
        oplot, results.PHOT_LAMBDA / (1+zpick), phot * 1d18 * off, psym = 8, symsize = 0.8, col = tscol
        oplot, results.FULL_LAMBDA / (1+zpick), model * off, col = tmcol, thick = 4
        oplot, results.PHOT_LAMBDA / (1+zpick), mphot * 1d18 * off, col = tmcol, psym = 4, symsize = 2, thick = 2
        cgtext, !X.CRANGE[1] * 0.99, model[value_locate(results.FULL_LAMBDA / (1+zpick), !X.CRANGE[1] * 0.95)] * off * 0.7, $
                /data, '!18<S/N>!X!Dpix!N='+string(sn, f = '(F4.1)'), charsize = 1, charthick = 6, $
                align = 1
        cgtext, !X.CRANGE[1] * 0.99, model[value_locate(results.FULL_LAMBDA / (1+zpick), !X.CRANGE[1] * 0.95)] * off * 0.7, $
                /data, '!18<S/N>!X!Dpix!N='+string(sn, f = '(F4.1)'), col = tscol, charsize = 1, charthick = 2, $
                align = 1
     endfor
  endelse
  io = innOff / bluDat.RE
  to = intOff / bluDat.RE
  oo = outOff / bluDat.RE
  legend, box = 0, /bottom, $;pos = [!x.crange[0], 10.^!y.crange[0]*1.5], /data, $
          ['!18r/r!X!De!N!X<'+string(io,f='(F4.2)'), 'fit', $
           string(io,f='(F4.2)')+'<!18r/r!X!De!N!X<'+string(to,f='(F4.2)'), 'fit', $
           string(to,f='(F4.2)')+'<!18r/r!X!De!N<'+string(oo,f='(F4.2)'), 'fit'], $
          col = [255,100,'00a500'x,'004400'x,'ffa500'x,'ff0000'x], $
          pspacing = 1, thick = [3,4,3,4,3,4], $
          charsize = 1, charthick = 2, linesty = replicate(0,6), /horiz, spacing = 0.2
  mstel    = alog10(total(10.^tres[where(tres.REGION ne 'OPTFL')].LMASS[0]))
  mstelerr = sqrt(total(0.5 * (tres[where(tres.REGION ne 'OPTFL')].LMASS[2] - tres[where(tres.REGION ne 'OPTFL')].LMASS[1])^2)) 
  legend, /top, /left, box = 0, $ ; pos = [0.17,8.4], /data
          [$                       ;'!18z!X!DGLASS!N='+string(bluDat.Z, f = '(F5.3)'), $
          '!18z!X='+string(zpick, f = '(F5.3)')+textoidl('\pm')+string(0.5*(zs[2,0] - zs[1,0]), f = '(F5.3)'), $
;           '!18z!X!Dfit!N='+string(optFlRes.ZFIT, f = '(F5.3)'), $
          '!18r!X!De!N='+string(hlrObs, f = '(F4.2)')+'"='+string(hlrKpc, f = '(F3.1)')+' kpc', $
          'F140W='+string(bluDat.MAG, f = '(F5.2)'), $
          'log!18M!X!D*!N='+string(mstel, f = '(F4.1)')+textoidl('\pm')+string(mstelerr, f = '(F4.2)')], $
          charsize = 1, charthick = 3

  device, /close
  spawn, 'gv '+outputName+' &'
  set_plot, 'X'

;  stop
  
end

pro test

  plotPyspecResultsSummary, $
     '00900_1_pyspecfitLogNormalResults/', './', $
     /domasked, /combine, /color, tiff_name = '00900_1_rgb.tiff'

  plotPyspecResultsSummary, $
     '00660_2_pyspecfitLogNormalResults/', './', $
     /domasked, /combine, /color, tiff_name = '00660_2_rgb.tiff'

  plotPyspecResultsSummary, $
     '00451_2_pyspecfitLogNormalResults/', './', $
     /domasked, /combine, /color, tiff_name = '00451_2_rgb.tiff'
  
end
