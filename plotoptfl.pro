function getlinelist
  
  waves = [3727.092, $
;           3729.875, $
           4102.89, $
           4862.68, $ 	
;           4960.295, $ 	
           5008.240, $ 	
           6564.61, $ 	
;           6718.29, $ 	
           6732.67, $ 	
;           3934.777, $ 		 
           3969.588, $ 		 
           4305.61, $ 		 
           5176.7, $
           5895.6];, $
;           6192, $
;           7095]
  names = ['[OII]', $
;           '[OII]', $
           'H'+greek('delta'), $   
           'H'+greek('beta'), $   
;           '[OIII]', $
           '[OIII]', $
           'H'+greek('alpha'), $   
;           '[SII]', $ 
           '[SII]', $
;           'CaK', $  
           'CaH+K', $  
           'G+H'+greek('gamma'), $  
           'Mg', $
           'Na'];, $
;           'TiO', $
;           'TiO'] 
  absorpFlag = [0, $
;                0, $
                1, $
                1, $
;                0, $
                0, $
                1, $
;                0, $
                0, $
;                1, $
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
pro plotOptfl, resultsDir, fitsDir, $
               resultsDir2 = resultsdir2, $
               DOMASKED = domasked, $
               MASS = mass, $
               MAG = mag, $
               EMAG = emag, $
               STRETCH = stretch, $
               CROP = crop, $
               NEG = neg, $
               COLOR = color, $
               TIFF_NAME = tiff_name, $
               CLUSTER = cluster, $
               OUTPUT = output

  if NOT keyword_set(DOMASKED) then $
     msksfx = '' $
  else $
     msksfx = '.masked'
  if NOT keyword_set(MASS) then mass = -99
  if NOT keyword_set(MAG) then mag = 1.
  if NOT keyword_set(EMAG) then emag = 0
  if NOT keyword_set(STRETCH) then stretch = 7
;  if NOT keyword_set(HEIGHT) then height = 2. ;; in asec
  if NOT keyword_set(CROP) then crop = 15 ;; in pix
  if NOT keyword_set(NEG) then neg = 0 else neg = 1
  if NOT keyword_set(COMBINE) then combine = 0 else combine = 1
  if NOT keyword_set(COLOR) then color = 0 else color = 1

  print, mag, 1./alog(10) * emag / mag
  
  mass -= alog10(mag)

  objno = strmid(resultsDir, 0, 5)
  pa    = strmid(resultsDir, 6, 1)

  if color AND NOT keyword_set(TIFF_NAME) then $
     tiff_name = objNo+'_'+pa+'_rgb.tiff'
  
  spawn, 'pwd > foo.dat'
  readcol, 'foo.dat', workingDir, f = 'A', /silent, /quick
  field = strmid(workingDir, 7, /rev)

  obj = field+' '+objno+'_'+pa
  
  ;; Read results
  regions  = ['OPTFL']
  regions += msksfx
  nregions = n_elements(regions)

  ;; Read raw data
  bluDat = mrdfits(fitsDir+'/'+objNo+'_'+pa+'_B.fits', 1, /silent)
  redDat = mrdfits(fitsDir+'/'+objNo+'_'+pa+'_R.fits', 1, /silent)

  width = median(bluDat.outerup - bluDat.outerdn)
  outer = median(bluDat.outerup - bluDat.interUp)
  inter = median(bluDat.interup - bluDat.innerUp)
  inner = median(bluDat.innerup - bluDat.innerDn)
    
  zs = fltarr(3, nregions)
  if keyword_set(DOMASKED) then $
     analyzetofits, resultsDir, objNo+'_'+pa+'_pyspecSummary.fits', $
                    fitsfile = fitsDir+'/'+objNo+'_'+pa+'_B.fits', sfh = 'LogNormal' $
  else $
     analyzetofits, resultsDir, objNo+'_'+pa+'_pyspecSummary.fits', /notmasked, $
                    fitsfile = fitsDir+'/'+objNo+'_'+pa+'_B.fits', sfh = 'LogNormal'
  tres = mrdfits(objNo+'_'+pa+'_pyspecSummary.fits', 1)
  
  for ii = 0, nregions - 1 do begin

     treg = repstr(regions[ii], '.masked', '')
     res = tres[where(tres.REGION eq treg)]

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
          
     if keyword_set(resultsDir2) then begin
        prefix2 = resultsDir2+'/'+regions[ii]
        readcol, prefix2+'.bestmodel', l2, m2, $
                 /silent, /quick, comment = '#'
        newm2 = interpol(m2, l2, newl)
        newm2 = convol(newm2, k, /edge_truncate, /center) 

        sd2 = {OTHER_MODEL: newm2}
        results = struct_addtags(results, sd2)
     endif
     
     case treg of
        'OPTFL': optFlRes = results
        'INNFL': innFlRes = results
        'INTER': interRes = results
        'OUTER': outerRes = results
        'INTUP': intUpRes = results
        'INTDN': intDnRes = results
        'OUTUP': outUpRes = results
        'OUTDN': outDnRes = results
     endcase

  endfor  
    
  bDI = bluDat.DIRECT_IM
  rDI = redDat.DIRECT_IM
  s = size(bDI, /dim)

  zplt = zs[0,0]
  zpick = zplt

  ;; Get spatial & wavelegth scales
  foo = mrdfits(bluDat.FILENAME, 'sci', head, /silent)
  pscale = sxpar(head, 'CD2_2') ;; It's downsampled by ~2x from WFC3 plate scale
  lscale = sxpar(head, 'CD1_1')
  pivot  = sxpar(head, 'CRVAL1')
  tpixz  = lscale / pivot * 10 * (1+zpick) ;; 10 pix offset in z space
  run = (findgen(s[0]) - (s[0]+1)/2) * pscale
  run += pscale/2
  
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

  stats = '!18z!X=!N'+string(zpick, f = '(F4.2)')+$
          '!C!18r!X!De!N='+string(hlrObs, f = '(F3.1)')+'"='+$
          string(hlrkpc, f = '(F3.1)')+' kpc'

  ;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;
  ;; Actually plot ;;
  ;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;

  set_plot, 'PS'
  if NOT keyword_set(OUTPUT) then $
     outputName = objno+'_'+pa+'_fixedSolar_Optfl.eps' $
  else $
     outputName = output
  xs = 10.
  device, filename = outputName, $
          /col, /encap, /decomp, bits_per_pixel = 8, $
          xsize = xs, ysize = 3, /in

  ar = 3./xs
  
  xstop = 0.95
  xs1 = crop
  xs2 = n_elements(run) - crop
  
  ;;;;;;;;;;;;;;;;;;
  ;; DIRECT IMAGE ;;
  ;;;;;;;;;;;;;;;;;;

  plot, findgen(s[0]), findgen(s[1]), /xsty, /ysty, /nodat, $
;        xtitle = greek('delta')+'!18x!X ["]', $
;        ytitle = greek('delta')+'!18y!X ["]', /iso, $
        xtickname = replicate(' ', 60), $
        ytickname = replicate(' ', 60), $
        pos = [0.08, 0.25, 0.08 + 0.6 * ar, 0.85], $
        xthick = 4, ythick = 4, charthick = 4, charsize = 1.25, $
        xtickint = 1, ytickint = 1, $;, title = 'direct image', $
        xminor = 1, yminor = 1, title = cluster+' '+objNo+'_'+pa
  if NOT color then begin
     cgloadct, 0
     cgimage, bdi, stretch = stretch, /over ;, /neg
  endif else begin     
     c2 = read_tiff(tiff_name, orient = tiffp)
     cim = c2
     for rr = 0, 2 do $
        cim[rr,*,*] = rotate(reform(c2[rr,*,*]), 7)
     cgimage, cim, stretch = 2, /over
  endelse
  yctr = s[1]/2
  x = replicate(55,2)
;  x = findgen(101) * (2*!pi) / 100.
;  if NOT color then $
;     rcol = '00a5ff'x $
;  else $
;     rcol = 'ffffff'x ;'e1cf9d'x
;  oplot, hlrobs * cos(x), hlrobs * sin(x), col = rcol
;  oplot, bludat.RE_SEX * cos(x) * pscale, bluDat.RE_SEX * sin(x) * pscale, col = rcol, linesty = 1
;  oplot, (bctroids - bctroids[2]) * pscale, bdots * pscale, psym = 1, col = rcol, $
;         symsize = 0.5
  oplot, replicate(10,2), 5+[0, hlrObs / pscale], thick = 8, col = 'ffffff'x
  oplot, replicate(15,2), 5+[0, 5/kpcPix], thick = 8, col = 'ffffff'x
  oplot, x, yctr+width*[-0.5,0.5], thick = 8, col = 'ffffff'x
  oplot, x, yctr-width/2.+[0,outer], thick = 8, col = 'ff0000'x
  oplot, x, yctr+width/2.-[0,outer], thick = 8, col = 'ffa500'x
  oplot, x, yctr-width/2.+outer+[0,inter], thick = 8, col = '00aa00'x
  oplot, x, yctr+width/2.-outer-[0,inter], thick = 8, col = '00ff00'x
  oplot, x, yctr-width/2.+outer+inter+[0,inner], thick = 8, col = 255
  if ii eq 0 then begin
     cgtext, 6, 2+hlrObs/2./pscale, $
             '!18r!X!De!N', align = 0, charsize = 0.9, /data, orien = 90, $
             charthick = 3, col = 'ffffff'x
     cgtext, 22, 3, '5 kpc', align = 0, charsize = 0.9, /data, orien = 90, $
             charthick = 3, col = 'ffffff'x
  endif
  cgtext, !X.CRANGE[0]+2, !Y.CRANGE[1]-7, /data, $
          stats, charsize = 0.9, charthick = 3, col = 'ffffff'x
;  cgtext, !X.CRANGE[1]-1, !Y.CRANGE[0]+5, /data, $
;          fitsDir+'/'+objNo+'_'+pa, charsize = 0.8, charthick = 3, col = 'ffffff'x, align = 1
  plot, (findgen(s[0]) - s[0]/2) * pscale, (findgen(s[1]) -s[1]/2) * pscale, $
        /xsty, /ysty, xminor = 2, yminor = 2, /nodat, $
        pos = [!X.WINDOW[0], !Y.WINDOW[0], !X.WINDOW[1], !Y.WINDOW[1]], /noer, $
        xtickint = 1, ytickint = 1, $
        xthick = 4, ythick = 4, charsize = 1.25, charthick = 4, $
        xtitle = greek('delta')+'!18x!X ["]', $
        ytitle = greek('delta')+'!18y!X ["]'
  plot, (findgen(s[0]) - s[0]/2) * pscale, (findgen(s[1]) -s[1]/2) * pscale, $
        /xsty, /ysty, xminor = 2, yminor = 2, /nodat, $
        pos = [!X.WINDOW[0], !Y.WINDOW[0], !X.WINDOW[1], !Y.WINDOW[1]], /noer, $
        xtickint = 1, ytickint = 1, $
        xthick = 5, ythick = 5, col = 'ffffff'x, $
        xticklen = 0.03, yticklen = 0.03, $
        xtickname = replicate(' ', 60), $
        ytickname = replicate(' ', 60)
  
  st = !X.WINDOW[1]+0.1

  fullrun = xstop - st - 0.05

  ;;;;;;;;;;;;;;;;;
  ;; Extractions ;;
  ;;;;;;;;;;;;;;;;;

  plotsym, 0, /fill
  yttl = '!18f!X!D'+greek('lambda')+'!N [10!E-18!N erg s!E-1!N cm!E-2!N '+textoidl('\AA')+']'
  tp = optFlRes.PHOTOMETRY * 1d18
  tp = tp[where(tp gt 0)]
  yr = [0, (2*max(tp)) > (2 * max(optFlRes.SPEC, /nan))];0.3*min(tp) < min(optFlRes.ERR)        ; 1.1 * max(innFlres.FULL_MODEL * offs[0])]

  plot, lambda, optFlRes.SPEC, $
        xran = [2000,18000] / (1+zpick), yran = yr, xsty = 8+1, $ ;;1.1 * max(optFlRes.SPEC * norms[0])
        xtitle = greek('lambda')+'!Drest!N ['+texToIDL('\AA')+']', $
;        xtickname = replicate(' ', 60), $
        pos = [st,!Y.WINDOW[0], xstop, !Y.WINDOW[1]], $
        ytitle = yttl, /nodat, /noer, $
        xthick = 4, ythick = 4, charthick = 4, charsize = 1.25, yminor = 2, /ysty;, /ylog
  axis, xaxis = 1, xran = [2000.,18000.], xtitle = greek('lambda')+'!Dobs!N ['+textoidl('\AA')+']', $
        xthick = 4, ythick = 4, charthick = 4, charsize = 1.25, /xsty
  lines = getLinelist()
  for ii = 0, n_elements(lines.WAVES) - 1 do begin
     if optFlRes.SPEC[value_locate(lambda / (1+zpick), lines.WAVES[ii])] gt 0 then begin
        if NOT lines.ABSORPTION[ii] then $
           oplot, replicate(lines.WAVES[ii], 2), !Y.CRANGE[1] * [0.8,1], col = '555555'x, linesty = 1 $
        else begin
           oplot, replicate(lines.WAVES[ii], 2), !Y.CRANGE[1] * [0.8,1], col = '555555'x
        endelse
        cgtext, lines.WAVES[ii] - 15, !Y.CRANGE[1] - 0.2, align = 1, /data, $
                lines.NAMES[ii], charsize = 1, charthick = 2, col = '555555'x, orien = 90
     endif
  endfor

  tscol = '777777'x
  tmcol = 0
  results = optFlRes 
  spec  = results.SPEC
  model = results.FULL_MODEL
  err   = results.ERR
  phot  = results.PHOTOMETRY
  ephot = results.EPHOT
  mphot = results.MODELPHOT
  sn = mean(spec/err, /nan)
  lambda = results.LAMBDA
  oploterror, lambda / (1+zpick), spec, err, psym = 3, col = tscol, $
              errcol = tscol, errthick = 4, /nohat
  oplot, lambda / (1+zpick), smooth(spec, 12, /nan, /edge_truncate), col = '444444'x, thick = 3
;  oplot, lambda / (1+zpick), err, col = 200, thick = 2
  oploterror, results.PHOT_LAMBDA / (1+zpick), phot * 1d18, ephot * 1d18, $
              errthick = 2, psym = 8, symsize = 1.35, /nohat ;, $
;                    col = tscol, errcol = tscol
  oplot, results.PHOT_LAMBDA / (1+zpick), phot * 1d18, psym = 8, symsize = 0.8, col = tscol
  oplot, results.FULL_LAMBDA / (1+zpick), model, col = tmcol, thick = 4
  if keyword_set(resultsDir2) then begin
     model2 = results.OTHER_MODEL
     oplot, results.FULL_LAMBDA / (1+zpick), model2, col = tmcol, thick = 3, linesty = 2
  endif
  oplot, results.PHOT_LAMBDA / (1+zpick), mphot * 1d18, col = tmcol, psym = 4, symsize = 2, thick = 2
  cgtext, !X.CRANGE[1] * 0.99, model[value_locate(results.FULL_LAMBDA / (1+zpick), !X.CRANGE[1] * 0.95)] * 0.7, $
          /data, string(sn, f = '(F4.1)'), charsize = 1, charthick = 6, $
          align = 1
  cgtext, !X.CRANGE[1] * 0.99, model[value_locate(results.FULL_LAMBDA / (1+zpick), !X.CRANGE[1] * 0.95)] * 0.7, $
          /data, string(sn, f = '(F4.1)'), col = tscol, charsize = 1, charthick = 2, $
          align = 1
  if strpos(resultsDir, 'LogNormal') then begin
     sfhtag = 'LogNormal'
     sfhtag2 = 'Delayed Exp.'
  endif else begin
     sfhtag = 'Delayed Exp.'
     sfhtag2 = 'LogNormal'
  endelse
  keys = ['Opt. Extract', $
          sfhtag]
  kcols = ['777777'x, 0]
  stys = replicate(0, n_elements(keys))
  if keyword_set(resultsDir2) then begin
     keys = [keys, sfhtag2]
     kcols = [kcols, 0]
     stys = [stys, 2]
  endif
  legend, box = 0, /bottom, /right, $;/clear, $   ;pos = [!x.crange[0], 10.^!y.crange[0]*1.5], /data, $
          keys, col = kcols, $
          pspacing = 1, thick = 5, $
          charsize = 0.8, charthick = 3, $
          linesty = stys, spacing = 0, /horiz
;  mstel    = alog10(total(10.^tres[where(tres.REGION ne 'OPTFL')].LMASS[0])) - alog10(mag)

  mstel = tres[where(tres.REGION eq 'OPTFL')].LMASS[0] - alog10(mag) 
  mstelerr = sqrt(total(0.5 * (tres[where(tres.REGION ne 'OPTFL')].LMASS[2] - tres[where(tres.REGION ne 'OPTFL')].LMASS[1])^2)) 
;  mstelerr = 0.5 * (tres[where(tres.REGION eq 'OPTFL')].LMASS[2] -
;  tres[where(tres.REGION eq 'OPTFL')].LMASS[1]) << -- This is way too small...
  mstelerr = sqrt(mstelerr^2 + (1./alog(10) * emag / mag)^2)
  legend, pos = [!X.CRANGE[0]+100, !Y.CRANGE[1]-0.2], /data, box = 0, $;, /left, box = 0, $ ; pos = [0.17,8.4], /data
          [$                       ;'!18z!X!DGLASS!N='+string(bluDat.Z, f = '(F5.3)'), $
          '!18z!X='+string(zpick, f = '(F5.3)')+textoidl('\pm')+string(0.5*(zs[2,0] - zs[1,0]), f = '(F5.3)'), $
;           '!18z!X!Dfit!N='+string(optFlRes.ZFIT, f = '(F5.3)'), $
          '!18r!X!De!N='+string(hlrObs, f = '(F4.2)')+'"='+string(hlrKpc, f = '(F3.1)')+' kpc', $
;          'F140W='+string(bluDat.MAG, f = '(F5.2)'), $
          'log!18M!X!D*!N='+string(mstel, f = '(F4.1)')+textoidl('\pm')+string(mstelerr, f = '(F4.2)')], $
          charsize = 1, charthick = 4

  device, /close
  spawn, 'gv '+outputName+' &'
  set_plot, 'X'

;  stop
  
end

pro doGood

  sum = loadsummarydata('./')
  
  cd, 'MACS1149/'  
  plotOptfl, $
     '00900_1_pyspecfitLogNormalResults/', './', $
     cluster = 'MACS1149', $
     mag = sum[where(sum.ID eq '00900_1')].MU, $
     emag = sum[where(sum.ID eq '00900_1')].EMU, $
     resultsDir2 = '00900_1_pyspecfitPhotResults/', $
     /domasked, /color, tiff_name = '00900_1_rgb.tiff' ;, output = 'interOuterTest.eps'

  cd, '../MACS0744/'
  spawn, 'ln -s ../plotoptfl.pro .'
  plotOptfl, $
     '00660_2_pyspecfitLogNormalResults/', './', $
     cluster = 'MACS0744', $
     mag = sum[where(sum.ID eq '00660_2')].MU, $
     emag = sum[where(sum.ID eq '00660_2')].EMU, $
     resultsDir2 = '00660_2_pyspecfitPhotResults/', $
     /domasked, /color, tiff_name = '00660_2_rgb.tiff'

  cd, '../MACS2129/'
  spawn, 'ln -s ../plotoptfl.pro .'
  plotOptfl, $
     '00451_2_pyspecfitLogNormalResults/', './', $
     cluster = 'MACS2129', $
     mag = sum[where(sum.ID eq '00451_2')].MU, $
     emag = sum[where(sum.ID eq '00451_2')].EMU, $
     resultsDir2 = '00451_2_pyspecfitPhotResults/', $
     /domasked, /color, tiff_name = '00451_2_rgb.tiff';, output = 'interOuterTest.eps'

  cd, '../MACS1423/'
  spawn, 'ln -s ../plotoptfl.pro .'
  plotOptfl, $
     '01916_2_pyspecfitLogNormalResults/', './', $
     cluster = 'MACS1423', $
     mag = sum[where(sum.ID eq '01916_2')].MU, $
     emag = sum[where(sum.ID eq '01916_2')].EMU, $
     resultsDir2 = '01916_2_pyspecfitPhotResults/', $
     /domasked, /color, tiff_name = '01916_2_rgb.tiff';, output = 'interOuterTest.eps'
  cd, '..'
  
end
