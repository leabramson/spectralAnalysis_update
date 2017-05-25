function getGeoZoff, struct, region

  struct = mrdfits(struct, 1, /silent)
  
  ctr = struct.LSF_INNER_PARS[1]
  dl  = struct.LAMBDA[1] - strunct.LAMBDA[0]
  ml  = median(struct.LAMBDA)
  
  case REGION of
     'INTUP': c2 = struct.LSF_INTER_UP_PARS[1]
     'INTDN': c2 = struct.LSF_INTER_DN_PARS[1]
     'OUTUP': c2 = struct.LSF_OUTER_UP_PARS[1]
     'OUTDN': c2 = struct.LSF_OUTER_DN_PARS[1]
  endcase

  loff = (c2 - ctr) * dl
  zoff = dl / ml
  
  RETURN, zoff
end

;;
;;
;;

pro feedtopyspecfit, field, $
                     SOURCELIST = SOURCELIST, $
                     DOMASKED = domasked, $
                     ZMAX = zmax, $
                     ZLOCK = zlock ;; for second passes
  
  if KEYWORD_SET(DOMASKED) then $
     msksfx = '.masked' $
  else $
     msksfx = ''

  if NOT keyword_set(ZMAX) then zmax = 3.
  
  ;; Go to the data
  cd, field
 
  ;; Default to the "good source list"
  if NOT keyword_set(SOURCELIST) then begin
     check = file_info('pyspecWorthy_sourceInfo.list')
     if check.exists then $
        sourcelist = 'pyspecWorthy_sourceInfo.list' $
     else $
        sourcelist = 'fileInfo.list'
  end
  
  ;; Check to see the relevant programs are here
  check1 = file_info('fitspec.py')
  check2 = file_info('plot1res.pro')
  check3 = file_info('unifypyspecinputs.pro')
  check4 = file_info('plotfittingregion.pro')
  check5 = file_info('setzprior.pro')
  check6 = file_info('pyspecfitonetrace.pro')
  check7 = file_info('plotpyspecresults.pro')
  check8 = file_info('easyspecvis.pro')
  
  if NOT check1.EXISTS then spawn, 'ln -s ../fitspec.py .'
  if NOT check2.EXISTS then spawn, 'ln -s ../plot1res.pro .'
  if NOT check3.EXISTS then spawn, 'ln -s ../unifypyspecinputs.pro .'
  if NOT check4.EXISTS then spawn, 'ln -s ../plotfittingregion.pro .'
  if NOT check5.EXISTS then spawn, 'ln -s ../setzprior.pro .'
  if NOT check6.EXISTS then spawn, 'ln -s ../pyspecfitonetrace.pro .'
  if NOT check7.EXISTS then spawn, 'ln -s ../plotpyspecresults.pro .'
  if NOT check8.EXISTS then spawn, 'ln -s ../easyspecvis.pro'
  
  ;; Read in the info for the sources to fit
  print, ''
  print, 'GETTING SOURCE INFO FROM : '+sourcelist
  print, ''
  readcol, sourcelist, $
           FILENAME, $
           ID, $
           PA, $
           GRISM, $
           Z_PHOT, $
           Z_SPEC, $
           Z_SPEC_QUAL, $
           f = 'A,A,A,A,F,F,F', $
           comment = '#', /silent
  nfiles = n_elements(id)
  print, ''
  
  uid = id+pa  
  uid = uid[sort(uid)]
  uid = uid[UNIQ(uid)]
  ntorun = n_elements(uID)
  
  print, ' >>> WILL SEND '+string(nfiles, f = '(I4)')+' GRISM SPECTRA TO PYSPECFIT... '
  usrcs = id[sort(id)]
  usrcs = usrcs[UNIQ(usrcs)]
  print, ' >>> ('+string(n_elements(usrcs), f = '(I4)')+' unique objects)'

  counter = floor(findgen(nfiles) * 0.5)
  for ii = 0, ntorun - 1 do begin

     bdata = strmid(uID[ii], 0,5)+'_'+strmid(uID[ii], 5,1)+'_B.fits' 
     rdata = strmid(uID[ii], 0,5)+'_'+strmid(uID[ii], 5,1)+'_R.fits'
     
     easyspecvis, bdata, rdata
     
     if n_elements(zlock) eq 0 then begin
     
        ;; Find a "redshift"
        ;; weight them using the gigz quality flag
        tpz = z_phot[counter[ii]]
        tz  = z_spec[counter[ii]]
        tzq = z_spec_qual[counter[ii]]
        
        if tpz le 0 OR tpz gt zmax then begin
           zstart = tz
           zsran = 0.5
        endif else if tz le 0 then begin
           zstart = tpz
           zsran = 0.5
        endif else if tz gt 0 AND tpz gt 0 then begin
           ws = tzq / 4.
           wp = 1 - tzq / 4. 
           zstart = wp * tpz + ws * tz
           zsran  = (sqrt(0.5) * abs(tpz - tz)) < 0.5
        endif
        
        ;; Lock it into a reasonable redshift range
        zprior  = setzprior(1.5, 0.5);, /flat)
        zkernel = setzprior(zstart, zsran)
        zout    = zprior.P * zkernel.P
        ztot = total(zout, /cum) / total(zout)
        hits = value_locate(ztot, [0.16,0.84])
        foo  = max(zout, zpk)
        gout = [1, zprior.Z[zpk], mean(abs(zprior.Z[hits] - zprior.z[zpk]))]
;     gz      = gaussfit(zprior.Z, zout, nterms = 3, gout)
        
        ozstart = zstart
        ozsran  = zsran
        zstart  = gout[1]                ;; Max of the resultant PDF assuming the ZPRIOR prior
        zsran   = (2.35 * gout[2]) > 0.1 ;; FWHM of the resultant PDF or 0.1, w/e is bigger

          print, ''
          if n_elements(zlock) eq 0 then $
             print, ' >>> Beginning at z_guess = '+string(zstart, f = '(F4.2)'); $
;          else $
;             print, ' >>> REDSHIFT IS LOCKED TO '+string(zlock, f = '(F4.2)')
          
        window, 2, xsize = 400, ysize = 400
        plot, zprior.Z, zprior.P / total(zprior.P), $
              xran = [0,2 * zstart], $
              xtitle = '!18z!X', ytitle = '!18P!X(!18z!X)'
        oplot, zkernel.Z, zkernel.P / total(zkernel.p), thick = 3, col = 255
        oplot, zkernel.Z, zout / total(zout), thick = 2, col = '00ff00'x
        oplot, replicate(zstart, 2), !Y.CRANGE, col = '777777'x, thick = 1
        oplot, replicate(zstart - zsran/2, 2), !Y.CRANGE, $
               col = '777777'x, thick = 1, linesty = 1
        oplot, replicate(zstart + zsran/2, 2), !Y.CRANGE, $
               col = '777777'x, thick = 1, linesty = 1
        legend, /top, /left, $
                ['!18z!X!Din!N='+$
                 string(ozstart, f = '(F4.2)')+$
                 textoidl('\pm')+string(ozsran, f = '(F4.2)'), $
                 '!18z!X!Dout!N='+$
                 string(zstart, f = '(F4.2)')+$
                 textoidl('\pm')+string(zsran, f = '(F4.2)')], $
                pspacing = 0.1, col = [255, '00ff00'x], charsize = 1.1
        window, 0     
        
     endif else begin

        zstart = zlock
        zsran = 0

     endelse
     
     tid = strmid(uid[ii], 0, 5)
     tpa = strmid(uid[ii], 5, 1)
     
     regions = ['OPTFL'+msksfx, $  ;; Spec extractions
                'INNFL'+msksfx, $
                'INTUP'+msksfx, $
                'INTDN'+msksfx, $
                'OUTUP'+msksfx, $
                'OUTDN'+msksfx  ]
     
     for jj = 0, n_elements(regions) - 1 do begin

        region = regions[jj]
        
        print, ' >>> >>> RUNNING PYSPECFIT ON '+uid[ii]+'_'+pa[ii]+' '+region+' REGION !!! <<< '
        print, ''
        wait, 2.
        
        ;; Look at what you're doing here
        unifypyspecinputs, $
           tid+'_'+tpa+'_B_'+region+'.pyspec', $
           tid+'_'+tpa+'_R_'+region+'.pyspec', $
           tid+'_'+tpa+'_'+region+'.unified.dat'
        if jj eq 0 then $
           plotfittingregion, tid+'_'+tpa+'_'+region+'.unified.dat', zstart $
        else $
           plotfittingregion, tid+'_'+tpa+'_'+region+'.unified.dat', zfit

        ;; Do the fit
;        if (n_elements(zlock) eq 0 AND zsran gt 0.1) then begin
        if (n_elements(zlock) eq 0 AND zsran gt 0) then begin  ;; LEA changed 170329 -- CAN'T LOCK Z

;           stop
           
           if jj eq 0 then $
              pyspecfitonetrace, tid, tpa, region, zstart, zsran, $
                                 ncores = 10, outdir = region $
           else $
              pyspecfitonetrace, tid, tpa, region, zfit, 0.1, $
                                 ncores = 10, outdir = region           

        endif else if (n_elements(zlock) ne 0 OR zsran le 0.1) then begin
           print, ' *** *** *** *** *** *** *** *** *** *** *** *** *** '
           print, ' *** *** *** *** *** *** *** *** *** *** *** *** *** '
           print, ' >>> REDSHIFT IS LOCKED TO '+string(zstart, f = '(F4.2)')
           print, ' *** *** *** *** *** *** *** *** *** *** *** *** *** '
           print, ' *** *** *** *** *** *** *** *** *** *** *** *** *** '
           print, ''

;           stop
           
           pyspecfitonetrace, tid, tpa, region, zstart, zsran, $
                              ncores = 10, /zlock, outdir = region
        endif
        
        outname = tid+'_'+tpa+'_pyspecfitResults'
        check = file_info(outname)
        if NOT check.exists then $
           spawn, 'mkdir '+outname
        spawn, 'mv '+region+'* '+outname+'/'

        if jj eq 0 then begin
           spawn, 'cat '+outname+'/'+region+'.analyze | grep "z" > tz.txt'
           readcol, 'tz.txt', zfit, f = 'X,F', /silent
           dz = zfit - zstart
           print, ''
           print, '>>> OUTPUT Z DIFFERS FROM INITIAL GUESS BY : '+string(dz, f = '(F5.2)')+' <<<'
           print, ''
           print, ''
        endif

        print, ''
        print, '!!! '+region+' EXTRACTION FIT RESULTS !!!'
        print, ''
        spawn, 'cat '+outname+'/'+region+'.analyze'
        print, ''
        print, ''
        if jj eq 0 then $
           print, 'Input redshift : ', zstart $
        else $
           print, 'Input redshift : ', zfit

        plot1res, outname+'/'+region+'.bestspec0', outname+'/'+region+'.bestspec1', zfit
        wait, 2.0
        print, ''
        for kk = 0, jj do $
           spawn, 'cat '+ outname+'/'+regions[kk]+'.analyze | grep "z"'
        wait, 2.0
        
     endfor

     ;; Make a nice looking plot.
     if msksfx eq '' then $
        plotpyspecresults, outname, './'  $
     else $
        plotpyspecresults, outname, './',  /domasked     

     if ii lt ntorun - 1 then $
        print, ' >>> >>> ALL DONE FOR '+tid+'_'+tpa+' || MOVING TO NEXT SOURCE ... ' $
     else $
        print, ' >>> >>> >>> ALL DONE <<< <<< <<< '
     print, ''
     print, ''
     wait, 2.0
     
  endfor

  cd, '..'
  
end

;;
;;
;;

pro doGood

;  FEEDTOPYSPECFIT, 'MACS0717', sourcelist = '01266_tofit.list',  /domasked
;  FEEDTOPYSPECFIT, 'MACS0717', sourcelist = '01236_tofit.list', /domasked

;  FEEDTOPYSPECFIT, 'MACS1149', sourcelist = '01931_tofit.list', /domasked, zlock = 1.41 
;  FEEDTOPYSPECFIT, 'MACS1423', /domasked
  

;  FEEDTOPYSPECFIT, 'MACS2129', sourcelist = '01547_tofit.list', /domasked, $
;                   zlock = 1.158
  ;FEEDTOPYSPECFIT, 'MACS0744', sourcelist = '01266_2_tofit.list', /domasked

;  feedtopyspecfit, 'ABEL2744', /domasked
;  feedtopyspecfit, 'MACS0717', /domasked
  feedtopyspecfit, 'MACS0744', /domasked
;  feedtopyspecfit, 'MACS1149', /domasked
;  feedtopyspecfit, 'MACS1423', /domasked
  feedtopyspecfit, 'MACS2129', /domasked
;  feedtopyspecfit, 'RXJC1347', /domasked
;  feedtopyspecfit, 'RXJC2248', /domasked

                                ;  feedtopyspecfit, 'MACS2129',
                                ;  /domasked, sourcelist =
                                ;  '01547_tofit.list'
;  feedtopyspecfit, 'MACS1149', sourcelist = '00753_1_tofit.list'
  feedtopyspecfit, 'MACS1149', sourcelist = '00900_1_tofit.list'
  
end
