;function getGeoZoff, struct, region
;
;  struct = mrdfits(struct, 1, /silent)
;  
;  ctr = struct.LSF_INNER_PARS[1]
;  dl  = struct.LAMBDA[1] - struct.LAMBDA[0]
;  ml  = median(struct.LAMBDA)
;  
;  case REGION of
;     'INTUP': c2 = struct.LSF_INTER_UP_PARS[1]
;     'INTDN': c2 = struct.LSF_INTER_DN_PARS[1]
;     'OUTUP': c2 = struct.LSF_OUTER_UP_PARS[1]
;     'OUTDN': c2 = struct.LSF_OUTER_DN_PARS[1]
;  endcase
;
;  loff = (c2 - ctr) * dl
;  zoff = dl / ml
;  
;  RETURN, zoff
;end

;;
;;
;;

pro feedtopyspecfit2, field, $
                      SOURCELIST = SOURCELIST, $
                      DOMASKED = domasked, $
                      ZMAX = zmax, $
                      ZLOCK = zlock, $
                      SUFFIX = suffix;; for second passes

  if NOT keyword_set(SUFFIX) then $
     suffix = '.masked'
  
  if KEYWORD_SET(DOMASKED) then $
     msksfx = suffix $
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
  if NOT (file_info('fitspec*.py')).exists            then $
     spawn, 'ln -s ../fitspec*.py .'
  if NOT (file_info('plot1res.pro')).exists           then $
     spawn, 'ln -s ../plot1res.pro .'
  if NOT (file_info('unifypyspecinputs.pro')).exists  then $
     spawn, 'ln -s ../unifypyspecinputs.pro .'
  if NOT (file_info('plotfittingregion.pro')).exists  then $
     spawn, 'ln -s ../plotfittingregion.pro .'
  if NOT (file_info('setzprior.pro')).exists          then $
     spawn, 'ln -s ../setzprior.pro .'
  if NOT (file_info('pyspecfitonetrace.pro')).exists  then $
     spawn, 'ln -s ../pyspecfitonetrace.pro .'
  if NOT (file_info('plotpyspecresults.pro')).exists  then $
     spawn, 'ln -s ../plotpyspecresults.pro .'
  if NOT (file_info('easyspecvis.pro')).exists        then $
     spawn, 'ln -s ../easyspecvis.pro'
  if NOT (file_info('getsrcimage.pro')).exists        then $
     spawn, 'ln -s ../getsrcimage.pro'
  if NOT (file_info('matchphot.pro')).exists          then $
     spawn, 'ln -s ../matchphot.pro'
  if NOT (file_info('printphotforpyspec.pro')).exists then $
     spawn, 'ln -s ../printphotforpyspec.pro'
  if NOT (file_info('plot1resps.pro')).exists then $
     spawn, 'ln -s ../plot1resps.pro'
  if NOT (file_info('maketifffromstack.pro')).exists then $
     spawn, 'ln -s ../maketifffromstack.pro'
  if NOT (file_info('foldspecphot.pro')).exists then $
     spawn, 'ln -s ../foldspecphot.pro'
  if NOT (file_info('getgeozoff.pro')).exists then $
     spawn, 'ln -s ../getgeozoff.pro'
  spawn, 'ln -s ../*.sex .'
  spawn, 'ln -s ../*.param .'
  spawn, 'ln -s ../*.conf .'
  
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

     ;; Grab the resolved photometry
     pname = repstr(bdata, '.fits', '_allbandImages.fits')     
     if NOT (file_info(pname)).EXISTS then $
        getsrcimage, field, master = bdata, /ds9
     maketifffromstack, pname, repstr(pname, 'B_allbandImages.fits', 'rgb.tiff')
     matchphot, master = bdata, multiband = pname, $
                dustfile = 'GALACTIC_EXTINCTION.dat'
     
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
     zprior  = setzprior(1.5, 0.5) ;, /flat)
     zkernel = setzprior(zstart, zsran)
     zout    = zprior.P * zkernel.P
     ztot = total(zout, /cum) / total(zout)
     hits = value_locate(ztot, [0.16,0.84])
     foo  = max(zout, zpk)
     gout = [1, zprior.Z[zpk], mean(abs(zprior.Z[hits] - zprior.z[zpk]))]
;     gz      = gaussfit(zprior.Z, zout, nterms = 3, gout)
     
     ozstart = zstart
     ozsran  = zsran
     zstart  = gout[1]                   ;; Max of the resultant PDF assuming the ZPRIOR prior
     zsran   = (2.35 * gout[2]) > 0.1    ;; FWHM of the resultant PDF or 0.1, w/e is bigger
     
     print, ''
     print, ' >>> Beginning at z_guess = '+string(zstart, f = '(F4.2)') ; $
          
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
     
     tid = strmid(uid[ii], 0, 5)
     tpa = strmid(uid[ii], 5, 1)

     basic_regions = ['OPTFL', $
                      'INNFL', $
                      'INTER', $ ;; Added LEA 2017 06 09         
                      'OUTER', $ ;; Mean of INTUP/DN and OUTUP/DN
                      'INTUP', $
                      'INTDN', $
                      'OUTUP', $
                      'OUTDN']   
     
     regions = basic_regions+msksfx  ;; Spec extractions
     
     if NOT (file_info(tid+'_'+tpa+'_OUTER.pyspec')).EXISTS then begin  
        foldspecphot, tid+'_'+tpa, 'INTER' ;; added LEA 2017 06 09                                           
        foldspecphot, tid+'_'+tpa, 'OUTER' ;; Produces the spec/phot files (masked only) for INTER OUTER fits
     endif
     
     ;; Setup the output directory
     outname = tid+'_'+tpa+'_pyspecfitPhotResults'
     check = file_info(outname)
     if NOT check.exists then $
        spawn, 'mkdir '+outname
     
     for jj = 0, n_elements(regions) - 1 do begin
        
        bdata = strmid(uID[ii], 0,5)+'_'+strmid(uID[ii], 5,1)+'_B.fits' 
        rdata = strmid(uID[ii], 0,5)+'_'+strmid(uID[ii], 5,1)+'_R.fits'
        
        region = regions[jj]
        
        print, f = '(%"\n >>> >>> RUNNING PYSPECFIT ON %s_%i %s REGION !!! <<<\n ")', $
               strmid(uID[ii], 0,5), pa[ii], region
        wait, 1.
        
        ;; Look at what you're doing here
        unifypyspecinputs, $
           tid+'_'+tpa+'_B_'+region+'.pyspec', $
           tid+'_'+tpa+'_R_'+region+'.pyspec', $
           tid+'_'+tpa+'_'+region+'.unified.dat'
        if jj eq 0 then $
           plotfittingregion, tid+'_'+tpa+'_'+region+'.unified.dat', zstart $
        else $
           plotfittingregion, tid+'_'+tpa+'_'+region+'.unified.dat', zfit
        
        ;; Print the photometry
        pfile = tid+'_'+tpa+'_'+basic_regions[jj]+'.sed'
        if region ne 'INTER' and region ne 'OUTER' then $
           printphotforpyspec, repstr(bdata, '.fits', '_resolvedSED.fits'), $
                               region = basic_regions[jj], $
                               output = pfile
        
        ;; Do the fit. Find a baseline solution using the optimal
        ;; extraction and optimal profile weighted photometry. Then,
        ;; use those results as inputs for the other regions, fixing
        ;; the redshift to the INNER trace if necessary using the LSF
        ;; centroids to establish wavelength zeropoint/redshift
        ;; offsets to limit the zrange examined by the fitter.
        if jj eq 0 then begin
;        if jj le 3 then begin
           
           pyspecfitonetrace, tid, tpa, region, zstart, zsran, $
                              photfile = pfile, $
                              ncores = 12, outdir = outname+'/'+region
           spawn, 'cat '+outname+'/'+region+'.analyze | grep z > tz.txt'
           readcol, 'tz.txt', zfit, zfl, zfh, f = 'X,F,F,F', /silent
           dz = zfit - zstart
           zerr = 0.5 * (zfh - zfl)

           print, f = '(%"\n\n>>> OUTPUT Z DIFFERS FROM INITIAL GUESS BY : %f <<<\n\n")', dz

           zfit2 = zfit
           
        endif else begin
           if basic_regions[jj] eq 'INNFL' then begin

              zfit2 = zfit
              pyspecfitonetrace, tid, tpa, region, zfit2, zerr, $
                                 photfile = pfile, $
                                 ncores = 12, outdir = outname+'/'+region

           endif else begin

              if basic_regions[jj] ne 'INTER' $
                 AND $
                 basic_regions[jj] ne 'OUTER' $
              then begin
                 zoffb = getGeoZoff(bdata, basic_regions[jj])
                 zoffr = getGeoZoff(rdata, basic_regions[jj])
              endif else if basic_regions[jj] eq 'INTER' then begin
                 zoffb1 = getGeoZoff(bdata, 'INTUP')
                 zoffr1 = getGeoZoff(rdata, 'INTUP')
                 zoffb2 = getGeoZoff(bdata, 'INTDN')
                 zoffr2 = getGeoZoff(rdata, 'INTDN')
                 zoffb  = 0.5 * (zoffb1 + zoffb2)
                 zoffr  = 0.5 * (zoffr1 + zoffr2)
              endif else if basic_regions[jj] eq 'OUTER' then begin
                 zoffb1 = getGeoZoff(bdata, 'OUTUP')
                 zoffr1 = getGeoZoff(rdata, 'OUTUP')
                 zoffb2 = getGeoZoff(bdata, 'OUTDN')
                 zoffr2 = getGeoZoff(rdata, 'OUTDN')
                 zoffb  = 0.5 * (zoffb1 + zoffb2)
                 zoffr  = 0.5 * (zoffr1 + zoffr2)
              endif
              zoff = 0.5 * (zoffb + zoffr)
              
              if zoff lt 2*zerr then begin
                 zerr2 = 0
                 zoff  = 0
                 print, f = '(%"\n\n >>> >>> Extaction centroid in-line w/ central zone <<< <<<\n'+$
                        ' >>> >>> forcing z = z_ctr \n\n")'
                 print, ' *** *** *** *** *** *** *** *** *** *** *** *** *** '
                 print, ' *** *** *** *** *** *** *** *** *** *** *** *** *** '
                 print, ' >>> REDSHIFT IS LOCKED TO '+string(zfit2, f = '(F5.3)')
                 print, ' *** *** *** *** *** *** *** *** *** *** *** *** *** '
                 print, ' *** *** *** *** *** *** *** *** *** *** *** *** *** '
              endif else $
                 zerr2 = zerr
                            
              pyspecfitonetrace, tid, tpa, region, (zfit2 + zoff), zerr2, $
                                 photfile = pfile, $
                                 ncores = 12, $
                                 outdir = outname+'/'+region
              
           endelse

           ;; Get some stats
           if zerr ne 0 then begin
              spawn, 'cat '+outname+'/'+region+'.analyze | grep z > tz.txt'
              readcol, 'tz.txt', zfit2, zfl2, zfh2, f = 'X,F,F,F', /silent
           endif
           
        endelse
        
;        if (n_elements(zlock) ne 0 OR zsran le 0.1) then begin
;           print, ' *** *** *** *** *** *** *** *** *** *** *** *** *** '
;           print, ' *** *** *** *** *** *** *** *** *** *** *** *** *** '
;           print, ' >>> REDSHIFT IS LOCKED TO '+string(zstart, f = '(F5.3)')
;           print, ' *** *** *** *** *** *** *** *** *** *** *** *** *** '
;           print, ' *** *** *** *** *** *** *** *** *** *** *** *** *** '
;           print, ''           
;           pyspecfitonetrace, tid, tpa, region, zstart, 0, $
;                              ncores = 10, /zlock, outdir = region
;        endif

        print, f = '(%"\n!!! %s EXTRACTION FIT RESULTS !!!\n")', region
        spawn, 'cat '+outname+'/'+region+'.analyze'
        print, ''
        print, ''
        if jj eq 0 then $
           print, 'Input redshift : ', zstart $
        else $
           print, 'Input redshift : ', zfit2
        
        ;; Plot some stuff
        plot1resps, outname+'/'+region+'.bestspec0', outname+'/'+region+'.bestspec1', $
                    redshift = zfit2, photfile = outname+'/'+region+'.bestphot', $
                    output = outname+'/'+region+'_fit.eps'
        spawn, 'gv '+outname+'/'+region+'_fit.eps &'
        
        plot1res, outname+'/'+region+'.bestspec0', outname+'/'+region+'.bestspec1', redshift = zfit2
        wait, 2.0
        print, ''
        for kk = 0, jj do $
           spawn, 'cat '+ outname+'/'+regions[kk]+'.analyze | grep z'
        wait, 2.0
        
     endfor

     ;; Make a nice looking plot.
;     if msksfx eq '' then $
;        plotpyspecresults, outname, './'  $
;     else $
;        plotpyspecresults, outname, './',  /domasked     

     if ii lt ntorun - 1 then $
        print, f = '(%"\n >>> >>> ALL DONE FOR %s_%i || MOVING TO NEXT SOURCE ... \n")', $
               tid, tpa $
     else $
        print, f = '(%"\n\n >>> >>> >>> ALL DONE <<< <<< <<< \n\n")'

     wait, 2.0
     
  endfor

  cd, '..'
  
end

pro doGood

;  feedtopyspecfit2, 'ABEL0370', /domasked
;  feedtopyspecfit2, 'ABEL2744', /domasked
  feedtopyspecfit2, 'MACS0717', /domasked
  feedtopyspecfit2, 'MACS0744', /domasked
  feedtopyspecfit2, 'MACS0744', /domasked, sourcelist = '00660_tofit.list'
  feedtopyspecfit2, 'MACS1149', /domasked
  feedtopyspecfit2, 'MACS2129', /domasked
  feedtopyspecfit2, 'MACS2129', /domasked, sourcelist = '00451_2_tofit.list'
  feedtopyspecfit2, 'MACS1423', /domasked

;  feedtopyspecfit2, 'RXJC1347', /domasked
;  feedtopyspecfit2, 'RXJC2248', /domasked

  feedtopyspecfitlognormal, 'MACS0744', /domasked, sourcelist = 'forLognormalFitting.list'
  feedtopyspecfitlognormal, 'MACS1149', /domasked, sourcelist = 'forLognormalFitting.list'
  feedtopyspecfitlognormal, 'MACS2129', /domasked, sourcelist = 'forLognormalFitting.list'
  feedtopyspecfitlognormal, 'MACS2129', /domasked, sourcelist = '00451_2_tofit.list'
  feedtopyspecfitlognormal, 'MACS0717', /domasked, sourcelist = 'forLognormalFitting.list'
  feedtopyspecfitlognormal, 'MACS1423', /domasked, sourcelist = 'forLognormalFitting.list'

  
end

pro check01266

  feedtopyspecfit2, 'MACS0717', /domasked, suffix = '.masked2', $
                    sourcelist = '01266_1_toRefit.list'
  feedtopyspecfitlognormal, 'MACS0717', /domasked, suffix = '.masked2', $
                            sourcelist = '01266_1_toRefit.list'
  
  
end

;;
;;
;;

;pro doGood;

;  FEEDTOPYSPECFIT, 'MACS0717', sourcelist = '01266_tofit.list',  /domasked

;  FEEDTOPYSPECFIT, 'MACS1149', sourcelist = '01931_tofit.list', /domasked, zlock = 1.41 
;  FEEDTOPYSPECFIT, 'MACS1423', /domasked
  

;  FEEDTOPYSPECFIT, 'MACS2129', sourcelist = '01547_tofit.list', /domasked, $
;                   zlock = 1.158
  ;FEEDTOPYSPECFIT, 'MACS0744', sourcelist = '01266_2_tofit.list', /domasked

;  feedtopyspecfit, 'ABEL2744', /domasked
;  feedtopyspecfit, 'MACS0717', /domasked
;  feedtopyspecfit2, 'MACS0744', /domasked
;  feedtopyspecfit2, 'MACS0717', sourcelist = '01236_tofit.list', /domasked
;  feedtopyspecfit2, 'MACS1149', /domasked
;  feedtopyspecfit2, 'MACS1423', /domasked
;  feedtopyspecfit2, 'MACS2129', /domasked
;  feedtopyspecfit, 'RXJC1347', /domasked
;  feedtopyspecfit, 'RXJC2248', /domasked

                                ;  feedtopyspecfit, 'MACS2129',
                                ;  /domasked, sourcelist =
                                ;  '01547_tofit.list'
;  feedtopyspecfit, 'MACS1149', sourcelist = '00753_1_tofit.list'
;  feedtopyspecfit2, 'MACS1149', sourcelist = '00900_1_tofit.list'
  
;end
