function getGeoZoff, struct, region

  struct = mrdfits(struct, 1, /silent)
  
  ctr = struct.LSF_INNER_PARS[1]
  dl  = struct.LAMBDA[1] - struct.LAMBDA[0]
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

pro feedtopyspecfitlognormal, field, $
                              SOURCELIST = SOURCELIST, $
                              DOMASKED = domasked
  
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
  if NOT (file_info('pyspecfitonetracelgnml.pro')).exists  then $
     spawn, 'ln -s ../pyspecfitonetracelgnml.pro .'
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

     if NOT (file_info(tid+'_'+tpa+'_OUTER.sed')).EXISTS then begin
        foldspecphot, tid+'_'+tpa, 'INTER' ;; added LEA 2017 06 09                                           
        foldspecphot, tid+'_'+tpa, 'OUTER' ;; Produces the spec/phot files (masked only) for INTER OUTER fits
     endif
     
     ;; Setup the output directory
     outname = tid+'_'+tpa+'_pyspecfitLogNormalResults'
     check = file_info(outname)
     if NOT check.exists then $
        spawn, 'mkdir '+outname

     for jj = 0, n_elements(regions) - 1 do begin

        bdata = strmid(uID[ii], 0,5)+'_'+strmid(uID[ii], 5,1)+'_B.fits' 
        rdata = strmid(uID[ii], 0,5)+'_'+strmid(uID[ii], 5,1)+'_R.fits'
        
        region = regions[jj]
        
        print, f = '(%"\n >>> >>> RUNNING PYSPECFIT ON %s_%i %s REGION !!! <<<\n ")', $
               strmid(uID[ii], 0,5), pa[ii], region

        if jj eq 0 then begin
           spawn, 'cat '+tid+'_'+tpa+'_pyspecfitPhotResults/'+regions[jj]+'.analyze | grep z > tmpZ.txt'
           readcol, 'tmpZ.txt', zfit, f = 'X,F'
        endif
           
        print, ' *** *** *** *** *** *** *** *** *** *** *** *** *** '
        print, ' *** *** *** *** *** *** *** *** *** *** *** *** *** '
        print, ' >>> REDSHIFT IS LOCKED TO '+string(zfit, f = '(F5.3)')
        print, ' *** *** *** *** *** *** *** *** *** *** *** *** *** '
        print, ' *** *** *** *** *** *** *** *** *** *** *** *** *** '

        wait, 1.
        
        ;; Look at what you're doing here
        unifypyspecinputs, $
           tid+'_'+tpa+'_B_'+region+'.pyspec', $
           tid+'_'+tpa+'_R_'+region+'.pyspec', $
           tid+'_'+tpa+'_'+region+'.unified.dat'
        plotfittingregion, tid+'_'+tpa+'_'+region+'.unified.dat', zfit

        ;; Grab the photometry
        pfile = tid+'_'+tpa+'_'+basic_regions[jj]+'.sed'
        
        ;; Do the fit. For the lognormal, the redshift must be
        ;; determined a priori. Hence, lock them all to the OPTIMAL
        ;; solution from the DelayedExponential fit.        
        pyspecfitonetracelgnml, tid, tpa, region, zfit, $
                                photfile = pfile, $
                                ncores = 12, outdir = outname+'/'+region

        print, f = '(%"\n!!! %s EXTRACTION FIT RESULTS !!!\n")', region
        spawn, 'cat '+outname+'/'+region+'.analyze'
        print, ''
        print, ''

        ;; Plot some stuff
        plot1resps, outname+'/'+region+'.bestspec0', outname+'/'+region+'.bestspec1', $
                    redshift = zfit, photfile = outname+'/'+region+'.bestphot', $
                    output = outname+'/'+region+'_fit.eps'
        spawn, 'gv '+outname+'/'+region+'_fit.eps &'

        plot1res, outname+'/'+region+'.bestspec0', outname+'/'+region+'.bestspec1', redshift = zfit
        wait, 2.0
        print, ''
        for kk = 0, jj do $
           spawn, 'cat '+ outname+'/'+regions[kk]+'.analyze | grep tau'
        wait, 1.0
        
     endfor

     if ii lt ntorun - 1 then $
        print, f = '(%"\n >>> >>> ALL DONE FOR %s_%i || MOVING TO NEXT SOURCE ... \n")', $
               tid, tpa $
     else $
        print, f = '(%"\n\n >>> >>> >>> ALL DONE <<< <<< <<< \n\n")'

     wait, 1.0
     
  endfor

  cd, '..'
  
end

pro doGood;

;  feedtopyspecfitlognormal, 'ABEL0370', /domasked, sourcelist = 'forLognormalFitting.list'
  feedtopyspecfitlognormal, 'MACS0717', /domasked, sourcelist = 'forLognormalFitting.list'
  feedtopyspecfitlognormal, 'MACS0744', /domasked, sourcelist = 'forLognormalFitting.list'
  feedtopyspecfitlognormal, 'MACS1149', /domasked, sourcelist = 'forLognormalFitting.list'
  feedtopyspecfitlognormal, 'MACS1423', /domasked, sourcelist = 'forLognormalFitting.list'
  feedtopyspecfitlognormal, 'MACS2129', /domasked, sourcelist = 'forLognormalFitting.list'
 
end
