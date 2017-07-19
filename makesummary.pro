;; Combine all the relevant fitting-derived info for a source

pro makesummary, expDir, lgnDir, $
                 FITSFILE = fitsfile, $
                 SOURCEID = sourceID, $
                 OUTDIR = outdir

  check = file_info(outdir)
  if not check.exists then $
     spawn, 'mkdir '+outdir
  
  regions = ['INNFL', 'INTER', 'OUTER', $
             'INTUP', 'INTDN', 'OUTUP', 'OUTDN', $
             'OPTFL']
 
  ;; Read the DelayedExponential Results
  analyzetofits, expdir, $
                 outdir+'/'+sourceID+'_exp_bestfit.fits', $
                 FITSFILE = fitsfile
  for ii = 0, n_elements(regions) - 1 do $
     reconstructhist, expDir, regions[ii], $
                      output = outdir+'/'+sourceID+'_'+regions[ii]+'_exp_recons.fits'

  combinereconstructions, outdir+'/'+sourceID+'_XXXXX_exp_recons.fits', $
                             output = outdir+'/'+sourceID+'_exp_recons.fits'

;  for ii = 0, n_elements(regions) - 1 do begin
;     regionEvo = mrdfits(outdir+'/'+sourceID+'_'+regions[ii]+'_exp_recons.fits', 1)
;     if ii eq 0 then $
;        regionsEvo = regionEvo $
;     else $
;        regionsEvo = struct_addtags(regionsEvo, regionEvo)
;  endfor
;  mwrfits, regionsEvo, outdir+'/'+sourceID+'_exp_recons.fits', /create
  
  ;; Read the LogNormal Results
  check = file_info(lgnDir+'/OUTDN.masked.analyze')
  if check.exists then begin
     analyzetofits, lgndir, outdir+'/'+sourceID+'_lgn_bestfit.fits', $
                    sfh = 'LogNormal', $
                    FITSFILE = fitsfile
     for ii = 0, n_elements(regions) - 1 do $
        reconstructhist, lgnDir, regions[ii], $
                         output = outdir+'/'+sourceID+'_'+regions[ii]+'_lgn_recons.fits', $
                         sfh = 'LogNormal'

     combinereconstructions, outdir+'/'+sourceID+'_XXXXX_lgn_recons.fits', $
                             output = outdir+'/'+sourceID+'_lgn_recons.fits'

;     master = mrdfits(outdir+'/'+sourceID+'_lgn_bestfit.fits', 1)
;     for ii = 0, n_elements(regions) - 1 do begin
;        regionEvo = mrdfits(outdir+'/'+sourceID+'_'+regions[ii]+'_lgn_recons.fits', 1)
;        if ii eq 0 then $
;           regionsEvo = regionEvo $
;        else $
;           regionsEvo = struct_addtags(regionsEvo, regionEvo)
;     endfor
;     mwrfits, regionsEvo, outdir+'/'+sourceID+'_lgn_recons.fits', /create
  endif
  
end

pro makeAllsummaries

  readcol, '20170530_finalSampleForAnalysis.list', files, $
           f = 'A', comment = '#'
  nfiles = n_elements(files)

  for ii = 0, nfiles - 1 do begin
     expDir = files[ii]+'_pyspecfitPhotResults/'
     lgnDir = files[ii]+'_pyspecfitLogNormalResults/'
     sourceID = strmid(files[ii], 9, strlen(files[ii]))
     makesummary, expDir, lgnDir, $
                  SOURCEID = sourceID, OUTDIR = 'resultsSummaries', $
                  FITSFILE = '/home/labramson/Projects/GLASS/spectralAnalysis/'+$
                  files[ii]+'_B.fits'

;     mwrfits, rads, 'resultsSummaries/'+sourceID+'_radii.fits', /create
  endfor
  
end

;pro add660
;
;  expdir = 'MACS0744/00660_2_pyspecfitPhotResults/'
;  lgndir = 'MACS0744/00660_2_pyspecfitLogNormalResults/'
;  sourceID = '00660_2'
;  makesummary, expDir, lgnDir, $
;               SOURCEID = sourceID, OUTDIR = 'resultsSummaries'
;  
;end

