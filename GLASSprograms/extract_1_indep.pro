;; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
;; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
;; Perform an extraction for 1 object allowing for geometric effects
;; by not stacking the UPPER and LOWER radial extractions.
;; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
;; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

function extract_1_indep, infile, z, $
                          Z_QUAL  = z_qual, $
                          Z_PHOT  = z_phot, $
                          PA      = pa, $
                          DOMASK  = domask, $
                          REFIT   = refit, $
                          INSPECT = inspect, $
                          FCFRACT = fcfract, $
                          ACS     = acs

  if N_ELEMENTS(z) eq 0       then z = -99
  if NOT keyword_set(Z_QUAL)  then z_qual = -99
  if NOT keyword_set(Z_PHOT)  then z_phot = -99
  if NOT keyword_SET(PA)      then pa = -99
  if NOT keyword_set(DOMASK)  then domask = 0 else domask = 1 ;; mask rather than subtract contamination
  if NOT keyword_set(REFIT)   then refit = 0 else refit = 1    ;; refit continuum centerline
  if NOT keyword_set(INSPECT) then inspect = 0 else inspect = 1
  if NOT keyword_set(FCFRACT) then fcFract = 0.01
  if NOT keyword_set(ACS)     then acs = 0 else acs = 1
  
  ;; READ CRAP
  null    = mrdfits(infile, 0, thead, /silent)
  spec    = mrdfits(infile, 'sci', head, /silent)
  model   = mrdfits(infile, 'model', /silent)
  contam  = mrdfits(infile, 'contam', /silent)
  sens    = mrdfits(infile, 'sens', /silent)
  trace   = mrdfits(infile, 'ytrace', /silent)
  err     = mrdfits(infile, 'wht', /silent)
  di      = mrdfits(infile, 'dsci', diHead, /silent)

  zmask   = spec
  zmask[*,*] = 0
  zmask[where(spec eq 0)] = 1
  fcmap   = contam / model ;; added 5/9/17
  
  badflag = 0
  if n_elements(spec[*,0]) ne n_elements(trace) then begin
     print, ''
     print, '!!! INCOMPLETE SPECTRUM !!!'
     trace = trace[0:n_elements(spec[*,0])-1]
     badflag = 1
  endif
  
  ;; BASIC INFO
  exptime = sxpar(thead, 'EXPTIME')   ;; SECONDS
  mag     = sxpar(thead, 'MAG')       ;; Probably F140W or IRTOT?

  
  ;; TO POISSON UNITS  
  spec   *= exptime
  contam *= exptime
  err[where(err le 0)] = median(err) ;; to avoid NaNs..
  var     = (err * exptime)^2

  
  ;; MAKE A WAVELENGTH AXIS
  l0      = sxpar(head, 'CRVAL1')
  dl      = sxpar(head, 'CD1_1')
  run     = sxpar(head, 'NAXIS1')
  lambda  = findgen(run) * dl + l0

  
  ;; MAKE A SPATIAL AXIS
  yctr    = sxpar(head, 'CRPIX2') - 1
  tx      = findgen(sxpar(head, 'NAXIS2'))
  height  = n_elements(tx)

  
  ;; MAKE A MASK IMAGE
  mask       = contam
  mask[where(contam ge 2 * sqrt(var), compl = good)] = 1
  mask[good] = 0                                      ;; Make the mask
  mask       = smooth(mask, [15,5], /edge_wrap, /nan) ;; Grow the mask
  mask[where(mask ge 1d-4 AND fcmap ge fcFract, compl = good)] = 1 ;; make sure you're masking things w/ significant fluxfractions, else ignore
  mask[good] = 0
  mask = mask > zmask ;; LEA 2017 26 May

  
  ;; MAKE A SCIENCE IMAGE
  sci   = (spec - contam)


  ;; GET SPATIAL PROFILE
  prof    = getProfile(di, acs = acs)
;  profile = prof.LSF       
  re      = ceil(prof.RE)        ;; Returns a pseudo "half-light radius"
  params  = prof.SPATIAL_PARAMS

  
  ;; FIND A NEW CTRLINE IF YOU MUST
  if NOT refit then begin
     ctroid = fltarr(2,run)
     ctroid[0,*] = trace
     ctroid[1,*] = re
  endif else begin
     ctroid = fltarr(4,run)
     for jj = 0, run - 1 do begin
        if 20 * mean(reform(sci[jj,*]) / sqrt(reform(var[jj,*]))) gt 10 then begin  ;;  "S/N in central region"
           g = gaussfit(tx, reform(sci[jj,*]), $
                        nterms = 3, out, $
                        sig = sig)
           ctroid[*,jj] = [out[1:2], sig[1:2]]
        endif
     endfor
     use = where(abs(ctroid[0,*] - height/2) le 10, nuse)
     tofit = use[where(abs(ctroid[0,use] - mean(ctroid[0,use])) $
                       le 2 * stddev(ctroid[0,use]))]
     cline = poly_fit(lambda[tofit], ctroid[0,tofit], $
                      measure_err = ctroid[2, tofit], 1)
     trace = cline[0] + lambda * cline[1]
  endelse

    
  ;; MAKE A MODEL IMAGE FOR OPTIMAL EXTRACTION
  modIm = sci  ;; Will contain the profile for optimal extraction
  for jj = 0, run - 1 do begin
     u = (tx - trace[jj]) / params[2]
     modIm[jj,*] = 1. / (u^2 + 1)^params[3]
     modIm[jj,*] /= total(modIm[jj,*])
  endfor
  extrIm = sci

  opti         = dblarr(run) ;; Optimal extractions will go here
  opti_Var     = dblarr(run)
  optiMask     = dblarr(run)
  optiMask_Var = dblarr(run)

  
  ;;  DEFINE EXTRACTION RADIAL REGIMES
  ;;  Do it by RE if RE is large enough, else just take rows of pix.
  if 0.25 * re ge 1 then begin
     innerup_bound = ceil ( trace   + 0.25 * re )                      ;;        R/Re < 0.25 
     innerdn_bound = floor( trace   - 0.25 * re )  
     interup_bound = ceil ( innerup_bound + 0.50 * re )                ;; 0.25 < R/Re < 0.75
     interdn_bound = floor( innerdn_bound - 0.50 * re )
     outerup_bound = ceil ( interup_bound + 0.75 * re ) < (height - 1) ;; 0.75 < R/Re < 1.50
     outerdn_bound = floor( interdn_bound - 0.75 * re ) > 0
  endif else begin
     innerup_bound = ceil ( trace         + 1 )
     innerdn_bound = floor( trace         - 1 )  
     interup_bound = ceil ( innerup_bound + 1 )                
     interdn_bound = floor( innerdn_bound - 1 )
     outerup_bound = ceil ( interup_bound + 1 ) < (height - 1) 
     outerdn_bound = floor( interdn_bound - 1 ) > 0
  endelse

  if max(outerup_bound) eq height - 1 $
     OR $
     min(outerdn_bound) eq 0 $
  then $
     badflag = 1
  
  ;; GET AN OPTIMALLY WIDE "OUTER APERTURE"
  optiWide = getoptizone(sci, var, $
                         clineUp = innerUp_bound, $
                         clineDn = innerDn_bound)
  optiWideUp_bound = innerUp_bound + optiWide
  optiWideDn_bound = innerDn_bound - optiWide

  
  ;; DEFINE THE REGIONS FOR THE 5 INDEPENDENT RADIAL EXTRACTIONS
  skies        = dblarr(run)
  integ        = dblarr(run)  ;;  total flux (w/in OUTER limits)
  outerUP      = dblarr(run)  ;;  upper OUTER region
  outerDN      = dblarr(run)  ;;  lower OUTER region
  interUP      = dblarr(run)
  interDN      = dblarr(run)
  inner        = dblarr(run)
  optiOutUP    = dblarr(run)
  optiOutDN    = dblarr(run)
 
  ivar         = dblarr(run)
  outvarUP     = dblarr(run)
  outvarDN     = dblarr(run)
  intvarUP     = dblarr(run)
  intvarDN     = dblarr(run)
  innvar       = dblarr(run)
  optiOutVarUP = dblarr(run)
  optiOutVarDN = dblarr(run)
  
  nAlls        = intarr(run)
  nInns        = intarr(run)
  nIntsUP      = intarr(run)
  nIntsDN      = intarr(run)
  nOutsUP      = intarr(run)
  nOutsDN      = intarr(run)
  nOptiUP      = intarr(run)
  nOptiDN      = intarr(run)

  outerUP_mask = intarr(run)
  outerDN_mask = intarr(run)
  interUP_mask = intarr(run)
  interDN_mask = intarr(run)
  inner_mask   = intarr(run)

  ;; Place to store the leight weighted extraction centroids
  innCtr   = fltarr(run)
  intCtrUP = fltarr(run)
  intCtrDN = fltarr(run)
  outCtrUP = fltarr(run)
  outCtrDN = fltarr(run)
  optCtrUP = fltarr(run)
  optCtrDN = fltarr(run)
  
  ;;  DO THE EXTRACTIONS
  for jj = 0, run - 1 do begin

     ;; Do the optimal extraction a la Horne 86, first w/o, then with masking
     
     opti[jj]     = total(modim[jj,*] * sci[jj,*] / var[jj,*]) $
                    / total(modim[jj,*]^2 / var[jj,*])

     opti_var[jj] = total(modim[jj,*]^2 / var[jj,*])^(-1)

     goods        = where(~mask[jj,*], ng)
     optiMask[jj] = total(modim[jj,goods]*spec[jj,goods]/var[jj,goods]) $
                    / total(modim[jj,goods]^2/var[jj,goods])

     optiMask_Var[jj] = total(modim[jj,goods]^2/var[jj,goods])^(-1)

     ;; Do the binned extractions
     all          = where( tx ge outerdn[jj] AND tx le outerup[jj], nall)
     integ[jj]    = total(sci[jj,all]) / nall
     ivar[jj]     = total(var[jj,all]) / nall^2
     nAlls[jj]    = nall
     
     outUP        = where(tx gt interup_bound[jj] AND tx le outerup_bound[jj], noutUP) ;; Upper outer region
     outDN        = where(tx ge outerdn_bound[jj] AND tx lt interdn_bound[jj], noutDN) ;; Lower outer region
     outerUP[jj]  = total(sci[jj,outUP]) / noutUP                                      ;; Upper outer extraction    
     outerDN[jj]  = total(sci[jj,outDN]) / noutDN                                      ;; Lower outer extraction    
     outvarUP[jj] = total(var[jj,outUP]) / noutUP^2                                    ;; Error^2 on upper outer extraction
     outvarDN[jj] = total(var[jj,outDN]) / noutDN^2                                    ;; Error^2 on lower outer extraction
     nOutsUP[jj]  = noutUP
     nOutsDN[jj]  = noutDN
     outerUP_mask[jj] = total(mask[jj,outUP])
     outerDN_mask[jj] = total(mask[jj,outDN])
     outCtrUP[jj] = total(tx[outUP] * modIm[jj, outUP]) / total(modIm[jj, outUP])           ;; Light-weighted extraction center.
     outCtrDN[jj] = total(tx[outDN] * modIm[jj, outDN]) / total(modIm[jj, outDN])
     
     intUP        = where(tx gt innerup_bound[jj] AND tx le interup_bound[jj], nintUP) ;; Upper intermediate region
     intDN        = where(tx ge interdn_bound[jj] AND tx lt innerdn_bound[jj], nintDN) ;; Lower intermediate region
     interUP[jj]  = total(sci[jj,intUP]) / nintUP                                      ;; Upper intermediate extraction
     interDN[jj]  = total(sci[jj,intDN]) / nintDN                                      ;; Lower intermediate extraction   
     intvarUP[jj] = total(var[jj,intUP]) / nintUP^2                                    ;; Error^2 on upper inter extraction
     intvarDN[jj] = total(var[jj,intDN]) / nintDN^2                                    ;; Error^2 on lower inter extraction
     nIntsUP[jj]  = nintUP
     nIntsDN[jj]  = nintDN
     interUP_mask[jj] = total(mask[jj,intUP])
     interDN_mask[jj] = total(mask[jj,intDN])
     intCtrUP[jj] = total(tx[intUP] * modIm[jj, intUP]) / total(modIm[jj, intUP])
     intCtrDN[jj] = total(tx[intDN] * modIm[jj, intDN]) / total(modIm[jj, intDN])

     in           = where(tx ge innerdn_bound[jj] AND tx le innerup_bound[jj], nin)    ;; Inner region
     inner[jj]    = total(sci[jj,in ]) / nin                                           ;; Inner extraction
     innvar[jj]   = total(var[jj,in ]) / nin^2                                         ;; Error^2 on inner extraction
     nInns[jj]    = nin
     inner_mask[jj] = total(mask[jj,in])
     innCtr[jj]     = total(tx[in] * modIm[jj, in]) / total(modIm[jj, in])

     optiOuterUP      = where( tx gt innerUp_bound[jj] AND tx le optiWideUp_bound[jj], noUP) ;; Upper optimal SN region
     optiOuterDn      = where( tx ge optiWideDn_bound[jj] AND tx lt innerDn_bound[jj], noDN) ;; Lower optimal SN region
     optiOutUP[jj]    = total(sci[jj,optiOuterUP]) / noUP                                    ;; Upper optimal SN extraction
     optiOutDN[jj]    = total(sci[jj,optiOuterDN]) / noDN                                    ;; Lower optimal SN extraction
     optiOutVarUP[jj] = total(var[jj,optiOuterUP]) / noUP^2                                  ;; Error^2 upper
     optiOutVarDN[jj] = total(var[jj,optiOuterDN]) / noDN^2                                  ;; Error^2 lower
     nOptiUP[jj]      = noUP
     nOptiDN[jj]      = noDN
     optCtrUp[jj]     = total(tx[optiOuterUP] * modIm[jj, optiOuterUP]) / total(modIm[jj, optiOuterUP])
     optCtrDn[jj]     = total(tx[optiOuterDN] * modIm[jj, optiOuterDN]) / total(modIm[jj, optiOuterDN])
     
     sky = mean([sci[jj,0:10], sci[jj,-11:*]])
     skies[jj]      = sky
      
     extrIm[jj,[optiOuterUP, $
                optiOuterDN]] = 5
     extrIm[jj,[outUP, $
                outDN]]       = 10
     extrIm[jj,[intUP, $
                intDN]]       = 20
     extrIm[jj,in]            = 40
     
  endfor       

  interExtrUpOff =      mean(intCtrUP - innCtr)
  interExtrDnOff = -1 * mean(intCtrDN - innCtr)
  outerExtrUpOff =      mean(outCtrUP - innCtr)
  outerExtrDnOff = -1 * mean(outCtrDN - innCtr)
  optiExtrUpOff  =      mean(optCtrUP - innCtr)
  optiExtrDNOff  = -1 * mean(optCtrDN - innCtr)
  
;  stop

  ;; Find the damn object
  dims = size(di, /dim)
  mwrfits, di, 'tmp.fits', /create 
  spawn, 'sex tmp.fits -c locateObjectOnStamp.sex'
  cat = mrdfits('tmp.cat', 2)
  if n_elements(cat) lt 1 then begin
     print, ''
     print, '!!! FAINT OBJECT !!! '
     print, ''
     spawn, 'sex tmp.fits -c locateObjectOnStamp_faint.sex'
     cat = mrdfits('tmp.cat', 2)
  endif
  xx = cat.X_IMAGE - 1 ;; to idl coords
  yy = cat.Y_IMAGE - 1
  dx = xx - (dims[0] / 2 - 1)
  dy = yy - (dims[1] / 2 - 1)
  dpos = sqrt(dx^2 + dy^2)
  obj = min(dpos, hit)
  ctr = yy[hit]
  re_Sex = cat[hit].FLUX_RADIUS
  if obj gt 10 then begin ;; in case object is actually a piece of some much brighter galaxy...
     ctr = dims[1] / 2 - 1
     re_Sex = -99
     badflag = 1
  endif
     
  dYoutUPUP  = (mean(outerUP_bound - trace) + ctr) < (height - 1)
  dYoutUPDN  = (mean(interUP_bound - trace) + ctr)
  dYoutDNDN  = (mean(outerDN_bound - trace) + ctr) > 0
  dYoutDNUP  = (mean(interDN_bound - trace) + ctr)
  if NOT keyword_Set(ACS) then begin
     x1 = 20
     x2 = 40
  endif else begin
     x1 = 10
     x2 = 30
  endelse
  lsfOuterUP = getLSF(di, x1 = x1, x2 = x2, $
                      y1 = dYoutUPDN, y2 = dYoutUPUP, acs = acs);, /debug)
  lsfOuterDN = getLSF(di, x1 = x1, x2 = x2, $               
                      y1 = dYoutDNDn, y2 = dYoutDNUP, acs = acs);, /debug)

  dYintUPUP  = mean(interUP_bound - trace) + ctr
  dYintUPDN  = mean(innerUP_bound - trace) + ctr
  dYintDNDN  = mean(interDN_bound - trace) + ctr
  dYintDNUP  = mean(innerDN_bound - trace) + ctr
  lsfInterUP = getLSF(di, x1 = x1, x2 = x2, $
                      y1 = dYintUPDN, y2 = dYintUPUP, acs = acs);, /debug)
  lsfInterDN = getLSF(di, x1 = x1, x2 = x2, $               
                      y1 = dYintDNDN, y2 = dYintDNUP, acs = acs);, /debug)

  dYinnUP = mean(innerUP_bound - trace) + ctr
  dYinnDN = mean(innerDN_bound - trace) + ctr
  lsfInner = getLSF(di, x1 = x1, x2 = x2, $
                    y1 = dYinnDN, y2 = dYinnUP, acs = acs);, /debug)
                                                      
;  inspect = 1
  if inspect then begin
     window, 0, xsize = 800, ysize = 1200
     !P.MULTI = [0,1,3]
     plot, lambda, findgen(60), /nodat, /xsty
     cgimage, spec, stretch = 2, /over
     plot, lambda, findgen(60), /nodat, /xsty
     cgimage, contam, stretch = 2, /over
     plot, lambda, findgen(60), /nodat, /xsty
     cgimage, sci, stretch = 2, /over
     print, ''
     print, 'ID:', sxpar(thead, 'POINTING'), '', sxpar(thead, 'ID')
     print, 'Continue?'
     k = get_kbrd(1)
     print, ''
     !P.MULTI = 0
  endif

;  stop
  
  ;;  WRITE OUTPUT
  savedata = {LAMBDA:             lambda, $
              Z:                  z, $
              Z_QUAL:             z_qual, $
              Z_PHOT:             z_phot, $
              RA:                 sxpar(thead, 'RA'), $
              DEC:                sxpar(thead, 'DEC'), $
              MAG:                mag, $
              F_TOT:              integ, $
              F_OPTI:             opti, $
              F_OPTI_MASK:        optiMask, $
              F_INNER:            inner, $
              F_INTER_DN:         interDN, $
              F_INTER_UP:         interUP, $
              F_OUTER_DN:         outerDN, $
              F_OUTER_UP:         outerUP, $
              F_OPTI_OUTER_DN:    optiOutDN, $
              F_OPTI_OUTER_UP:    optiOutUP, $
              VAR_TOT:            ivar, $
              VAR_OPTI:           opti_var, $
              VAR_OPTI_MASK:      optiMask_Var, $
              VAR_INNER:          innvar, $
              VAR_INTER_DN:       intvarDN, $
              VAR_INTER_UP:       intvarUP, $
              VAR_OUTER_DN:       outvarDN, $
              VAR_OUTER_UP:       outvarUP, $
              VAR_OPTI_OUTER_DN:  optiOutVarDN, $
              VAR_OPTI_OUTER_UP:  optiOutVarUP, $
              NPIX_TOT:           nAlls, $
              NPIX_INNER:         nInns, $
              NPIX_INTER_DN:      nIntsDN, $
              NPIX_INTER_UP:      nIntsUP, $
              NPIX_OUTER_DN:      nOutsDN, $
              NPIX_OUTER_UP:      nOutsUP, $
              NPIX_OPTI_OUTER_DN: nOptiDN, $
              NPIX_OPTI_OUTER_UP: nOptiUP, $
              MASK_OUTER_UP:      outerUP_mask, $
              MASK_OUTER_DN:      outerDN_mask, $
              MASK_INTER_UP:      interUP_mask, $
              MASK_INTER_DN:      interDN_mask, $
              MASK_INNER:         inner_mask, $
              SENSITIVITY:        sens, $
              EXPTIME:            exptime, $
              FILENAME:           infile, $
              MODEL_TRACE:        modim, $
              EXTRACT_ZONES:      extrIM, $
              EXTR_MDPT_INNER:    innCtr, $
              EXTR_INTUP_OFF:     interExtrUpOff, $
              EXTR_INTDN_OFF:     interExtrDnOff, $
              EXTR_OUTUP_OFF:     outerExtrUpOff, $
              EXTR_OUTDN_OFF:     outerExtrDnOff, $
              EXTR_OPTUP_OFF:     optiExtrUpOff, $
              EXTR_OPTDN_OFF:     optiExtrDnOff, $              
              LSF_TOT:            prof.LSF, $
              LSF_INNER:          lsfInner.LSF, $
              LSF_INTER_UP:       lsfInterUP.LSF, $
              LSF_INTER_DN:       lsfInterDN.LSF, $
              LSF_OUTER_UP:       lsfOuterUP.LSF, $
              LSF_OUTER_DN:       lsfOuterDN.LSF, $
              LSF_TOT_PARS:       prof.LSF_PARAMS, $
              LSF_INNER_PARS:     lsfInner.LSF_PARAMS, $
              LSF_INTER_UP_PARS:  lsfInterUP.LSF_PARAMS, $
              LSF_INTER_DN_PARS:  lsfInterDN.LSF_PARAMS, $
              LSF_OUTER_UP_PARS:  lsfOuterUP.LSF_PARAMS, $
              LSF_OUTER_DN_PARS:  lsfOuterDN.LSF_PARAMS, $
              RE:                 prof.RE, $
              RE_SEX:             re_sex, $
              SPECORIG:           spec, $
              MODEL2D:            model, $
              COTAM_FRACT_MAP:    fcmap, $
              SPEC2D:             sci, $
              DIRECT_IM:          di, $
              TRACE:              trace, $
              INNERUP:            innerup_bound, $
              INNERDN:            innerdn_bound, $
              INTERUP:            interup_bound, $
              INTERDN:            interdn_bound, $
              OUTERUP:            outerup_bound, $
              OUTERDN:            outerdn_bound, $
              OPTI_OUTERUP:       optiWideUp_bound, $
              OPTI_OUTERDN:       optiWideDn_bound, $
              PA:                 pa, $
              MASK:               mask, $
              BADFLAG:            badflag}
;  mwrfits, savedata, pa+'_extractions/'+string(ii, f = '(I04)')+'_radStuff.fits', /create
  
  RETURN, savedata
end
