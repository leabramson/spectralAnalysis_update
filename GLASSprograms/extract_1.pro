function extract_1, infile, z, $
                    Z_PHOT  = z_phot, $
                    PA      = pa, $
                    DOMASK  = domask, $
                    REFIT   = refit, $
                    INSPECT = inspect

  if N_ELEMENTS(z) eq 0 then z = -99
  if NOT keyword_set(Z_PHOT) then z_phot = -99
  if NOT keyword_SET(PA) then pa = -99
  if NOT keyword_set(DOMASK) then domask = 0 else domask = 1 ;; mask rather than subtract contamination
  if NOT keyword_set(REFIT) then refit = 0 else refit = 1    ;; refit continuum centerline
  if NOT keyword_set(INSPECT) then inspect = 0 else inspect = 1
  
  ;; READ CRAP
  null   = mrdfits(infile, 0, thead, /silent)
  spec   = mrdfits(infile, 'sci', head, /silent)
  contam = mrdfits(infile, 'contam', /silent)
  sens   = mrdfits(infile, 'sens', /silent)
  trace  = mrdfits(infile, 'ytrace', /silent)
  err    = mrdfits(infile, 'wht', /silent)
  di     = mrdfits(infile, 'dsci', /silent)
  
  ;; BASIC INFO
  exptime = sxpar(thead, 'EXPTIME')   ;; SECONDS
  mag     = sxpar(thead, 'MAG')       ;; Probably F140W?

  ;; TO POISSON UNITS
  spec   *= exptime
  contam *= exptime
  err[where(err le 0)] = median(err) ;; to avoid NaNs..
  var     = (err * exptime)^2
  
  ;; GET SPATIAL PROFILE/LSF
  prof    = getProfile(di)
  profile = prof.LSF       
  re      = ceil(prof.RE)        ;; Returns a pseudo "half-light radius"
  params  = prof.SPATIAL_PARAMS
  
  ;; MAKE A WAVELENGTH AXIS
  l0  = sxpar(head, 'CRVAL1')
  dl  = sxpar(head, 'CD1_1')
  run = sxpar(head, 'NAXIS1')
  lambda = findgen(run) * dl + l0

  yctr = sxpar(head, 'CRPIX2') - 1
  tx = findgen(sxpar(head, 'NAXIS2'))
  height = n_elements(tx)

  ;; MAKE A MASK IMAGE
  mask = contam
  mask[where(contam ge 2 * sqrt(var), compl = good)] = 1
  mask[good] = 0                                ;; Make the mask
  mask = smooth(mask, [15,5], /edge_wrap, /nan) ;; Grow the mask
  mask[where(mask ge 1d-4, compl = good)] = 1
  mask[good] = 0
  
  ;; MAKE A SCIENCE IMAGE
  sci   = (spec - contam)

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
     tofit = use[where(abs(ctroid[0,use] - median(ctroid[0,use])) $
                       le 2 * stddev(ctroid[0,use]))]
     cline = poly_fit(lambda[tofit], ctroid[0,tofit], $
                      measure_err = ctroid[2, tofit], 1)
     trace = cline[0] + lambda * cline[1]
  endelse

  ;; MAKE A MODEL IMAGE FOR OPTIMAL EXTRACTION
  modIm = sci  ;; Will contain the profile for optimal extraction
  for jj = 0, n_elements(trace) - 1 do begin
     u = (tx - trace[jj]) / params[2]
     modIm[jj,*] = 1. / (u^2 + 1)^params[3]
     modIm[jj,*] /= total(modIm[jj,*])
  endfor
  extrIm = sci
  
  ;;  DEFINE EXTRACTION RADIAL REGIMES
  innerup = ceil ( trace   + 0.25 * re )
  innerdn = floor( trace   - 0.25 * re )
  interup = ceil ( innerup + 0.50 * re )
  interdn = floor( innerdn - 0.50 * re )
  outerup = ceil ( interup + 0.75 * re )
  outerdn = floor( interdn - 0.75 * re )

  ;; Get an optimally wide "outer aperture"
  optiWide = getoptizone(sci, var, $
                         clineUp = innerUp, $
                         clineDn = innerDn)
  optiWideUp = innerUp + optiWide
  optiWideDn = innerDn - optiWide
  
  skies  = dblarr(run)
  integ  = dblarr(run)
  outer  = dblarr(run)
  inter  = dblarr(run)
  inner  = dblarr(run)
  optiOut = dblarr(run)
  
  ivar   = dblarr(run)
  ovar   = dblarr(run)
  intvar = dblarr(run)
  innvar = dblarr(run)
  optiOutVar = dblarr(run)
  
  nAlls = intarr(run)
  nInns = intarr(run)
  nInts = intarr(run)
  nOuts = intarr(run)
  nOpti = intarr(run)
  
  ;;  DO THE EXTRACTIONS
  opti         = dblarr(run)
  opti_Var     = dblarr(run)
  optiMask     = dblarr(run)
  optiMask_Var = dblarr(run)
  for jj = 0, run - 1 do begin

     ;; Do the optimal extraction

;     if jj eq floor(run/3. * 2) then stop
     
     opti[jj]     = total(modim[jj,*] * sci[jj,*] / var[jj,*]) $
                    / total(modim[jj,*]^2 / var[jj,*])

     opti_var[jj] = total(modim[jj,*]^2 / var[jj,*])^(-1)

     goods    = where(~mask[jj,*], ng)
     optiMask[jj] = total(modim[jj,goods]*spec[jj,goods]/var[jj,goods]) $
                    / total(modim[jj,goods]^2/var[jj,goods])

     optiMask_Var[jj] = total(modim[jj,goods]^2/var[jj,goods])^(-1)

     ;; Do the binned extractions
     all = where( tx ge outerdn[jj] AND tx le outerup[jj], nall)
     out = where((tx gt interup[jj] AND tx le outerup[jj]) $
                 OR $
                 (tx ge outerdn[jj] AND tx lt interdn[jj]), nout)
     int = where((tx gt innerup[jj] AND tx le interup[jj]) $
                 OR $
                 (tx ge interdn[jj] AND tx lt innerdn[jj]), nint)
     in  = where(tx ge innerdn[jj] AND tx le innerup[jj], nin)

     optiOuter = where((tx gt innerUp[jj] AND tx le optiWideUp[jj]) $
                       OR $
                       (tx ge optiWideDn[jj] AND tx lt innerDn[jj]), no)
     
     sky = mean([sci[jj,0:10], sci[jj,-11:*]])

     skies[jj]      = sky
     integ[jj]      = total(sci[jj,all]) / nall
     outer[jj]      = total(sci[jj,out]) / nout
     inter[jj]      = total(sci[jj,int]) / nint
     inner[jj]      = total(sci[jj,in ]) / nin 
     optiOut[jj]    = total(sci[jj,optiOuter]) / no
     
     ivar[jj]       = total(var[jj,all]) / nall^2
     ovar[jj]       = total(var[jj,out]) / nout^2
     intvar[jj]     = total(var[jj,int]) / nint^2
     innvar[jj]     = total(var[jj,in ]) / nin^2
     optiOutVar[jj] = total(var[jj,optiOuter]) / no^2
     
     nAlls[jj]      = nall
     nInns[jj]      = nin
     nInts[jj]      = nint
     nOuts[jj]      = nout
     nOpti[jj]      = no

     extrIm[jj,optiOuter] = 5
     extrIm[jj,out]       = 10
     extrIm[jj,int]       = 20
     extrIm[jj,in ]       = 40
     
  endfor       

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
  savedata = {LAMBDA:         lambda, $
              Z:              z, $
              Z_PHOT:         z_phot, $
              RA:             sxpar(thead, 'RA'), $
              DEC:            sxpar(thead, 'DEC'), $
              MAG:            mag, $
              F_TOT:          integ, $
              F_OPTI:         opti, $
              F_INNER:        inner, $
              F_INTER:        inter, $
              F_OUTER:        outer, $
              F_OPTI_OUTER:   optiOut, $
              VAR_TOT:        ivar, $
              VAR_OPTI:       opti_var, $
              VAR_INNER:      innvar, $
              VAR_INTER:      intvar, $
              VAR_OUTER:      ovar, $
              VAR_OPTI_OUTER: optiOutVar, $
              NPIX_TOT:       nAlls, $
              NPIX_INNER:     nInns, $
              NPIX_INTER:     nInts, $
              NPIX_OUTER:     nOuts, $
              NPIX_OPTI_OUTER:nOpti, $
              SENSITIVITY:    sens, $
              EXPTIME:        exptime, $
              FILENAME:       infile, $
              MODEL_TRACE:    modim, $
              EXTRACT_ZONES:  extrIM, $
;              EXTRACT_MIDPTS: , $
              LSF:            prof.LSF, $
              RE:             prof.RE, $
              SPEC_ORIG:      spec, $
              SPEC2D:         sci, $
              DIRECT_IM:      di, $
              TRACE:          trace, $
              INNERUP:        innerup, $
              INNERDN:        innerdn, $
              INTERUP:        interup, $
              INTERDN:        interdn, $
              OUTERUP:        outerup, $
              OUTERDN:        outerdn, $
              OPTI_OUTERUP:   optiWideUp, $
              OPTI_OUTERDN:   optiWideDn, $
              PA:             pa, $
              MASK:           mask}
;  mwrfits, savedata, pa+'_extractions/'+string(ii, f = '(I04)')+'_radStuff.fits', /create
  
  RETURN, savedata
end
