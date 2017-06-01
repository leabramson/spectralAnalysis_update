pro writeforpyspecfit, inspec, output, $
                       CHIP = chip, $
                       REGION = region, $
                       MASKMORE = maskmore;, $
;                       ACS = acs

  if NOT keyword_Set(MASKMORE) then maskmore = 0 else maskmore = 1
;  if NOT keyword_Set(ACS) then acs = 0 else acs = 1
  
  data = mrdfits(inspec,1,/silent)
  lambda = data.LAMBDA
  case region of
     'inner': begin
        spec   = data.F_INNER / data.SENSITIVITY
        err    = sqrt(data.VAR_INNER) / data.SENSITIVITY
        msk    = data.MASK_INNER 
        lsf    = data.LSF_INNER
        lsfpar = data.LSF_INNER_PARS
     end
     'interUP': begin
        spec   = data.F_INTER_UP / data.SENSITIVITY
        err    = sqrt(data.VAR_INTER_UP) / data.SENSITIVITY
        msk    = data.MASK_INTER_UP
        lsf    = data.LSF_INTER_UP
        lsfpar = data.LSF_INTER_UP_PARS
     end
     'interDN': begin
        spec   = data.F_INTER_DN / data.SENSITIVITY
        err    = sqrt(data.VAR_INTER_DN) / data.SENSITIVITY
        msk    = data.MASK_INTER_DN
        lsf    = data.LSF_INTER_DN
        lsfpar = data.LSF_INTER_DN_PARS
     end
     'outerUP': begin
        spec   = data.F_OUTER_UP / data.SENSITIVITY
        err    = sqrt(data.VAR_OUTER_UP) / data.SENSITIVITY 
        msk    = data.MASK_OUTER_UP
        lsf    = data.LSF_OUTER_UP
        lsfpar = data.LSF_OUTER_UP_PARS
     end
     'outerDN': begin
        spec   = data.F_OUTER_DN / data.SENSITIVITY
        err    = sqrt(data.VAR_OUTER_DN) / data.SENSITIVITY 
        msk    = data.MASK_OUTER_DN
        lsf    = data.LSF_OUTER_DN
        lsfpar = data.LSF_OUTER_DN_PARS
     end
     'optiTOT': begin
        spec   = data.F_OPTI / data.SENSITIVITY
        err    = sqrt(data.VAR_OPTI) / data.SENSITIVITY
        msk    = data.MASK_OUTER_DN + data.MASK_OUTER_UP + $
                 data.MASK_INTER_DN + data.MASK_INTER_UP + $
                 data.MASK_INNER
        lsf    = data.LSF_TOT
        lsfpar = data.LSF_TOT_PARS
     end
  endcase

  ;; CENTER THE LSF
  lsf = [lsf,lsf[-1]]
;  lsf /= total(lsf)
;  foo = max(lsf, pk)
;  del = (n_elements(lsf) - 1) / 2 - pk
;  lsf = shift(lsf, del) ;; center the LSF
  run = indgen(n_elements(lsf)) - (n_elements(lsf) - 1) / 2
  u      = run / lsfpar[2] 
  newLSF = lsfpar[0] / (u^2 + 1)^lsfpar[3] 
  newLSF /= total(newLSF)  ;; new, centered LSF

  mask   = replicate(1, n_elements(lambda))
  
  if CHIP eq 'blue' then begin
     qui = where(lambda lt 8400 OR lambda gt 1.12d4)  ;;  hard limits from vis. inspection
     mask[qui] = 0
;     qui = where(lambda gt 1.12d4)
;     mask[qui] = 0
  endif else if CHIP eq 'red' then begin
     qui = where(lambda lt 11000 OR lambda gt 1.650d4)
     mask[qui] = 0
;     qui = where(lambda gt 1.650d4)
;     mask[qui] = 0
  endif else if chip eq 'ACS' then begin
     qui = where(lambda lt 5512 OR lambda gt 9522)
     mask[qui] = 0
  endif

  if MASKMORE then begin
     mask[where(msk gt 0)] = 0
     mask[where(spec le 0 OR err le 0)] = 0     
  endif
     
;  plot, lambda, spec, yran = [0,400]
;  oplot, lambda[where(mask)], spec[where(mask)], col = 255
;  stop
  
  close, 1
  openw, 1, output, width = 128
  printf, 1, '#WAVELENGTH FLUX ERROR MASK'
  for ii = 0, n_elements(lambda) - 1 do $
     printf, 1, lambda[ii], spec[ii], err[ii], mask[ii]
  close, 1
  
  close, 1
  openw, 1, repstr(output, '.pyspec', '.lsf'), width = 128
  printf, 1, '#CTR_OFFSET LSF'
  for ii = 0, n_elements(lsf) - 1 do $
     printf, 1, run[ii], lsf[ii]
  close, 1

  print, ''
  print, 'ONE-D SPECTRUM FOR '+inspec+', REGION '+region+' OUTPUT TO : '+output
  print, ''
  
end

;;
;;
;;

;pro doRoundOne, inlist, FIELD = field
;
;  readcol, inlist, files, ids, pas, grisms, $
;           f = 'A,A,A,A', comment = '#'
;  tid = ids[sort(ids)]
;  tid = tid[UNIQ(tid)]
;  ngals = n_elements(tid)
;
;  for ii = 0, ngals - 1 do begin
;     toget = where(ids eq tid[ii], nfiles) ;; Grab all the files for this galaxy
;     tfiles  = files[toget] 
;     tpas    = pas[toget]
;     tgrisms = grisms[toget]
;
;     for jj = 0, nfiles - 1 do begin
;        case tgrisms[jj] of
;           'B': chip = 'blue'
;           'R': chip = 'red'
;        endcase
;        ;; write them out
;;        outname =
;;        'spectralAnalysis/'+field+'/'+tid[ii]+'_'+tpas[jj]+'_'+tgrisms[jj] ;; Different dir. structure for faint objects...
;        outname = field+'/'+tid[ii]+'_'+tpas[jj]+'_'+tgrisms[jj]
;        writeForPyspecfit, tfiles[jj], outname+'_inner.pyspec', $
;                           chip = chip, region = 'inner'
;        writeForPyspecfit, tfiles[jj], outname+'_inter.pyspec', $
;                           chip = chip, region = 'inter'
;        writeForPyspecfit, tfiles[jj], outname+'_outer.pyspec', $
;                           chip = chip, region = 'outer'
;        print, ''
;        print, '1D spectra written to : '+outname+'_ ... '
;     endfor;
;
;  endfor
;     
;end
; DOROUNDONE, 'spectralAnalysis/M1149_extractions.list', field = 'M1149'
; DOROUNDONE, 'spectralAnalysis/M2129_extractions.list', field = 'M2129'
; DOROUNDONE, 'spectralAnalysis/M0717_extractions.list', field = 'M0717'

; DOROUNDONE, 'M1149_faint_extractions.list', field = 'M1149_faint'
; DOROUNDONE, 'M2129_faint_extractions.list', field = 'M2129_faint'
; DOROUNDONE, 'M0717_faint_extractions.list', field = 'M0717_faint'
; DOROUNDONE, 'R1347_faint_extractions.list', field = 'R1347_faint'
; DOROUNDONE, 'M1423_faint_extractions.list', field = 'M1423_faint'

