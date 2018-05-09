function calcbt, infile, mag, emag, $
                 NITER = niter

  if NOT keyword_set(NITER) then niter = 100

  data = infile
  
  uq = where(data.REGION eq 'INNFL' OR $
             data.REGION eq 'INTER' OR $
             data.REGION eq 'OUTER')

  q = where(data.REGION eq 'OPTFL')
  
  ro = data[0].RE_OBS
  rp = data[0].RE_PHYS
  r  = (mrdfits(data[0].FITSFILE, 1, /silent)).RE
  rs = (mrdfits(data[0].FITSFILE, 1, /silent)).RE_SEX
  
  pscale = ro / r
  kscale = rp / r
  
  ;; Correct for magnification
  data.LMASS -= alog10(mag)
  data.AREA_PHYS /= mag
  
  masses  = data.LMASS[0]
  emasses = sqrt(total((0.5 * (data[uq].LMASS[2] - data[uq].LMASS[1]))^2) + emag^2)
  
  areas  = alog10(data.AREA_PHYS)
  eareas = emag
  
  sizes  = alog10(rp / sqrt(mag))
  esizes = 0.5 * emag
  
  ssizes  = alog10(rs * kscale / sqrt(mag))
  essizes = 0.5 * emag
  
  ;; Estimate B/T
  ;; Mass w/in r_e over mass w/in 2 * re
  tdat = data[uq]
  tdat = tdat[sort(tdat.RNORM)]
  mmaster = data[q].LMASS[0]
  emmaster = 0.5 * (data[q].LMASS[2] - data[1].LMASS[1])
  m    = tdat.LMASS[0] - areas
  em   = 0.5 * (tdat.LMASS[2] - tdat.LMASS[1]);)^2 ;sqrt(+ eareas^2)
  tr   = tdat.RNORM
  dr = 0.02
  rad = findgen(2./dr+1) * dr
  mprofs = fltarr(n_elements(rad), niter)
  tbts   =  fltarr(niter,2) ;; LEA 20170801 added after tweaking plotsizemass.pro to try to combat profile-integral/total mass offsets
  for ll = 0, niter - 1 do begin
     tm = m + randomn(seed, 3) * em
     t = interpol(tm, tr, rad)     
     mprofs[*,ll] = total(2 * !pi * rad * dr * rp^2 * 10.^t / mag, /cum)
     mprofs[*,ll] /= mprofs[-1,ll] / 10.^(mmaster + randomn(seed, 1 ) * emmaster)[0]
     tbts[ll,*] = [mprofs[value_locate(rad * rp, 2.5),ll] / mprofs[-1,ll], $      ;; M(r < 2.5 kpc) from Morishita+15
                   mprofs[value_locate(rad * rp, 2.5),ll] / 10.^data[q].LMASS[0]]
  endfor

;  bts[ii,jj]  = 0.5 * (median(tbts[*,0]) + median(tbts[*,1]))
;  ebts[ii,jj] = sqrt((0.5 * abs(median(tbts[*,0]) - median(tbts[*,1])))^2 + $
;                     0.5 * (variance(tbts[*,0], /nan) + variance(tbts[*,1], /nan)))
  bts  = median(tbts[*,0])
  ebts = stddev(tbts[*,0], /nan)

  RETURN, [bts, ebts]
end
