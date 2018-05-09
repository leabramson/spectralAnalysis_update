;; Take output from cutOnZ and cutOnMag and extract spectra.

; do1field, 'uniqueAllCL_1.0_1.8_0_21.5.list', 'MACS1149', outdir = 'testM1149'
; do1field, 'uniqueAllCL_1.0_1.8_0_21.5.list', 'MACS0717', outdir = 'testM0717'
; do1field, 'uniqueAllCL_1.0_1.8_0_21.5.list', 'MACS1423', outdir = 'testM1423'
; do1field, 'uniqueAllCL_1.0_1.8_0_21.5.list', 'MACS2129', outdir = 'testM2129'
; do1field, 'uniqueAllCL_1.0_1.8_0_21.5.list', 'RXJ1347',  outdir = 'testR1347'
pro do1field, inlist, field, $
              OUTDIR = outdir

  if NOT keyword_set(OUTDIR) then $
     outdir = './' $
  else begin
     foo = file_test(outdir, /dir)
     if NOT foo then $
        spawn, 'mkdir '+outdir
  endelse
  
  toGet = dump1field(inlist, field)

  for ii = 0, n_elements(toGet.Z) - 1 do begin

     check = file_info(toget.PA1B[ii])
     if check.exists then begin
        toWrite1b = extract_1(toget.PA1B[ii], toGet.Z[ii], PA = toget.PA1) ;, /refit)
        toWrite2b = extract_1(toget.PA2B[ii], toGet.Z[ii], PA = toget.PA2) ;, /refit)
        toWrite1r = extract_1(toget.PA1R[ii], toGet.Z[ii], PA = toget.PA1) ;, /refit)
        toWrite2r = extract_1(toget.PA2R[ii], toGet.Z[ii], PA = toget.PA2) ;, /refit)

        mwrfits, toWrite1b, outdir+'/'+strcompress(string(toGet.ID[ii])+'_1_B.fits', /rem), /create
        mwrfits, toWrite2b, outdir+'/'+strcompress(string(toGet.ID[ii])+'_2_B.fits', /rem), /create
        mwrfits, toWrite1r, outdir+'/'+strcompress(string(toGet.ID[ii])+'_1_R.fits', /rem), /create
        mwrfits, toWrite2r, outdir+'/'+strcompress(string(toGet.ID[ii])+'_2_R.fits', /rem), /create
     endif

  endfor

end

pro writeOut, inspec, output, CHIP = chip, $
              region = region
  
  data = mrdfits(inspec,1)
  lambda = data.LAMBDA
  case region of
     'inner': begin
        spec   = data.F_INNER / data.SENSITIVITY         ;; S/N = data.F_INNER / sqrt(data.VAR_INNER) in Poisson units
        err    = sqrt(data.VAR_INNER) / data.SENSITIVITY ;; spec / err maintains this ratio
     end
     'inter': begin
        spec   = data.F_INTER / data.SENSITIVITY
        err    = sqrt(data.VAR_INTER) / data.SENSITIVITY
     end
     'outer': begin
        spec   = data.F_OUTER / data.SENSITIVITY
        err    = sqrt(data.VAR_OUTER) / data.SENSITIVITY
     end
  endcase
  mask   = replicate(1, n_elements(lambda))

  if CHIP eq 'blue' then begin
     qui = where(lambda lt 8150)
     mask[qui] = 0
     qui = where(lambda gt 1.12d4)
     mask[qui] = 0
  endif else begin
     qui = where(lambda lt 11000)
     mask[qui] = 0
     qui = where(lambda gt 1.675d4)
     mask[qui] = 0
  endelse

  plot, lambda, spec, yran = [0,400]
  oplot, lambda[where(mask)], spec[where(mask)], col = 255

;  stop
  
  close, 1
  openw, 1, output, width = 128
  printf, 1, '#WAVELENGTH FLUX ERROR MASK'
  for ii = 0, n_elements(lambda) - 1 do $
     printf, 1, lambda[ii], spec[ii], err[ii], mask[ii]
  close, 1

  lsf = data.LSF
  lsf = [lsf,lsf[-1]]
  lsf /= total(lsf)
  run = indgen(n_elements(lsf)) - (n_elements(lsf) - 1) / 2
  
  close, 1
  openw, 1, 'LSF_'+output, width = 128
  printf, 1, '#CTR_OFFSET LSF'
  for ii = 0, n_elements(lsf) - 1 do $
     printf, 1, run[ii], lsf[ii]
  close, 1
  

end

pro writeAll
  writeOut, '753_1_B.fits', '753_1_B_INNER_FOR_DREW.dat', $
            chip = 'blue', region = 'inner'
  writeOut, '753_1_R.fits', '753_1_R_INNER_FOR_DREW.dat', $
            chip = 'red', region = 'inner'

  writeOut, '753_1_B.fits', '753_1_B_INTER_FOR_DREW.dat', $
            chip = 'blue', region = 'inter'
  writeOut, '753_1_R.fits', '753_1_R_INTER_FOR_DREW.dat', $
            chip = 'red', region = 'inter'

  writeOut, '753_1_B.fits', '753_1_B_OUTER_FOR_DREW.dat', $
            chip = 'blue', region = 'outer'
  writeOut, '753_1_R.fits', '753_1_R_OUTER_FOR_DREW.dat', $
            chip = 'red', region = 'outer'
end
