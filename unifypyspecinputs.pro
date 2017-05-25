pro unifypyspecinputs, blueside, redside, output

  readcol, blueside, $
           lb, fb , eb, mb, $
           f = 'F,F,F,F', $
           /quick, /silent, comment = '#'

  readcol, redside, $
           lr, fr, er, mr, $
           f = 'F,F,F,F', $
           /quick, /silent, comment = '#'

  dl = 2.
  lrun   = ceil((max(lr) - min(lb)) / dl)
  lambda = findgen(lrun) * dl + min(lb)
  bs     = where(lambda le max(lb))
  rs     = where(lambda ge min(lr))

  ub = where(mb)
  ur = where(mr)
  
  fb = interpol(fb, lb, lambda)
  eb = interpol(eb, lb, lambda)
  mb = interpol(mb, lb, lambda)
  
  fr = interpol(fr, lr, lambda)
  er = interpol(er, lr, lambda)
  mr = interpol(mr, lr, lambda)

  fb[where(lambda gt 12000)] = 0
;  fr[where(lambda lt min(lr[ur]) OR lambda gt 18000)] = 0
  
  outf = fltarr(lrun)
  oute = fltarr(lrun)
  outm = mb + mr

  for ii = 0, lrun - 1 do begin
     outf[ii] = (fb[ii] / eb[ii]^2 + fr[ii] / er[ii]^2) $
                   / (1. / eb[ii]^2 + 1. / er[ii]^2)
     oute[ii] = sqrt(1. / (1. / eb[ii]^2 + 1. / er[ii]^2))
  end

  outc = outm
  outm[where(outm gt 0)] = 1
  
  savedata = {LAMBDA:   lambda, $
              FLUX:     outf, $
              STDDEV:   oute, $
              COVERAGE: outc, $
              MASK:     outm }
  mwrfits, savedata, output, /create
  
end
