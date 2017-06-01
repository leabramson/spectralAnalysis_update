pro pyspecfitonetracelgnml, id, pa, region, z, $
                            PHOTFILE = photfile, $
                            ncores = ncores, $
                            OUTDIR = outdir, $
                            PYFILE = pyfile

  on_error, 2
  
  if NOT keyword_Set(NCORES) then ncores = 10
  if NOT keyword_set(PYFILE) then pyfile = 'fitspecLgnml.py'
  if NOT keyword_set(OUTDIR) then outdir = ''
  if NOT keyword_set(PHOTFILE) then photfile = 'None'
  
  pyprefix = '/home/labramson/Projects/GLASS/spectralAnalysis/'
  pyfile = pyprefix+pyfile

  instructions = strcompress("mpirun -np "+string(ncores, f = '(I02)')+$
                             " python "+pyfile+" "+$
                             string(id, f = '(I05)')+"_"+string(pa, f = '(I1)')+" "+$
                             string(region)+" "+$
                             string(z)+" "+$
;                             string(getage(z) * 1d9)+" "+$
                             photfile+" "+outdir)
  
  print, f = '(%"\n\n\nSending %s spectrum of source %i PA %i to PYSPECFIT ... \n")', $
  region, id, pa
  print, ' >>> '+instructions
  print, ''
  
  spawn, instructions

  
end

;;
;;
;;

;;
;;
;;

pro plot1res, blueRes, redRes, redshift

  readcol, blueRes, $
           lambdaB, innTraceB, innEerrB, innPolyB, innModelB, innUsedB
  readcol, redRes, $
           lambdaR, innTraceR, innEerrR, innPolyR, innModelR, innUsedR

  innSpecB = innTraceB * innPolyB
  innSpecB[where(~innUsedB)] = !values.F_NAN

  innErrB = innEerrB * innPolyB
  innErrB[where(~innUsedB)] = !values.F_NAN
  
  innModelB *= innPolyB
  innModelB[where(~innUsedB)] = !values.F_NAN

  innSpecR = innTraceR * innPolyR
  innSpecR[where(~innUsedR)] = !values.F_NAN

  innErrR = innEerrR * innPolyR
  innErrR[where(~innUsedR)] = !values.F_NAN
  
  innModelR *= innPolyR
  innModelR[where(~innUsedR)] = !values.F_NAN

;  lambda = [lambdaB, lambdaR] / (1 + redshift)
;  inner  = [innSpecB, innSpecR]
;  innerErr = [innErrB, innErrR]
;  innMod = [innModelB, innModelR]

  plot, lambdaB  / (1 + redshift), innSpecB, yran = [0,400], xran = [3000,8000], /nodat
  oplot, lambdaB / (1 + redshift), innSpecB
  oplot, lambdaR / (1 + redshift), innSpecR
  oplot, lambdaB / (1 + redshift), innErrB, col = '777777'x
  oplot, lambdaR / (1 + redshift), innErrR, col = '777777'x
  oplot, lambdaB / (1 + redshift), innModelB, col = 255, thick = 1
  oplot, lambdaR / (1 + redshift), innModelR, col = 255, thick = 1
  oplot, replicate(4000, 2), !Y.CRANGE, col = '555555'x, linesty = 1
  oplot, replicate(4863, 2), !Y.CRANGE, col = '555555'x, linesty = 1
  oplot, replicate(5007, 2), !Y.CRANGE, col = '555555'x, linesty = 1
  oplot, replicate(6563, 2), !Y.CRANGE, col = '555555'x, linesty = 1
  oplot, replicate(6720, 2), !Y.CRANGE, col = '555555'x, linesty = 1
  
end
