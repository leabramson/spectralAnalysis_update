pro pyspecfitonetrace, id, pa, region, z, zran, $
                       PHOTFILE = photfile, $
                       ncores = ncores, $
                       zlock = zlock, $
                       OUTDIR = outdir, $
                       PYFILE = pyfile

  on_error, 2
  
  if NOT keyword_Set(NCORES) then ncores = 10
  if NOT keyword_Set(ZLOCK ) then $
     zlock  = 0 $
  else $
     zlock = 1
  if NOT keyword_set(PYFILE) then pyfile = 'fitspecSolarPhot.py'
  if NOT keyword_set(OUTDIR) then outdir = ''
  if NOT keyword_set(PHOTFILE) then photfile = 'None'
  
  pyprefix = '/home/labramson/Projects/GLASS/spectralAnalysis/'
  pyfile = pyprefix+pyfile

  if zran eq 0 then zlock = 1
  
  if NOT zlock then $
     instructions = strcompress("mpirun -np "+string(ncores, f = '(I02)')+$
                                " python "+pyfile+" "+$
                                string(id, f = '(I05)')+"_"+string(pa, f = '(I1)')+" "+$
                                string(region)+" "+$
                                string(z)+" "+$
                                string(zran)+" "+photfile+" "+outdir) $
  else $
     instructions = strcompress("mpirun -np "+string(ncores, f = '(I02)')+$
                                " python "+pyfile+" "+$
                                string(id, f = '(I05)')+"_"+string(pa, f = '(I1)')+" "+$
                                string(region)+" "+$
                                string(z)+" "+$
                                string(0)+" "+photfile+" "+outdir)

  print, ''
  print, ''
  print, ''
  print, "Sending "+region+" spectrum of source "+id+" PA "+pa+" to PYSPECFIT ... "
  print, ''
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
