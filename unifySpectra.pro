pro unifySpectra, blueSpecFile, redSpecFile, $
                  OUTPUT = output

  bluData = mrdfits(blueSpecFile, 1)
  redData = mrdfits(redSpecFile, 1)

  bl            = bluData.LAMBDA
  binner        = bluData.F_INNER        / bluData.SENSITIVITY
  binterUP      = bluData.F_INTER_UP     / bluData.SENSITIVITY
  binterDN      = bluData.F_INTER_DN     / bluData.SENSITIVITY
  bouterUP      = bluData.F_OUTER_UP     / bluData.SENSITIVITY
  bouterDN      = bluData.F_OUTER_DN     / bluData.SENSITIVITY
;  bouter_opti   = bluData.F_OPTI_OUTER   / bluData.SENSITIVITY
  binn_err      = bluData.VAR_INNER      / bluData.SENSITIVITY^2
  bint_errUP    = bluData.VAR_INTER_UP   / bluData.SENSITIVITY^2
  bint_errDN    = bluData.VAR_INTER_DN   / bluData.SENSITIVITY^2
  bout_errUP    = bluData.VAR_OUTER_UP   / bluData.SENSITIVITY^2
  bout_errDN    = bluData.VAR_OUTER_DN   / bluData.SENSITIVITY^2
;  bout_opti_err = bluData.VAR_OPTI_OUTER / bluData.SENSITIVITY^2
  
  rl            = redData.LAMBDA
  rinner        = redData.F_INNER        / redData.SENSITIVITY
  rinterUP      = redData.F_INTER_UP     / redData.SENSITIVITY
  rinterDN      = redData.F_INTER_DN     / redData.SENSITIVITY
  routerUP      = redData.F_OUTER_UP     / redData.SENSITIVITY
  routerDN      = redData.F_OUTER_DN     / redData.SENSITIVITY
;  router_opti   = redData.F_OPTI_OUTER   / redData.SENSITIVITY
  rinn_err      = redData.VAR_INNER      / redData.SENSITIVITY^2
  rint_errUP    = redData.VAR_INTER_UP   / redData.SENSITIVITY^2
  rint_errDN    = redData.VAR_INTER_DN   / redData.SENSITIVITY^2
  rout_errUP    = redData.VAR_OUTER_UP   / redData.SENSITIVITY^2
  rout_errDN    = redData.VAR_OUTER_DN   / redData.SENSITIVITY^2
;  rout_opti_err = redData.VAR_OPTI_OUTER / redData.SENSITIVITY^2
  
  dl = 2.
  lrun   = ceil((max(rl) - min(bl)) / dl)
  lambda = findgen(lrun) * dl + min(bl)
  bs     = where(lambda le max(bl))
  rs     = where(lambda ge min(rl))

  binner        = interpol(binner       , bl, lambda)
  binterUP      = interpol(binterUP     , bl, lambda)
  binterDN      = interpol(binterDN     , bl, lambda)
  bouterUP      = interpol(bouterUP     , bl, lambda)
  bouterDN      = interpol(bouterDN     , bl, lambda)
;  bouter_opti   = interpol(bouter_opti  , bl, lambda)
  binn_err      = interpol(binn_err     , bl, lambda)
  bint_errUP    = interpol(bint_errUP   , bl, lambda)
  bint_errDN    = interpol(bint_errDN   , bl, lambda)
  bout_errUP    = interpol(bout_errUP   , bl, lambda)
  bout_errDN    = interpol(bout_errDN   , bl, lambda)
;  bout_opti_err = interpol(bout_opti_err, bl, lambda)
  
  rinner        = interpol(rinner       , rl, lambda)
  rinterUP      = interpol(rinterUP     , rl, lambda)
  rinterDN      = interpol(rinterDN     , rl, lambda)
  routerUP      = interpol(routerUP     , rl, lambda)
  routerDN      = interpol(routerDN     , rl, lambda)
;  router_opti   = interpol(router_opti  , rl, lambda)
  rinn_err      = interpol(rinn_err     , rl, lambda)
  rint_errUP    = interpol(rint_errUP   , rl, lambda)
  rint_errDN    = interpol(rint_errDN   , rl, lambda)
  rout_errUP    = interpol(rout_errUP   , rl, lambda)
  rout_errDN    = interpol(rout_errDN   , rl, lambda)
;  rout_opti_err = interpol(rout_opti_err, rl, lambda)
  
  cut = where(lambda gt max(bl))
  binner[cut] = 0
  binter[cut] = 0
  bouter[cut] = 0
  binn_err[cut] = 1d6
  bint_err[cut] = 1d6
  bout_err[cut] = 1d6

  cut = where(lambda lt min(rl))
  rinner[cut] = 0
  rinter[cut] = 0
  router[cut] = 0
  rinn_err[cut] = 1d6
  rint_err[cut] = 1d6
  rout_err[cut] = 1d6

  inner          = dblarr(lrun)
  inner_var      = dblarr(lrun)
  inter          = dblarr(lrun)
  inter_var      = dblarr(lrun)
  outer          = dblarr(lrun)
  outer_var      = dblarr(lrun)
  opti_outer     = dblarr(lrun)
  opti_outer_var = dblarr(lrun)
  for ii = 0, lrun - 1 do begin

     inner[ii] = (binner[ii] / binn_err[ii] + rinner[ii] / rinn_err[ii]) $
                 / (1. / binn_err[ii] + 1. / rinn_err[ii])
     inner_var[ii] = 1. / (1. / binn_err[ii] + 1. / rinn_err[ii])

     inter[ii] = (binter[ii] / bint_err[ii] + rinter[ii] / rint_err[ii]) $
                 / (1. / bint_err[ii] + 1. / rint_err[ii])
     inter_var[ii] = 1. / (1. / bint_err[ii] + 1. / rint_err[ii])

     outer[ii] = (bouter[ii] / bout_err[ii] + router[ii] / rout_err[ii]) $
                 / (1. / bout_err[ii] + 1. / rout_err[ii])
     outer_var[ii] = 1. / (1. / bout_err[ii] + 1. / rout_err[ii])

     opti_outer[ii] = (bouter_opti[ii] / bout_opti_err[ii] + router_opti[ii] / rout_opti_err[ii]) $
                      / (1. / bout_opti_err[ii] + 1. / rout_opti_err[ii])
     opti_outer_var[ii] = 1. / (1. / bout_opti_err[ii] + 1. / rout_opti_err[ii])
     
  endfor
  
  savedata = {LAMBDA: lambda, REDSHIFT: bluData.Z, $
              Z_PHOT: bluData.Z_PHOT, $
              INNER: inner, VAR_INNER: inner_var, $
              INTER: inter, VAR_INTER: inter_var, $
              OUTER: outer, VAR_OUTER: outer_var, $
              OPTI_OUTER: opti_outer, VAR_OPTI_OUTER: opti_outer_var}

  mwrfits, savedata, output, /create
  
end
