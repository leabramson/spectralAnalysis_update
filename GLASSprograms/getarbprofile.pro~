function getProfile, di

  ;; SPECTRAL DIRECTION
  profile = total(di[20:40,20:40], 2)
  s       = size(di, /dim)
  lsf     = mpfitpeak(findgen(n_elements(profile)), profile, /moffat, out)
  tx      = findgen(s[0]) - (out[1] + 20)
  u       = tx / out[2]
  lsf     = out[0] / (u^2 + 1)^out[3]
  lsf /= total(lsf)
  
  ;; SPATIAL DIRECTION
  spatial = total(di[20:40,20:40], 1)
  prof    = mpfitpeak(findgen(n_elements(spatial)), spatial, /moffat, out)
  tx      = findgen(s[1]) - (out[1] + 20)
  u       = tx/out[2]
  profile = out[0] / (u^2 + 1)^out[3]
  tot     = total(profile,/cum) / total(profile)
  q1 = value_locate(tot, 0.25)
  q2 = value_locate(tot, 0.75)
  re = 0.5 * (q2 - q1)  ;;  Half the light's in here

  out[1]  = 0
  dat = {SPATIAL_PARAMS: out, RE: re, LSF: LSF}
  
  RETURN, dat; profile
end
