pro readtakaphotoz, incat, output

  ;id z_spec z_a z_m1 chi_a z_p chi_p z_m2 odds l68 u68  l95 u95  l99 u99  nfilt q_z z_peak peak_prob z_mc

  readcol, incat    , $
           id       , $
           z_spec   , $
           z_a      , $
           z_m1     , $
           chi_a    , $
           z_p      , $
           chi_p    , $
           z_m2     , $
           odds     , $
           l68      , $
           u68      , $
           l95      , $
           u95      , $
           l99      , $
           u99      , $
           nfilt    , $
           q_z      , $
           z_peak   , $
           peak_prob, $
           z_mc     , $
           f = $
           'L,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,',$
           /silent
  nsrcs = n_elements(id)
  
  savedata = {id_taka  : 0L, $
              z_spec   : 0., $
              z_a      : 0., $
              z_m1     : 0., $
              chi_a    : 0., $
              z_p      : 0., $
              chi_p    : 0., $
              z_m2     : 0., $
              odds     : 0., $
              l68      : 0., $
              u68      : 0., $
              l95      : 0., $
              u95      : 0., $
              l99      : 0., $
              u99      : 0., $
              nfilt    : 0., $
              q_z      : 0., $
              z_peak   : 0., $
              peak_prob: 0., $
              z_mc     : 0.}
  savedata =replicate(savedata, nsrcs)

  for ii = 0, nsrcs - 1 do begin
     savedata[ii].id_taka   = id[ii]     
     savedata[ii].z_spec    = z_spec[ii]       
     savedata[ii].z_a       = z_a[ii]      
     savedata[ii].z_m1      = z_m1[ii]     
     savedata[ii].chi_a     = chi_a[ii]     
     savedata[ii].z_p       = z_p[ii]     
     savedata[ii].chi_p     = chi_p[ii]     
     savedata[ii].z_m2      = z_m2[ii]     
     savedata[ii].odds      = odds[ii]     
     savedata[ii].l68       = l68[ii]     
     savedata[ii].u68       = u68[ii]     
     savedata[ii].l95       = l95[ii]     
     savedata[ii].u95       = u95[ii]     
     savedata[ii].l99       = l99[ii]     
     savedata[ii].u99       = u99[ii]     
     savedata[ii].nfilt     = nfilt[ii]     
     savedata[ii].q_z       = q_z[ii]     
     savedata[ii].z_peak    = z_peak[ii]     
     savedata[ii].peak_prob = peak_prob[ii]     
     savedata[ii].z_mc      = z_mc[ii]     
  endfor

  mwrfits, savedata, output, /create

end
