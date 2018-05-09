pro makederivedtable

  master = loadsummarydata('..')
  ngals = n_elements(master)

  clust  = strarr(ngals)
  id     = strarr(ngals)
  
  mstel  = fltarr(ngals)
  emstel = fltarr(ngals)

  mstel2  = fltarr(ngals) ;; delexps
  emstel2 = fltarr(ngals)
    
  sfr    = fltarr(ngals)
  esfr   = fltarr(ngals)

  sfr2   = fltarr(ngals) ;; delexps
  esfr2  = fltarr(ngals)

  t0     = fltarr(ngals)
  et0    = fltarr(ngals)

  tau    = fltarr(ngals)
  etau   = fltarr(ngals)

  z      = fltarr(ngals)
  ez     = fltarr(ngals)

  av     = fltarr(ngals)
  eav    = fltarr(ngals)

  av2     = fltarr(ngals) ;; delexps
  eav2    = fltarr(ngals)

  re     = fltarr(ngals)
  ere    = fltarr(ngals)
  resex  = fltarr(ngals)
  eresex = fltarr(ngals)
  
  age    = fltarr(ngals)
  eage   = fltarr(ngals)

  tauexp  = fltarr(ngals)
  etauexp = fltarr(ngals)
  
  limit  = bytarr(ngals)

  for ii = 0, ngals - 1 do begin
     clust[ii] = master[ii].FIELD
     id[ii]    = master[ii].ID

     mag  = master[ii].MU
     emag = master[ii].EMU
     lemag = emag / mag / alog(10)
     
     lgndata    = mrdfits(id[ii]+'_lgn_bestfit.fits', 1)

     tinn = where(lgndata.REGION eq 'INNFL')
     tint = where(lgndata.REGION eq 'INTER')
     tout = where(lgndata.REGION eq 'OUTER')

     tme1 = (lgndata[tinn].LMASS[2] - lgndata[tinn].LMASS[1]) / 2.
     tme2 = (lgndata[tint].LMASS[2] - lgndata[tint].LMASS[1]) / 2.
     tme3 = (lgndata[tout].LMASS[2] - lgndata[tout].LMASS[1]) / 2.

     ttot = sqrt(tme1^2 + 2 * tme2^2 + 2 * tme3^2 + lemag^2)

     tsme1 = (lgndata[tinn].SFR[2] - lgndata[tinn].SFR[1]) / 2.
     tsme2 = (lgndata[tint].SFR[2] - lgndata[tint].SFR[1]) / 2.
     tsme3 = (lgndata[tout].SFR[2] - lgndata[tout].SFR[1]) / 2.

     tstot = (sqrt(tsme1^2 + 2 * tsme2^2 + 2 * tsme3^2) / $
              lgndata[where(lgndata.REGION eq 'OPTFL')].SFR[0] $
;     (lgndata[tinn].SFR[0] + 2*lgndata[tint].SFR[0] + 2*lgndata[tout].SFR[0]) $
              / alog(10))^2 + lemag^2
     tstot = sqrt(tstot)

     print, 'lgn', id[ii], ttot, tstot

     lgndata    = lgndata[where(lgndata.REGION eq 'OPTFL')]

     z[ii]      = lgndata.Z[0]
     ez[ii]     = 0.5 * (lgndata.Z[2] - lgndata.Z[1])

     av[ii]     = lgndata.A_V[0]
     eav[ii]    = 0.5 * (lgndata.A_V[2] - lgndata.A_V[1])
     
     tsfr       = lgndata.SFR[0]
     etsfr      = tsfr - lgndata.SFR[1]
     if tsfr lt (2 * etsfr) then begin
        limit[ii] = 1
        sfr[ii]   = alog10(2 * etsfr / mag)
        esfr[ii]  = 0
     endif else begin
        sfr[ii]   = alog10(tsfr / mag)
        etsfr     = 0.5/alog(10) * (lgndata.SFR[2] - lgndata.SFR[1]) / tsfr
        esfr[ii]  = sqrt(etsfr^2 + lemag^2)
     endelse

     mebar = 0.5 * (lgndata.LMASS[2] - lgndata.LMASS[1])
     mstel[ii]  = lgndata.LMASS[0] - alog10(mag)
     emstel[ii] = sqrt(mebar^2 + lemag^2)

     t0[ii]  = alog(lgndata.T0[0] / 1d9)
     et0[ii] = 0.5 * (lgndata.T0[2] - lgndata.T0[1]) / lgndata.T0[0]

     tau[ii] = lgndata.TAU[0]
     etau[ii] = 0.5 * (lgndata.TAU[2] - lgndata.TAU[1])

     imdat = mrdfits('../'+clust[ii]+'/'+id[ii]+'_B.fits', 1)
     foo = mrdfits(imdat.FILENAME, 'sci', head)
     pixscale = sxpar(head, 'CD2_2') ;; asec/pix
     kpc = zang(1.0, z[ii], /silent) ;; asec/kpc
          
     re[ii]     = imdat.RE / sqrt(mag) * pixscale / kpc
     ere[ii]    = (sqrt(1 + emag / mag) - 1) * re[ii]
     resex[ii]  = imdat.RE_SEX / sqrt(mag) * pixscale / kpc
     eresex[ii] = (sqrt(1 + emag / mag) - 1) * resex[ii]
          
     expdata = mrdfits(id[ii]+'_exp_bestfit.fits', 1)

     tinn = where(expdata.REGION eq 'INNFL')
     tint = where(expdata.REGION eq 'INTER')
     tout = where(expdata.REGION eq 'OUTER')

     tme1 = (expdata[tinn].LMASS[2] - expdata[tinn].LMASS[1]) / 2.
     tme2 = (expdata[tint].LMASS[2] - expdata[tint].LMASS[1]) / 2.
     tme3 = (expdata[tout].LMASS[2] - expdata[tout].LMASS[1]) / 2.

     ttot = sqrt(tme1^2 + 2 * tme2^2 + 2 * tme3^2 + lemag^2)

     tsme1 = (expdata[tinn].SFR[2] - expdata[tinn].SFR[1]) / 2.
     tsme2 = (expdata[tint].SFR[2] - expdata[tint].SFR[1]) / 2.
     tsme3 = (expdata[tout].SFR[2] - expdata[tout].SFR[1]) / 2.

     tstot = (sqrt(tsme1^2 + 2 * tsme2^2 + 2 * tsme3^2) / $
              expdata[where(expdata.REGION eq 'OPTFL')].SFR[0] $
;     (expdata[tinn].SFR[0] + 2*expdata[tint].SFR[0] + 2*expdata[tout].SFR[0]) $
              / alog(10))^2 + lemag^2
     tstot = sqrt(tstot)
     
     print, 'exp', id[ii], ttot, tstot
     
     expdata = expdata[where(expdata.REGION eq 'OPTFL')]

     age[ii]  = expdata.AGE[0] / 1d9                         
     eage[ii] = 0.5 * (expdata.AGE[2] - expdata.AGE[1]) / 1d9

     tauexp[ii]  = expdata.TAU[0] / 1d9                         
     etauexp[ii] = 0.5 * (expdata.TAU[2] - expdata.TAU[1]) / 1d9

     av2[ii]     = expdata.A_V[0]
     eav2[ii]    = 0.5 * (expdata.A_V[2] - expdata.A_V[1])
     
     tsfr       = expdata.SFR[0]
     etsfr      = tsfr - expdata.SFR[1]
     if tsfr lt (2 * etsfr) then begin
        limit[ii] = 1
        sfr2[ii]   = alog10(2 * etsfr / mag)
        esfr2[ii]  = 0
     endif else begin
        sfr2[ii]   = alog10(tsfr / mag)
        etsfr     = 0.5/alog(10) * (expdata.SFR[2] - expdata.SFR[1]) / tsfr
        esfr2[ii]  = sqrt(etsfr^2 + lemag^2)
     endelse
     print, sfr2[ii], esfr2[ii]
     
     mebar = 0.5 * (expdata.LMASS[2] - expdata.LMASS[1])
     mstel2[ii]  = expdata.LMASS[0] - alog10(mag)
     emstel2[ii] = sqrt(mebar^2 + lemag^2)
     
  endfor

  close, 1
  openw, 1, 'derivedParameters.tbl'
  for ii = 0, ngals - 1 do begin
     if NOT limit[ii] then begin 
        printf, 1, f = '(%"%s & %s & %5.2f$\\pm$%5.2f & %5.2f$\\pm$%4.2f & %4.2f$\\pm$%4.2f & %4.2f$\\pm$%4.2f & %4.2f$\\pm$%4.2f & $(%4.2f\\pm%4.2f,\\,%4.2f\\pm%4.2f)$ & $(%4.2f\\pm%4.2f,\\,%4.2f\\pm%4.2f)$ ")', $
                clust[ii], id[ii], mstel[ii], emstel[ii], sfr[ii], esfr[ii], re[ii], ere[ii], resex[ii], eresex[ii], av[ii], eav[ii], t0[ii], et0[ii], tau[ii], etau[ii], age[ii], eage[ii], tauexp[ii], etauexp[ii]
     endif else begin
        printf, 1, f = '(%"%s & %s & %5.2f$\\pm$%5.2f & $<$%5.2f & %4.2f$\\pm$%4.2f & %4.2f$\\pm$%4.2f & %4.2f$\\pm$%4.2f & $(%4.2f\\pm%4.2f,\\,%4.2f\\pm%4.2f)$ & $(%4.2f\\pm%4.2f,\\,%4.2f\\pm%4.2f)$ ")', $
                clust[ii], id[ii], mstel[ii], emstel[ii], sfr[ii], re[ii], ere[ii], resex[ii], eresex[ii], av[ii], eav[ii], t0[ii], et0[ii], tau[ii], etau[ii], age[ii], eage[ii], tauexp[ii], etauexp[ii]
     endelse
  endfor
  close, 1

  close, 1
  openw, 1, 'derivedExpParameters.tbl'
  for ii = 0, ngals - 1 do begin
     if NOT limit[ii] then begin 
        printf, 1, f = '(%"%s & %5.2f$\\pm$%5.2f & %5.2f$\\pm$%4.2f & %4.2f$\\pm$%4.2f & $(%4.2f\\pm%4.2f,\\,%4.2f\\pm%4.2f)$ ")', $
                id[ii], mstel2[ii], emstel2[ii], sfr2[ii], esfr2[ii], av2[ii], eav2[ii], age[ii], eage[ii], tauexp[ii], etauexp[ii]
     endif else begin
        printf, 1, f = '(%"%s & %5.2f$\\pm$%5.2f & $<$%5.2f & %4.2f$\\pm$%4.2f & $(%4.2f\\pm%4.2f,\\,%4.2f\\pm%4.2f)$ ")', $
                id[ii], mstel2[ii], emstel2[ii], sfr2[ii], av2[ii], eav2[ii], tau[ii], etau[ii], age[ii], eage[ii], tauexp[ii], etauexp[ii]
     endelse
  endfor
  close, 1
end
