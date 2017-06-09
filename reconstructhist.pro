pro reconstructhist, dir, region, $
                     SFH = sfh, $
                     NOTMASKED = notmasked, $
                     OUTPUT = output

  if NOT keyword_set(SFH) then $
     sfh = 'DelayedExponential'

  if NOT keyword_set(NOTMASKED) then $
     notmasked = 0 $
  else $
     notmasked = 1
  
  prefix = dir+'/'+region+'.masked'
  if notmasked then prefix = repstr(prefix, '.masked', '')

  ;; Get the "master" info
  master = prefix+'.analyze'
  spawn, 'cat '+master+' | grep logmstar > tmpMass.dat'
  spawn, 'cat '+master+' | grep sfr > tmpSFR.dat'
  spawn, 'cat '+master+' | grep z > tmpZ.dat'
  
  readcol, 'tmpMass.dat', $
           lmass, lmasslo, lmasshi, f = 'X,F,F,F', /silent 
  readcol, 'tmpSFR.dat', $
           sfr, sfrlo, sfrhi, f = 'X,F,F,F', /silent 
  if file_lines('tmpZ.dat') gt 0 then begin
     readcol, 'tmpZ.dat', $
              z, zlo, zhi, f = 'X,F,F,F', /silent
     zcheck = 1
  endif else begin
     if sfh eq 'LogNormal' then $
        tdir = repstr(dir, 'LogNormal', 'Phot') $
     else $
        tdir = dir
     if NOT notmasked then $
        spawn, 'cat '+tdir+'/INNFL.masked.analyze | grep z > tmpZ.dat' $
     else $
        spawn, 'cat '+tdir+'/INNFL.analyze | grep z > tmpZ.dat'
     readcol, 'tmpZ.dat', $
              z, zlo, zhi, f = 'X,F,F,F', /silent
     zcheck = 0
  endelse
  z_obs = z[0]
  
  mass = 10.^lmass[0]
  masshi = 10.^(lmasshi[0] - lmass[0]) * mass
  masslo = 10.^(lmass[0] - lmasslo[0]) * mass

  sfr = sfr[1]
  sfrhi = sfrhi[1]
  sfrlo = sfrlo[1]

  ;; Set clocks
  t0 = getage(20)
  t1 = getage(0)
  dt = 0.05
  time = findgen((t1-t0)/dt+1)*dt + t0
  tobs = getage(z_obs)
  obsdex = value_locate(time, tobs)
  
  ;; Read the chains. This is tricky because "z" might not always be a
  ;; parameter...
  chainFile = prefix+'post_equal_weights.dat'
  case sfh of
     'DelayedExponential': begin
        ;; SII6718_ew z SII6733_ew OIII_ew
        ;; H_ew A_V age ltau
        ;; logmstar logssfr sfr
        ;; Ha_flux Hb_flux OIII_flux
        ;; SII6718_flux SII6733_flux
        if zcheck then $
           readcol, chainFile, $
                    siiew1, z, siiew2, oiiiew, $
                    hew, av, age, ltau, $
                    chispec, chiphot, lmass, $
                    lssfr, sfr, $
                    haf, hbf, oiiif, $
                    siif1, siif2, /quick, /silent $
        else $
           readcol, chainFile, $
                    siiew1, hew, oiiiew, $
                    siiew2, av, age, ltau, $
                    chispec, chiphot, lmass, $
                    lssfr, sfr, $
                    haf, hbf, oiiif, $
                    siif1, siif2, /quick, /silent

        ndata = n_elements(age)
        massHist = fltarr(n_elements(time), ndata)
        sfrHist  = fltarr(n_elements(time), ndata)
        ssfrHist = fltarr(n_elements(time), ndata)

        mass = 10.^lmass
        tau  = 10.^ltau/1d9
        age /= 1d9

        for ii = 0, ndata - 1 do begin
           t = time - (tobs - age[ii])
           startdex = value_locate(time, tobs - age[ii]) > 0
           if startdex eq obsdex then startdex -= 1
           sfh = t/tau[ii]^2 * exp(-t/tau[ii])

           ;; Normalize the SFH to the inferred SFR
           sfrHist[*,ii] = sfh / sfh[obsdex] * sfr[ii] 

           ;; Tabulate the mass growth history assuming the modern retention factor is historically accurate
           norm = mass[ii] / (tsum(time[startdex:obsdex], sfh[startdex:obsdex]) * 1d9)
           masshist[0:startdex-1,ii] = 0
           masshist[startdex:*,ii] = total(sfh[startdex:*] * norm * dt * 1d9, /cum) 

           ssfrHist[*,ii] = sfrHist[*,ii] / masshist[*,ii]
           
        endfor
        
     end
     'LogNormal': begin
        prefix = repstr(prefix, 'Phot', 'LogNormal')
        chainFile = prefix+'post_equal_weights.dat'
        readcol, chainFile, $
                 siiew1, hew, oiiiew, $
                 hbew, siiew2, tau, t0, av, $
                 chispec, chiphot, lmass, $
                 lssfr, sfr, $
                 oiiif, siif1, siif2, haf, hbf, $
                 /quick, /silent

        ndata = n_elements(t0)
        massHist = fltarr(n_elements(time), ndata)
        sfrHist  = fltarr(n_elements(time), ndata)
        ssfrHist = fltarr(n_elements(time), ndata)
        
        mass = 10.^lmass
        t0 /= 1d9
        
        for ii = 0, ndata - 1 do begin
           sfh = exp(-0.5 * alog(time/t0[ii])^2/tau[ii]^2) / time

           ;; Normalize the SFH to the inferred SFR
           sfrHist[*,ii] = sfh / sfh[obsdex] * sfr[ii] 

           ;; Tabulate the mass growth history assuming the modern retention factor is historically accurate
           norm = mass[ii] / (tsum(time[0:obsdex], sfh[0:obsdex]) * 1d9)
           masshist[*,ii] = total(sfh * norm * dt * 1d9, /cum) 

           ssfrHist[*,ii] = sfrHist[*,ii] / masshist[*,ii]
           
        endfor
        
     end
  endcase

  mgh  = fltarr(n_elements(time), 3)
  sfh  = fltarr(n_elements(time), 3)
  ssfh = fltarr(n_elements(time), 3)
  for ii = 0, n_elements(time) - 1 do begin
     tmp = masshist[ii,sort(masshist[ii,*])]
     mgh[ii,*] = [tmp[0.16 * ndata], $
                  tmp[0.50 * ndata], $
                  tmp[0.84 * ndata]]
     tmp = sfrHist[ii,sort(sfrHist[ii,*])]
     sfh[ii,*] = [tmp[0.16 * ndata], $
                  tmp[0.50 * ndata], $
                  tmp[0.84 * ndata]]
     tmp = ssfrHist[ii,sort(ssfrHist[ii,*])]
     ssfh[ii,*] = [tmp[0.16 * ndata], $
                   tmp[0.50 * ndata], $
                   tmp[0.84 * ndata]]
  endfor

  com = $
     'savedata = {REGION: region, '+$
     'TIME    : time, '+$
     'REDSHIFT: getredshift(time), '+$
     'TOBS    : tobs, '+$
     'SFH     : sfh, '+$
     'MGH     : mgh, '+$
     'SSFH    : ssfh}'
  void = execute(com)

  mwrfits, savedata, output, /create
  
end

;;
;;
;;

pro test1, dir

;  dir = 'MACS1149/00900_1_pyspecfitPhotResults'
;  dir = 'ABEL0370/03762_1_pyspecfitPhotResults'

  regions = ['INNFL', 'INTER', 'OUTER', $
             'INTUP', 'INTDN', 'OUTUP', 'OUTDN', $
             'OPTFL']
  cols    = [255, '00aa00'x, 'ff5500'x, $
             '00ff00'x, '00aa00'x, 'ffa500'x, 'ff0000'x, $
             '777777'x]

  reconstructhist, dir, regions[0], $
                   output = 'test_'+regions[0]+'_reconstruct.fits'
  tdat = mrdfits('test_'+regions[0]+'_reconstruct.fits', 1)

  sfhs = fltarr(n_elements(tdat.TIME), n_elements(regions), 3)
  mghs = fltarr(n_elements(tdat.TIME), n_elements(regions), 3)

  sfhs[*,0,*] = tdat.SFH
  mghs[*,0,*] = tdat.MGH
  
  for ii = 1, n_elements(regions) - 1 do begin
     reconstructhist, dir, regions[ii], $
                      output = 'test_'+regions[ii]+'_reconstruct.fits', $
                      sfh = 'LogNormal'
     tdat = mrdfits('test_'+regions[ii]+'_reconstruct.fits', 1)
     sfhs[*,ii,*] = tdat.SFH
     mghs[*,ii,*] = tdat.MGH
  endfor

  time = tdat.TIME
  tobs = tdat.TOBS

  xxx = [time, reverse(time)]  

  !p.multi = [0,2,0]
  
  plot, time, alog10(sfhs[*,-1,1]), /nodat, $
        yran = [0,3], $
        xtitle = '!18t!X [Gyr]', $
        ytitle = 'log !18SFR!X [M!D'+sunsymbol()+'!N yr!E-1!N]', $
        xran = [-2,2] + tobs, /xsty
  for ii = 0, n_elements(regions) - 2 do begin
     yyy = alog10([sfhs[*,ii,0], reverse(sfhs[*,ii,2])])
     polyfill, xxx, yyy > !Y.CRANGE[0], col = cols[ii], $
               /line_fill, thick = 1, spacing = 0.05, orien = ii * 10
     oplot, time, alog10(sfhs[*,ii,1]), thick = 4
     oplot, time, alog10(sfhs[*,ii,1]), col = cols[ii]
  endfor
  oplot, replicate(tobs, 2), !Y.CRANGE, linesty = 5

  plot, time, alog10(mghs[*,-1,1]), /nodat, $
        yran = [-1, 0.7] + alog10(mghs[value_locate(time, tobs),-2,1])  , $
        xtitle = '!18t!X [Gyr]', $
        ytitle = 'log !18M!X!D*!N(!18t!X) [M!D'+sunsymbol()+'!N]', $
        xran = [-2,2] + tobs, /xsty, /ysty
  for ii = 0, n_elements(regions) - 2 do begin
     yyy = alog10([mghs[*,ii,0], reverse(mghs[*,ii,2])])
     polyfill, xxx, yyy > !Y.CRANGE[0], col = cols[ii], $
               /line_fill, thick = 1, spacing = 0.05, orien = ii * 10
     oplot, time, alog10(mghs[*,ii,1]), thick = 4
     oplot, time, alog10(mghs[*,ii,1]), col = cols[ii]
  endfor
  oplot, replicate(tobs, 2), !Y.CRANGE, linesty = 5

  !p.multi = 0  
    
end

;;
;;
;;

pro testevo

  innfl = mrdfits('test_INNFL_reconstruct.fits', 1)
  intup = mrdfits('test_INTUP_reconstruct.fits', 1)
  intdn = mrdfits('test_INTDN_reconstruct.fits', 1)
  outup = mrdfits('test_OUTUP_reconstruct.fits', 1)
  outdn = mrdfits('test_OUTDN_reconstruct.fits', 1)

  optfl = mrdfits('test_OPTFL_reconstruct.fits', 1)
  
  tobs = innfl.TOBS

  !p.multi = [0,2,0]

  plotsym, 8, /fill
  plot, [8,11], [-11,-7], /nodat, $
        xtitle = 'log M(t)', ytitle = 'log sSFR(t)', /iso
  epochs = where(optfl.mgh ge 1d9 AND optfl.time le tobs+1)
  xxx = alog10([optfl.MGH[epochs,0], reverse(optfl.MGH[epochs,0])])
  yyy = alog10([optfl.SFH[epochs,1]/optfl.MGH[epochs,0], $
                reverse(optfl.SFH[epochs,2]/optfl.MGH[epochs,0])])
  polyfill, xxx, yyy, col = '777777'x

  regions = ['INNFL', 'INTUP', 'INTDN', 'OUTUP', 'OUTDN', 'OPTFL']
  cols    = [255, '00ff00'x, '005500'x, 'ffa500'x, 'ff0000'x, '777777'x]

  plot, [8,11], [-11,-7], /nodat, $
        xtitle = 'log M(t)', ytitle = 'log sSFR(t)', /iso
  for ii = 0, 4 do begin
     case ii of
        0: begin
           xxx = alog10([innfl.MGH[epochs,0], reverse(innfl.MGH[epochs,0])]) > !x.crange[0]
           yyy = alog10([innfl.SFH[epochs,1]/innfl.MGH[epochs,0], $
                         reverse(innfl.SFH[epochs,2]/innfl.MGH[epochs,0])]) > !Y.crange[0]
        end
        1: begin
           xxx = alog10([intup.MGH[epochs,0], reverse(intup.MGH[epochs,0])]) > !x.crange[0]
           yyy = alog10([intup.SFH[epochs,1]/intup.MGH[epochs,0], $
                         reverse(intup.SFH[epochs,2]/intup.MGH[epochs,0])]) > !Y.crange[0]
         end
        2: begin
           xxx = alog10([intdn.MGH[epochs,0], reverse(intdn.MGH[epochs,0])]) > !x.crange[0]
           yyy = alog10([intdn.SFH[epochs,1]/intdn.MGH[epochs,0], $
                         reverse(intdn.SFH[epochs,2]/intdn.MGH[epochs,0])]) > !Y.crange[0]
        end
        3: begin
           xxx = alog10([outup.MGH[epochs,0], reverse(outup.MGH[epochs,0])]) > !x.crange[0]
           yyy = alog10([outup.SFH[epochs,1]/outup.MGH[epochs,0], $
                         reverse(outup.SFH[epochs,2]/outup.MGH[epochs,0])]) > !Y.crange[0]
        end
        4: begin
           xxx = alog10([outdn.MGH[epochs,0], reverse(outdn.MGH[epochs,0])]) > !x.crange[0]
           yyy = alog10([outdn.SFH[epochs,1]/outdn.MGH[epochs,0], $
                         reverse(outdn.SFH[epochs,2]/outdn.MGH[epochs,0])]) > !Y.crange[0]
        end
     endcase
     polyfill, xxx, yyy, col = cols[ii]
  endfor

  !p.multi = 0
  
  time = optfl.TIME
  niter = n_elements(time)
  distm = fltarr(niter, 5)
  dists = fltarr(niter, 5)
  for ii = 0, niter - 1 do begin
     dm = alog10([innfl.MGH[ii,0], intUP.MGH[ii,0], intDn.MGH[ii,0], outUp.MGH[ii,0], outDN.MGH[ii,0]]) - alog10(optfl.MGH[ii,0])
     ds = alog10([innfl.SFH[ii,0], intUP.SFH[ii,0], intDn.SFH[ii,0], outUp.SFH[ii,0], outDN.SFH[ii,0]]) - alog10(optfl.SFH[ii,0])
     distm[ii,*] = dm
     dists[ii,*] = ds
  endfor

  ms = fltarr(niter, 5)
  ss = fltarr(niter, 5)
  for ii = 0, niter - 1 do begin
     ms[ii,*] = alog10([innfl.MGH[ii,0], intUP.MGH[ii,0], intDn.MGH[ii,0], outUp.MGH[ii,0], outDN.MGH[ii,0]])
     ss[ii,*] = alog10([innfl.SFH[ii,0], intUP.SFH[ii,0], intDn.SFH[ii,0], outUp.SFH[ii,0], outDN.SFH[ii,0]])
  endfor

  
  
  stop
  
end
