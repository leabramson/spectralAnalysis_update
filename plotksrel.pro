pro plotksrel, lgnlist, explist

  readcol, lgnlist, lgnfiles, f = 'A'  
  readcol, explist, expfiles, f = 'A'

  ngals = n_elements(lgnfiles)

  ms  = fltarr(ngals, 3, 2)
  ems = fltarr(ngals, 3, 2)
  ss  = fltarr(ngals, 3, 2)
  ess = fltarr(ngals, 3, 2)

  for ii = 0, ngals - 1 do begin

     for jj = 0, 1 do begin
        case jj of
           0: begin
              file = lgnfiles[ii]
           end
           1: begin
              file = expfiles[ii]
           end
        endcase
        
        data = mrdfits(file, 1)

        q = value_locate(data.REGION, $
                         ['INNFL', 'INTER', 'OUTER'])
        print, data[q].REGION

        ms[ii,*,jj]  = data[q].LMASS[0] - alog10(data[q].AREA_PHYS)
        ems[ii,*,jj] = sqrt(0.5 * (data[q].LMASS[2] - data[q].LMASS[1])^2 + 0.065^2)

        lsfr = data[q].SFR[0] - data[q].SFR[1]
        nsig = 2.
        qui = where(data[q].SFR[0] gt nsig * lsfr)
        uls = where(data[q].SFR[0] le nsig * lsfr)
        ss[ii,qui,jj] = data[q[qui]].SFR[0] / data[q[qui]].AREA_PHYS
        ess[ii,qui,jj] = 0.5 * (data[q[qui]].SFR[2] - data[q[qui]].SFR[1]) / data[q[qui]].AREA_PHYS
        if min(uls) ge 0 then $
           ess[ii,uls,jj] = nsig * lsfr[uls] / data[q[uls]].AREA_PHYS

     endfor
  endfor

  ;; fit relation (LEA 20170822)
  niter = 200
  m = fltarr(2,niter)
  b = fltarr(2,niter)
  s = fltarr(2,niter)
  pivot = fltarr(2)
  for jj = 0, 1 do begin
     tms = reform(ms[*,q,jj])
     tems = reform(ems[*,q,jj])
     tss = reform(ss[*,q,jj])
     tes = reform(ess[*,q,jj])
     use = where(alog10(tss) gt -2, nuse)
     tserr = tes[use]/tss[use]/alog(10)
     pivot[jj] = median(tms[use])
     for ii = 0, niter - 1 do begin
        xx = tms[use] + randomn(seed, nuse) * tems[use] - pivot[jj]
        yy = alog10(tss[use]) + randomn(seed, nuse) * tes[use]
        foo = poly_fit(xx, yy, 1)
        b[jj,ii] = foo[0]
        m[jj,ii] = foo[1]
        s[jj,ii] = sqrt(variance(alog10(tss[use]) - (foo[0] + foo[1] * (tms[use] - pivot[jj])), /nan) $
                        - mean(tes[use]^2, /nan))
     endfor
  endfor
  m0 = 7.6
  m1 = 9.9
  dm = 0.1
  tmasses = findgen((m1-m0)/dm+1) * dm + m0
  fit1 = fltarr(n_elements(tmasses), niter)
  for ii = 0, niter - 1 do $
     fit1[*,ii] = b[0,ii] + m[0,ii] * (tmasses - pivot[0])
  tfit1 = fltarr(n_elements(tmasses),2)
  tfit1[*,0] = median(fit1, dim = 2)
  tfit1[*,1] = stddev(fit1, dim = 2, /nan)
  
  fit1 = fltarr(n_elements(tmasses), niter)
  for ii = 0, niter - 1 do $
     fit1[*,ii] = b[1,ii] + m[1,ii] * (tmasses - pivot[1])
  tfit2 = fltarr(n_elements(tmasses),2)
  tfit2[*,0] = median(fit1, dim = 2)
  tfit2[*,1] = stddev(fit1, dim = 2, /nan)
  
;  stop

;  cd16rel = -7.63 + 0.68 * tmasses
  
  cols = [255,'00a500'x,'ff5500'x]
  set_plot, 'PS'
  device, filename = 'ksrel.eps', $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = 7, ysize = 4, /in
  plot, ms[*,*,0], alog10(ss[*,*,0]), psym = 1, $
        /nodat, $
        xtitle = 'log '+greek('Sigma')+'!D!18M!X!L*!N [M!D'+sunsymbol()+'!N kpc!E-2!N]', $
        ytitle = 'log '+greek('Sigma')+'!D!18SFR!X!N [M!D'+sunsymbol()+'!N yr!E-1!N kpc!E-2!N]', $
        xran = [7.5,10], yran = [-2,1], /ysty, $
        pos = [0.15,0.15,0.5,0.9], $
        xthick = 3, ythick = 3, charsize = 1.25, charthick = 4
  polyfill, [tmasses, reverse(tmasses)], $
            [tfit1[*,0]-tfit1[*,1], reverse(tfit1[*,0]+tfit1[*,1])], $
            col = 'cccccc'x
;  oplot, tmasses, tfit1[*,0], thick = 6
  oplot, tmasses, median(b[0,*]) + median(m[0,*]) * (tmasses - pivot[0]), thick = 6
;  oplot, tmasses, cd16rel, linesty = 4, thick = 2
;  for ii = 0, niter - 1 do $
;     oplot, !X.CRANGE, b[0,ii] + m[0,ii] * (!X.CRANGE - pivot[0]), col = 'cccccc'x, thick = 1
  for ii = 0, 2 do begin
     q = where(ss[*,ii,0] gt 0, compl = uls)
     plotsym, 0, /fill
     oploterror, [ms[q,ii,0]], alog10([ss[q,ii,0]]), $
                 [ems[q,ii,0]], [ess[q,ii,0]/ss[q,ii,0]/alog(10)], $
                 psym = 8, symsize = 1.25, /nohat, errthick = 4
     oplot, [ms[q,ii,0]], alog10([ss[q,ii,0]]), psym = 8, col = cols[ii]
     plotsym, 1, thick = 4
     if min(uls) ge 0 then begin
        for jj = 0, n_elements(uls) -1 do begin
           oplot, ms[uls[jj],ii,0] + ems[uls[jj],ii,0] * [-1,1], $
                  replicate(alog10(ess[uls[jj],ii,0])> (!Y.CRANGE[0] + 0.2), 2), $
                  col = cols[ii], thick = 4
           oplot, [ms[uls[jj],ii,0]], alog10([ess[uls[jj],ii,0]]) > [!Y.CRANGE[0] + 0.2], $
                  psym = 8, col = cols[ii], symsize = 2
        endfor
     endif
  endfor
  cgtext, !X.CRANGE[0]+0.1, !Y.CRANGE[1]-0.3, /data, $
          '!18y!X=('+string(median(m[0,*]), f = '(F4.2)')+textoidl('\pm')+$
          string(stddev(m[0,*], /nan), f = '(F4.2)')+')!18x!X-('+$
          string(abs(median(b[0,*])), f = '(F4.2)')+textoidl('\pm')+$
          string(stddev(b[0,*], /nan), f = '(F4.2)')+')', charsize = 0.95, charthick = 3
  cgtext, !X.CRANGE[0]+0.1, !Y.CRANGE[1]-0.5, /data, $
          greek('sigma')+'!Dint!N='+string(median(s[0,*]), f = '(F4.2)')+' dex', $
          charsize = 1, charthick = 4
  plotsym, 0, /fill
  legend, pos = [!X.CRANGE[0]+0.1, !Y.CRANGE[1]-0.65], /data, box = 1, $
          ['Inner', 'Middle', 'Outer'], $
          psym = 8, col = cols, symsize = [1,1,1], pspacing = 0.5, $
          charsize = 1.1, charthick = 3, spacing = 1.5
          
  plot, ms[*,*,0], alog10(ss[*,*,0]), psym = 1, $
        /nodat, $
        xtitle = 'log '+greek('Sigma')+'!D!18M!X!L*!N [M!D'+sunsymbol()+'!N kpc!E-2!N]', $
;        ytitle = 'log '+greek('Sigma')+'!D!18SFR!X!N [M!D'+sunsymbol()+'!N yr!E-1!N kpc!E-2!N]', $
        xran = [7.5+1d-6,10], yran = [-2,1], /xsty, $
        pos = [0.5,0.15,0.85,0.9], ysty = 8+1, ytickname = replicate(' ',60), /noer, $
        xthick = 3, ythick = 3, charsize = 1.25, charthick = 4
  axis, yaxis = 1, yran = !Y.CRANGE, /ysty, $
        ytitle = 'log '+greek('Sigma')+'!D!18SFR!X!N [M!D'+sunsymbol()+'!N yr!E-1!N kpc!E-2!N]', $
        ythick = 3, charsize = 1.25, charthick = 4
  polyfill, [tmasses, reverse(tmasses)], $
            [tfit2[*,0]-tfit2[*,1], reverse(tfit2[*,0]+tfit2[*,1])], $
            col = 'cccccc'x
;  oplot, tmasses, tfit2[*,0], thick = 6
  oplot, tmasses, median(b[1,*]) + median(m[1,*]) * (tmasses - pivot[1]), thick = 6
  for ii = 0, 2 do begin
     q = where(ss[*,ii,1] gt 0, compl = uls)
     plotsym, 0, /fill
     oploterror, [ms[q,ii,1]], alog10([ss[q,ii,1]]), $
                 [ems[q,ii,1]], [ess[q,ii,1]/ss[q,ii,1]/alog(10)], $
                 psym = 8, symsize = 1.25, /nohat, errthick = 4, $
                 col = '777777'x, errcol = '777777'x
     oplot, [ms[q,ii,1]], alog10([ss[q,ii,1]]), psym = 8, col = cols[ii]
     plotsym, 1, thick = 4
     if min(uls) ge 0 then begin
        for jj = 0, N_elements(uls) - 1 do begin
           oplot, ms[uls[jj],ii,1] + ems[uls[jj],ii,1] * [-1,1], $
                  replicate(alog10([ess[uls[jj],ii,1]]) > [!Y.CRANGE[0] + 0.2], 2), $
                  col = cols[ii], thick = 4
           oplot, [ms[uls[jj],ii,1]], alog10([ess[uls[jj],ii,1]]) > [!Y.CRANGE[0] + 0.2], $
                  psym = 8, col = cols[ii], symsize = 2
        endfor
     endif
  endfor
  cgtext, !X.CRANGE[0]+0.1, !Y.CRANGE[1]-0.3, /data, $
          '!18y!X=('+string(median(m[1,*]), f = '(F4.2)')+textoidl('\pm')+$
          string(stddev(m[1,*], /nan), f = '(F4.2)')+')!18x!X-('+$
          string(abs(median(b[1,*])), f = '(F4.2)')+textoidl('\pm')+$
          string(stddev(b[1,*], /nan), f = '(F4.2)')+')', charsize = 0.95, charthick = 3
  cgtext, !X.CRANGE[0]+0.1, !Y.CRANGE[1]-0.5, /data, $
          greek('sigma')+'!Dint!N='+string(median(s[1,*]), f = '(F4.2)')+' dex', $
          charsize = 1, charthick = 4
  device, /close
  spawn, 'gv ksrel.eps &'
  set_plot, 'X'
  print, pivot
  
  stop
  
end
;plotksrel, 'study5_lgn_Bestfits.list', 'study5_exp_Bestfits.list'
