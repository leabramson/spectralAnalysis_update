pro plotsedhasfrs, lgnlist, explist

  summaryData = mrdfits('objectSummaryData.fits', 1)
  
  readcol, lgnlist, lgnfiles, f = 'A'  
  readcol, explist, expfiles, f = 'A'

  ngals = n_elements(lgnfiles)

  regions = ['OPTFL', 'INNFL', 'INTER', 'OUTER']
  
  ms  = fltarr(ngals, 4, 2)
  ems = fltarr(ngals, 4, 2)
  ss  = fltarr(ngals, 4, 2)
  ess = fltarr(ngals, 4, 2)

  mags  = []
  emags = []
  z     = []
  d     = []
  for ii = 0, ngals - 1 do begin

     for jj = 0, 1 do begin
        case jj of
           0: begin
              file = lgnfiles[ii]
              nn = strmid(file, 0, 7)
;              mag = summaryData[where(summaryData.ID eq nn)].MU
;              emag = summaryData[where(summaryData.ID eq nn)].EMU / mag / alog(10)
              mags  = [mags, summaryData[where(summaryData.ID eq nn)].MU]
              emags  = [emags, summaryData[where(summaryData.ID eq nn)].EMU]
              z = [z, summaryData[where(summaryData.ID eq nn)].Z]
              d = [d, dluminosity(z, /cm)]
           end
           1: begin
              file = expfiles[ii]
           end
        endcase
        
        data = mrdfits(file, 1)
        
        for kk = 0, n_elements(regions) - 1 do begin

           td = data[where(data.REGION eq regions[kk])]

           ms[ii,kk,jj] = median(td.SFR / mags[ii])
           ems[ii,kk,jj] = stddev(td.SFR / mags[ii], /nan); / median(td.SFR / mags[ii]) / alog(10)

           if ii eq 0 then stop
           
           lsfrs = gethasfr(td.HA_FLUX, z[ii], hbf = td.HB_FLUX, /dustcor) * 2./3.
           
;           line = td.HA_FLUX / 1d18 * 0.67 / mags[ii] ;; Shapley+15 assume it's 1/3 NII
;           line_corr = line * 10.^(0.4 * td.AV * (1.9 - 0.15 * td.AV)) ;; Mason+17 from Wuyts+13
;           lsfrs = 4 * !pi * line_corr * d^2 / 1.26d41 ;; K98

           if jj eq 0 then begin
              ss[ii,kk,*]  = median(lsfrs / mags[ii])
              ess[ii,kk,*] = stddev(lsfrs / mags[ii], /nan) ; / median(lsfrs / mags[ii]) / alog(10)
           endif
        endfor
                   
     endfor
  endfor

  emags = emags / mags / alog(10)
  
;  plot, alog10(ms[*,*,0]), alog10(ss[*,*,0]), $
;        psym = 1, xran = [-0.5,2], yran = [-0.5,2], /iso, /nodat
;  oploterror, alog10(ms[*,*,0]), alog10(ss[*,*]), $
;              sqrt(ems[*,*,0]^2 + emags^2), sqrt(ess[*,*]^2 + emags^2), /nohat, $
;              col = 'ffa500'x, psym = 1, errcol = 'ffa500'x
;  oploterror, alog10(ms[*,*,1]), alog10(ss[*,*]), $
;              sqrt(ems[*,*,1]^2 + emags^2), sqrt(ess[*,*]^2 + emags^2), /nohat, $
;              col = '00a5ff'x, psym = 1, errcol = '00a5ff'x
;;  oplot, alog10(ms[*,*,1]), alog10(ss[*,*,1]), psym = 1, col = '00a5ff'x                  
;  one_one
  
;  stop

  ms = ms[0,*,*]
  ems = ems[0,*,*]
  ss = ss[0,*,*]
  ess = ess[0,*,*]
  
  nsig = 1

  cols = ['555555'x, 255,'00a500'x,'ff5500'x]
  set_plot, 'PS'
  device, filename = 'sedHaSFR.eps', $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = 7, ysize = 4, /in
  plot, alog10(ms[*,*,0]), alog10(ss[*,*,0]), psym = 1, $
        /nodat, $
        xtitle = 'log !18SFR!X!Dlgnml!N [M!D'+sunsymbol()+'!N yr!E-1!N]', $
        ytitle = 'log !18SFR!X!DH'+greek('alpha')+'!N [M!D'+sunsymbol()+'!N yr!E-1!N]', $
        xran = [0,2], yran = [0,2], /xsty, /ysty, $
        pos = [0.15,0.15,0.5,0.9], $
        xthick = 3, ythick = 3, charsize = 1.25, charthick = 4
  oplot, !X.CRANGE, !Y.CRANGE, linesty = 1
  for ii = 0, n_elements(regions) - 1 do begin
     q = where(ss[*,ii,0] gt nsig * ess[*,ii,0] $
               AND $
               ms[*,ii,0] gt nsig * ems[*,ii,0], compl = uls)
     plotsym, 0, /fill
     oploterror, alog10([ms[q,ii,0]]), alog10([ss[q,ii,0]]), $
                 [ems[q,ii,0]/ms[q,ii,0]/alog(10)], [ess[q,ii,0]/ss[q,ii,0]/alog(10)], $
                 psym = 8, symsize = 1.25, /nohat, errthick = 4
     oplot, alog10([ms[q,ii,0]]), alog10([ss[q,ii,0]]), psym = 8, col = cols[ii]
     plotsym, 1, thick = 4
     if min(uls) ge 0 then begin
        for jj = 0, n_elements(uls) -1 do begin
           oplot, alog10(ms[uls[jj],ii,0]) + ems[uls[jj],ii,0]/ms[uls[jj],ii,0]/alog(10) * [-1,1], $
                  replicate(alog10(nsig * ess[uls[jj],ii,0])> (!Y.CRANGE[0] + 0.2), 2), $
                  col = cols[ii], thick = 4
           oplot, alog10([ms[uls[jj],ii,0]]), alog10(nsig * [ess[uls[jj],ii,0]]) > [!Y.CRANGE[0] + 0.2], $
                  psym = 8, col = cols[ii], symsize = 2
        endfor
     endif
  endfor
;  cgtext, !X.CRANGE[0]+0.1, !Y.CRANGE[1]-0.3, /data, $
;          '!18y!X=('+string(median(m[0,*]), f = '(F4.2)')+textoidl('\pm')+$
;          string(stddev(m[0,*], /nan), f = '(F4.2)')+')!18x!X-('+$
;          string(abs(median(b[0,*])), f = '(F4.2)')+textoidl('\pm')+$
;          string(stddev(b[0,*], /nan), f = '(F4.2)')+')', charsize = 0.95, charthick = 3
;  cgtext, !X.CRANGE[0]+0.1, !Y.CRANGE[1]-0.5, /data, $
;          greek('sigma')+'!Dint!N='+string(median(s[0,*]), f = '(F4.2)')+' dex', $
;          charsize = 1, charthick = 4
  plotsym, 0, /fill
  legend, /bottom, /right, box = 1, $;pos = [!X.CRANGE[0]+0.1, !Y.CRANGE[1]-0.65], /data, box = 1, $
          ['Optimal', 'Inner', 'Middle', 'Outer'], $
          psym = 8, col = cols, symsize = [1,1,1,1], pspacing = 0.5, $
          charsize = 1.1, charthick = 3, spacing = 1.5
          
  plot, alog10(ms[*,*,1]), alog10(ss[*,*,1]), psym = 1, $
        /nodat, $
        xtitle = 'log !18SFR!X!DDelExp!N [M!D'+sunsymbol()+'!N yr!E-1!N]', $
;        ytitle = 'log '+greek('Sigma')+'!D!18SFR!X!N [M!D'+sunsymbol()+'!N yr!E-1!N kpc!E-2!N]', $
        xran = [1d-6,2], yran = [0,2], /xsty, $
        pos = [0.5,0.15,0.85,0.9], ysty = 8+1, ytickname = replicate(' ',60), /noer, $
        xthick = 3, ythick = 3, charsize = 1.25, charthick = 4
  axis, yaxis = 1, yran = !Y.CRANGE, /ysty, $
        ytitle = 'log !18SFR!X!DH'+greek('alpha')+'!N [M!D'+sunsymbol()+'!N yr!E-1!N]', $
        ythick = 3, charsize = 1.25, charthick = 4
  oplot, !X.CRANGE, !Y.CRANGE, linesty = 1
;  polyfill, [tmasses, reverse(tmasses)], $
;            [tfit2[*,0]-tfit2[*,1], reverse(tfit2[*,0]+tfit2[*,1])], $
;            col = 'cccccc'x
;  oplot, tmasses, tfit2[*,0], thick = 6
;  oplot, tmasses, median(b[1,*]) + median(m[1,*]) * (tmasses - pivot[1]), thick = 6
  for ii = 0, n_elements(regions) - 1 do begin
     q = where(ss[*,ii,1] gt nsig * ess[*,ii,1] $
               AND $
               ms[*,ii,1] gt nsig * ems[*,ii,1], compl = uls)
     plotsym, 0, /fill
     oploterror, alog10([ms[q,ii,1]]), alog10([ss[q,ii,1]]), $
                 [ems[q,ii,1]/ms[q,ii,1]/alog(10)], [ess[q,ii,1]/ss[q,ii,1]/alog(10)], $
                 psym = 8, symsize = 1.25, /nohat, errthick = 4, $
                 col = '777777'x, errcol = '777777'x
     oplot, alog10([ms[q,ii,1]]), alog10([ss[q,ii,1]]), psym = 8, col = cols[ii]
     plotsym, 1, thick = 4
     if min(uls) ge 0 then begin
        for jj = 0, N_elements(uls) - 1 do begin
           oplot, alog10(ms[uls[jj],ii,1]) + ems[uls[jj],ii,1]/ms[uls[jj],ii,1]/alog(10) * [-1,1], $
                  replicate(alog10(nsig * [ess[uls[jj],ii,1]]) > [!Y.CRANGE[0] + 0.2], 2), $
                  col = cols[ii], thick = 4
           oplot, alog10([ms[uls[jj],ii,1]]), alog10(nsig * [ess[uls[jj],ii,1]]) > [!Y.CRANGE[0] + 0.2], $
                  psym = 8, col = cols[ii], symsize = 2
        endfor
     endif
  endfor
;  cgtext, !X.CRANGE[0]+0.1, !Y.CRANGE[1]-0.3, /data, $
;          '!18y!X=('+string(median(m[1,*]), f = '(F4.2)')+textoidl('\pm')+$
;          string(stddev(m[1,*], /nan), f = '(F4.2)')+')!18x!X-('+$
;          string(abs(median(b[1,*])), f = '(F4.2)')+textoidl('\pm')+$
;          string(stddev(b[1,*], /nan), f = '(F4.2)')+')', charsize = 0.95, charthick = 3
;  cgtext, !X.CRANGE[0]+0.1, !Y.CRANGE[1]-0.5, /data, $
;          greek('sigma')+'!Dint!N='+string(median(s[1,*]), f = '(F4.2)')+' dex', $
;          charsize = 1, charthick = 4
  device, /close
  spawn, 'gv sedHaSFR.eps &'
  set_plot, 'X'
;  print, pivot
  
  stop
  
end
;plotsedhasfrs, 'study5_lgn_Chains.list', 'study5_exp_Chains.list'
