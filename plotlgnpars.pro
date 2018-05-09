function topdf, inarr

  s = size(inarr, /dim)

  outarr = fltarr(s[0], s[1])

  for ii = 0, s[0] - 1 do $
     for jj = 0, s[1] - 1 do $
        outarr[ii,jj] = total(inarr lt inarr[ii,jj]) / float(n_elements(inarr))
  
  RETURN, outarr
end

;;
;;
;;

function getssfr, t0, tau, z

  ts = getage(10)
  te = 25
  dt = 0.1
  time = findgen((te-ts)/dt+1) * dt + t0
  epoch = value_locate(time, getage(z))
  
  sfh = 1./sqrt(2*!pi*tau^2)/time * $
        exp(-0.5 * (alog(time) - t0)^2/tau^2)

  ssfr = sfh[epoch] / total(sfh[0:epoch] * 0.6)
  
  RETURN, ssfr/1d9
end

;;
;;
;;

pro plotlgnpars, lgnlist

  summaryData = mrdfits('objectSummaryData.fits', 1)

  readcol, lgnlist, files, f = 'A', /silent
  nfiles = n_elements(files)

  readcol, 'model_323.dat', mass, z, weight, t0, tau, $
           f = 'F,F,F,F,F', /silent
  mass = thinit(mass, 3)
  t0 = thinit(t0, 3)
  tau = thinit(tau, 3)
  
  cols = [255, '00a500'x, 'ff5500'x]
  plotsym, 0, /fill
  set_plot, 'PS'
  device, filename = 'lgnSFHcovars.eps', $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = 5, ysize = 5, /in
;  multiplot, [2,2], /square, $
;             mxtitle = '!18T!X!D0!N [ln(gyr)]', $
;             mytitle = greek('tau')+' [ln(gyr)]', $
;             mxtitsize = 1.25, mxtitoffset = 1, $
;             mytitsize = 1.25, mytitoffset = 1
  xm = 0.15
  ym = xm
  xw = (1-1.5*xm)/2.
  yw = (1-1.5*ym)/2.
  medt0 = fltarr(nfiles,3)
  emedt0 = fltarr(nfiles,3)
  medtau = fltarr(nfiles,3)
  emedtau = fltarr(nfiles,3)
  for ii = 0, nfiles - 1 do begin

     data = mrdfits(files[ii], 1)
     inn = where(data.REGION eq 'INNFL')
     int = where(data.REGION eq 'INTER')
     out = where(data.REGION eq 'OUTER')
     regions = [inn, int, out]

     for jj = 0, 2 do begin
        medt0[ii,jj] = median(data[regions[jj]].T0)
        emedt0[ii,jj] = stddev(data[regions[jj]].T0, /nan)
        medtau[ii,jj] = median(data[regions[jj]].TAU)
        emedtau[ii,jj] = stddev(data[regions[jj]].TAU, /nan)
     endfor

  endfor

  for ii = 0, nfiles - 1 do begin

     data = mrdfits(files[ii], 1)
     qui = where(summaryData.ID eq strmid(files[ii], 0, 7))   

     ssfrs = []
     for kk = 0, n_elements(t0) - 1 do $
        ssfrs = [ssfrs, getssfr(t0[kk], tau[kk], summaryData[qui].Z)]
     
     case ii of
        0: plot, medt0[ii,*] - emedt0[ii,*], medtau[ii,*] + emedtau[ii,*], $
                 xran = [0,3], $
                 yran = [0,2.5], $
                 /nodat, $
                 xthick = 3, ythick = 3, $
                 charthick = 4, charsize = 1.25, $
                 pos = [xm,ym+yw,xm+xw,ym+2*yw], $
                 xtickname = replicate(' ', 60), xtickint = 1
        1: plot, medt0[ii,*] - emedt0[ii,*], medtau[ii,*] + emedtau[ii,*], /noer, $
                 xran = [1d-6,3], /xsty, $
                 yran = [0,2.5], $
                 /nodat, $
                 xthick = 3, ythick = 3, $
                 charthick = 4, charsize = 1.25, $
                 pos = [xm+xw,ym+yw,xm+2*xw,ym+2*yw], $
                 xtickname = replicate(' ', 60), $
                 ytickname = replicate(' ', 60), xtickint = 1           
        2: plot, medt0[ii,*] - emedt0[ii,*], medtau[ii,*] + emedtau[ii,*], /noer, $
                 xran = [0,3], /xsty, $
;                 yran = [0,0.5], $
                 yran = [0,1.6], /ysty, ytickint = 0.5, $
                 /nodat, $
                 xthick = 3, ythick = 3, $
                 charthick = 4, charsize = 1.25, $
                 pos = [xm,ym,xm+xw,ym+yw], xtickint = 1
        3: plot, medt0[ii,*] - emedt0[ii,*], medtau[ii,*] + emedtau[ii,*], /noer, $
                 xran = [1d-6,3], /xsty, $
;                 yran = [0,0.5], $
                 yran = [0,1.6], /ysty, ytickint = 0.5, $
                 /nodat, $
                 xthick = 3, ythick = 3, $
                 charthick = 4, charsize = 1.25, $
                 pos = [xm+xw,ym,xm+2*xw,ym+yw], $
                 ytickname = replicate(' ', 60), xtickint = 1   
     endcase     
     for jj = 2, 0, -1 do begin
        h = hist_2d(data[regions[jj]].T0, data[regions[jj]].TAU, $
                    min1 = !X.CRANGE[0], max1 = !X.CRANGE[1], bin1 = 0.1, $
                    min2 = !Y.CRANGE[0], max2 = !Y.CRANGE[1], bin2 = 0.1)
;        oplot, data[regions[jj]].T0, data[regions[jj]].TAU, $
;               psym = 8, col = cols[jj], symsize = 0.5
        case jj of
           0: ccols = long(['00ccff'x, '0065ff'x, '0011ff'x])
           1: ccols = long(['00ff00'x, '00a500'x, '005500'x])
           2: ccols = long(['ffcc00'x, 'ff6500'x, 'ff0000'x])
        endcase

        test = topdf(h)
;        cgcontour, h / float(max(h)), nlevels = 3, levels = [0.16,0.50,0.84], label = 0, /onim, $
;                   c_col = ccols, /fill, /cell_fill
        cgcontour, test, nlevels = 3, levels = [0.10,0.50,0.90], label = 0, /onim, $
                   c_col = ccols, /fill, /cell_fill
        endfor
     oplot, t0, tau, col = '777777'x, psym = 8, symsize = 0.25

     thit = where(data.REGION eq 'OPTFL')
     tssfr = alog10(data[thit].SFR[0]) - data[thit].LMASS[0]
     thit2 = where(abs(alog10(ssfrs) - tssfr) le 0.5, NHITS)
;     thit2 = where(alog10(ssfrs) le -10, nhits)
     if nhits gt 0 then $
        oplot, t0[thit2], tau[thit2], psym = 8, symsize = 0.25, col = '333333'x

     for jj = 0, 2 do begin
        oploterror, [medt0[ii,jj]], [medtau[ii,jj]], $
                    emedt0[ii,jj], emedtau[ii,jj], $
                    psym = 8, symsize = 1.5, errthick = 4, /nohat
        oplot, [medt0[ii,jj]], [medtau[ii,jj]], $
               psym = 8, symsize = 1.25, col = cols[jj]
     endfor    

     case strmid(files[ii], 0, 7) of
        '00900_1': nn = 'SSF'
        '00451_2': nn = 'CSF'
        '01916_2': nn = 'PSB'
        '00660_2': nn = 'PAS'
     endcase
     
     cgtext, !X.CRANGE[1] - 0.1, !Y.CRANGE[1] - 0.25, $
             /data, align = 1, $
             nn, charsize = 1, charthick = 4
     oplot, replicate(alog(getage(summaryData[qui].Z)), 2), !Y.CRANGE, linesty = 1, thick = 5

     if ii eq 2 then begin
        legend, pos = [!X.CRANGE[0]+0.1, !Y.CRANGE[1] - 0.025], /data, box = 0, $;/top, /left, box = 0, $
                ['Inner', 'Middle', 'Outer', 'G13 all', greek('Delta')+'!Dlog!18sSFR!X!N<0.5'], $
                psym = 8, col = [cols,'777777'x, '333333'x], $
                pspacing = 0.5, charsize = 1, charthick = 4, $
                symsize = [1,1,1,0.5,0.5]
     endif else if ii eq 3 then begin
;        legend, pos = [!X.crange[0]+0.1, 2], /data, box = 0, /clear, $
;                'age of Uni.', $
;                linesty = 1, thick = 4, $
;                pspacing = 0.5, charsize = 1, charthick = 4
        cgtext, orien = 90, align = 1, $
                alog(getage(summaryData[qui].Z)) - 0.1, !Y.CRANGE[1] - 0.05, /data, $
                'age of Uni.', charsize = 1, charthick = 4;, col = 'aaaaaa'x
     endif     
     
  endfor
  cgtext, mean([xm, xm+2*xw]), 0.05, '!18T!X!D0!N [ln(Gyr)]', $
          charsize = 1.25, charthick = 4, align = 0.5, /norm
  cgtext, 0.025, mean([0.1, 0.9]), greek('tau')+' [ln(Gyr)]', $
          charsize = 1.25, charthick = 4, align = 0.5, orien = 90, /norm  

  device, /close
  set_plot, 'X'
  spawn, 'gv lgnSFHcovars.eps &'

  
;  stop

end
;plotlgnpars, 'study5_lgn_Chains.list'
