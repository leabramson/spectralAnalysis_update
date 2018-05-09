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

pro plotexppars, explist

  summaryData = mrdfits('objectSummaryData.fits', 1)

  readcol, explist, files, f = 'A', /silent
  nfiles = n_elements(files)

  cols = [255, '00a500'x, 'ff5500'x]
  plotsym, 0, /fill
  set_plot, 'PS'
  device, filename = 'expSFHcovars.eps', $
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
        medt0[ii,jj] = median(data[regions[jj]].AGE)
        emedt0[ii,jj] = stddev(data[regions[jj]].AGE, /nan)
        medtau[ii,jj] = median(data[regions[jj]].LTAU)
        emedtau[ii,jj] = stddev(data[regions[jj]].LTAU, /nan)
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
                 xran = [7,10], $
                 yran = [7,10], $
                 /nodat, $
                 xthick = 3, ythick = 3, $
                 charthick = 4, charsize = 1.25, $
                 pos = [xm,ym+yw,xm+xw,ym+2*yw], $
                 xtickname = replicate(' ', 60), xtickint = 1.
        1: plot, medt0[ii,*] - emedt0[ii,*], medtau[ii,*] + emedtau[ii,*], /noer, $
                 xran = [7+1d-6,10], /xsty, $
                 yran = [7,10], $
                 /nodat, $
                 xthick = 3, ythick = 3, $
                 charthick = 4, charsize = 1.25, $
                 pos = [xm+xw,ym+yw,xm+2*xw,ym+2*yw], $
                 xtickname = replicate(' ', 60), $
                 ytickname = replicate(' ', 60), xtickint = 1.      
        2: plot, medt0[ii,*] - emedt0[ii,*], medtau[ii,*] + emedtau[ii,*], /noer, $
                 xran = [7,10], /xsty, $
                 yran = [7,10], $
;                 yran = [0,1.6], /ysty, ytickint = 0.5, $
                 /nodat, $
                 xthick = 3, ythick = 3, $
                 charthick = 4, charsize = 1.25, $
                 pos = [xm,ym,xm+xw,ym+yw], xtickint = 1.
        3: plot, medt0[ii,*] - emedt0[ii,*], medtau[ii,*] + emedtau[ii,*], /noer, $
                 xran = [7+1d-6,10], /xsty, $
                 yran = [7,10], $
;                 yran = [0,1.6], /ysty, ytickint = 0.5, $
                 /nodat, $
                 xthick = 3, ythick = 3, $
                 charthick = 4, charsize = 1.25, $
                 pos = [xm+xw,ym,xm+2*xw,ym+yw], $
                 ytickname = replicate(' ', 60), xtickint = 1.   
     endcase     
     for jj = 2, 0, -1 do begin
        h = hist_2d(data[regions[jj]].AGE, data[regions[jj]].LTAU, $
                    min1 = !X.CRANGE[0], max1 = !X.CRANGE[1], bin1 = 0.1, $
                    min2 = !Y.CRANGE[0], max2 = !Y.CRANGE[1], bin2 = 0.1)
;        oplot, data[regions[jj]].AGE, data[regions[jj]].LTAU, $
;               psym = 8, col = cols[jj], symsize = 0.5
        case jj of
           0: ccols = long(['00ccff'x, '0065ff'x, '0011ff'x])
           1: ccols = long(['00ff00'x, '00a500'x, '005500'x])
           2: ccols = long(['ffcc00'x, 'ff6500'x, 'ff0000'x])
        endcase

        test = topdf(h)
;        cgcontour, h / float(max(h)), nlevels = 3, levels = [0.16,0.50,0.84], label = 0, /onim, $
;                   c_col = ccols, /fill, /cell_fill
        cgcontour, test, nlevels = 3, levels = [0.10,0.50,0.95], label = 0, /onim, $
                   c_col = ccols, /fill, /cell_fill
        endfor

     for jj = 0, 2 do begin
        oploterror, [medt0[ii,jj]], [medtau[ii,jj]], $
                    emedt0[ii,jj], emedtau[ii,jj], $
                    psym = 8, symsize = 1.5, errthick = 4, /nohat
        oplot, [medt0[ii,jj]], [medtau[ii,jj]], $
               psym = 8, symsize = 1.25, col = cols[jj]
     endfor    
;     if ii ne 2 then $

     case strmid(files[ii], 0, 7) of
        '00900_1': nn = 'SSF'
        '00451_2': nn = 'CSF'
        '01916_2': nn = 'PSB'
        '00660_2': nn = 'PAS'
     endcase

     cgtext, !X.CRANGE[0] + 0.1, !Y.CRANGE[1] - 0.25, $
             /data, align = 0, $
             nn, charsize = 1, charthick = 4; $
;     else $
;        cgtext, !X.CRANGE[1] - 0.1, !Y.CRANGE[0] + 0.25, $
;                /data, align = 1, $
;                strmid(files[ii], 0, 7), charsize = 1, charthick = 4 

     if ii eq 3 then begin
        legend, pos = [!X.CRANGE[0]+0.1, mean(!Y.CRANGE) + 0.25], /data, box = 0, $ ;/bottom, /right, box = 0, $;
                ['Inner', 'Middle', 'Outer'], $; 'G13 all', greek('Delta')+'!Dlog!18sSFR!X!N<0.5'], $
                psym = 8, col = [cols], $ ;,'777777'x, '333333'x], $
                pspacing = 0.5, charsize = 1, charthick = 4, $
                symsize = [1,1,1];,0.5,0.5]
     endif else if ii eq 3 then begin
;        legend, pos = [!X.crange[0]+0.1, 2], /data, box = 0, /clear, $
;                'age of Uni.', $
;                linesty = 1, thick = 4, $
;                pspacing = 0.5, charsize = 1, charthick = 4
;        cgtext, orien = 90, align = 1, $
;                alog(getage(summaryData[qui].Z)) - 0.1, !Y.CRANGE[1] - 0.05, /data, $
;                'age of Uni.', charsize = 1, charthick = 4;, col = 'aaaaaa'x
     endif     
     
  endfor
  cgtext, mean([xm, xm+2*xw]), 0.05, 'log age [yr]', $
          charsize = 1.25, charthick = 4, align = 0.5, /norm
  cgtext, 0.025, mean([0.1, 0.9]), 'log '+greek('tau')+'!Dexp!N [yr]', $
          charsize = 1.25, charthick = 4, align = 0.5, orien = 90, /norm  

  device, /close
  set_plot, 'X'
  spawn, 'gv expSFHcovars.eps &'

  
;  stop

end
;plotexppars, 'study5_exp_Chains.list'
