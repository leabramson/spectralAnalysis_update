pro plotuvj, lgnlist, explist, $
             PLOTBOTH = plotboth

  if NOT keyword_set(PLOTBOTH) then plotboth = 0 else plotboth = 1
  
  summaryData = mrdfits('objectSummaryData.fits', 1)
;  magerr = 2./alog(10) * summaryData.EMU/summaryData.MU
  
  readcol, explist, expfiles, f = 'A'
  readcol, lgnlist, lgnfiles, f = 'A'

  ;;read bc03 00900_1 summary (LEA 20170822)
  ;;made using the best-fit (T0,tau) pair and a_v
  ;; (optical depth = 2.5) w/ default ISM extinct frac.
  sum00900c = mrdfits('00900_1_dust_summary.fits',1)
  mvj = sum00900c.V - sum00900c.J ;; AB mags
  muv = sum00900c.U - sum00900c.V
  q1 = where(sum00900c.AGE le getage(1.019))
  q2 = where(sum00900c.AGE le getage(1.019)+0.5)
  mq = where(abs(sum00900c.AGE - getage(1.019)) le 1.0)
  
  ngals = n_elements(lgnfiles)

  test     = mrdfits(lgnfiles[0], 1)
  regions  = test.REGION
  qui = where(regions eq 'OPTFL' OR $
              regions eq 'INNFL' OR $
              regions eq 'INTER' OR $
              regions eq 'OUTER')
  test = test[qui]
  regions = test.REGION
  nregions = n_elements(regions)

  files = strarr(ngals)

  errs = fltarr(ngals, nregions)

  uv  = fltarr(ngals, nregions, 2)
  vj  = fltarr(ngals, nregions, 2)  
  avs = fltarr(ngals, nregions,2)  
  names = []
  for ii = 0, ngals - 1 do begin
     for jj = 0, 1 do begin
        case jj of
           0: begin
              file = lgnfiles[ii]
              files[ii] = strmid(file, 0, 7)
              case files[ii] of
                 '00900_1': nn = 'SSF'
                 '00451_2': nn = 'CSF'
                 '01916_2': nn = 'PSB'
                 '00660_2': nn = 'PAS'
              endcase
              names = [names, nn]
              qui = where(summaryData.ID eq files[ii])
              sumdata = summaryData[qui]
              modfile = '../'+sumData.FIELD+'/'+sumData.ID+'_pyspecfitLogNormalResults/'
           end
           1: begin
              file = expfiles[ii]
              modfile = '../'+sumData.FIELD+'/'+sumData.ID+'_pyspecfitPhotResults/'
           end
        endcase

        data = mrdfits(file, 1)
        z = data[0].Z[0]

        hit = []        
        for kk = 0, nregions - 1 do begin
           region = regions[kk]
           hit = [hit, where(data.REGION eq region)]
           tfile = modfile+region+'.masked.bestmodel'
           ptfile = modfile+region+'.masked.bestphot'
           readcol, ptfile, l,f,e, /silent
           use = where(l/(1+z) gt 3000)
           err = mean(2.5/alog(10)*e[use]/f[use]) < 0.1
           errs[ii,kk] = err ;sqrt(2) *
;           print, tfile
           foo = getuvj(tfile, z)
           uv[ii,kk,jj] = foo.UV
           vj[ii,kk,jj] = foo.VJ
        endfor

        avs[ii,*,jj] = data[hit].A_V[0]
        
     endfor
  endfor
  ;; set AV to the mean; there's a big systematic offset though
  ;; btween SFHS

  ncolors = 8
  cgloadct, 3, ncolors = ncolors, clip = [10,240], /rev

  avs = reform(avs[*,*,0]) ;mean(avs, dim = 3);mean(avs, dim = 3);
  avs = avs[*,where(regions eq 'OPTFL')]
  minavs = 0
  maxavs = 2
  davs = (maxavs - minavs) / float(ncolors - 1)
  tavs = findgen(ncolors)*davs+minavs
  tcols = []
  for ii = 0, ngals - 1 do begin
     hit = value_locate(tavs, avs[ii]) ; foo = min(abs(tavs - avs[ii]), hit)
     tcols = [tcols, hit[0]]
  endfor
  
  readcol, '../f1.cat', x1, y1, /silent
  readcol, '../g1.cat', x2, y2, /silent
  readcol, '../h1.cat', x3, y3, /silent
  xxx = [x1,x2,x3]
  yyy = [y1,y2,y3]

  slp = (y2[1] - y2[0])/(x2[1] - x2[0])
  cols = ['777777'x, 255, '00a500'x, 'ff5500'x]
  tr = ['OPTFL', 'INNFL', 'INTER', 'OUTER']
  
  set_plot, 'PS'
  xsize = 9
  ysize = 6
  ar = ysize / xsize
  device, filename = 'uvj.eps', $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = xsize, ysize = ysize, /in
  
  plotsym, 0, /fill
  plot, [0.75,1.6-1d-6], [1,2], /xsty, /ysty, /nodat, $
        pos = [0.1,0.20,0.475,0.90], $
        xtitle = '!18V!X-!18J!X', $
        ytitle = '!18U!X-!18V!X', $
        xthick = 3, ythick = 3, $
        charsize = 1.7, charthick = 4
  xxx = [!X.CRANGE[0], x3[1], !X.CRANGE[1] < x2[1], !X.CRANGE[1] < x2[1], !X.CRANGE[0]]
  yyy = [y3[0],        y3[0], y3[0] + slp * (xxx[2] - x3[1]), y1[1] < !Y.CRANGE[1], y1[1] < !Y.CRANGE[1]]
  oplot, x1, y1, thick = 4
  oplot, x2, y2, thick = 4
  oplot, x3, y3, thick = 4
  polyfill, xxx, yyy, $
            /line_fill, col = 'cccccc'x, thick = 1, $
            spacing = 0.05, orien = -15
  cgtext, !X.crange[1] - 0.05, !Y.crange[1] - 0.1, align = 1, $
          'Quiescent', charsize = 1.2, charthick = 4, col = long('0000cc'x)
  cgtext, !X.crange[0] + 0.05, !Y.crange[0] + 0.05, $
          'Starforming', charsize = 1.2, charthick = 4, col = 'aa0000'x
  tot = where(regions eq 'OPTFL')
  oploterror, reform(vj[*,tot,0]), reform(uv[*,tot,0]), $
              reform(errs[*,tot]), reform(errs[*,tot]), $
              psym = 8, symsize = 2.5, /nohat
  for ii = 0, ngals - 1 do begin
     oplot, reform([vj[ii,tot,0]]), reform([uv[ii,tot,0]]), psym = 8, $
            col = cgcolor(string(tcols[ii])), symsize = 2 ;'777777'x,
     cgtext,reform([vj[ii,tot,0]])+0.05, reform([uv[ii,tot,0]])+0.05, $
            names[ii], charsize = 1.25, charthick = 5
  endfor
;  for jj = 0, ngals - 1 do begin
;     for ii = 0, n_elements(tr) - 1 do begin
;        qui = where(regions eq tr[ii])
;        oploterror, reform(vj[jj,qui,0]), reform(uv[jj,qui,0]), $
;                    reform(errs[jj,qui]), reform(errs[jj,qui]), $
;                    psym = 8, symsize = 1.4, /nohat 
;        oplot, reform(vj[jj,qui,0]), reform(uv[jj,qui,0]), psym = 8, col = cols;[ii]
;     endfor
;  endfor
  
  dp =  [0.05,0,0.05,0]
  cgtext, 0.75-dp[0], 0.06, '!18V!X-!18J!X', $
          charsize = 1.7, charthick = 4, /norm, align = 0.5
  cgtext, 0.985, 0.55, '!18U!X-!18V!X', $
          charsize = 1.7, charthick = 4, /norm, orien = 90, align = 0.5
  plot, [0.75,1.6-1d-6], [1+1d-6,2-1d-6], /xsty, /ysty, /nodat, $
        pos = [0.55,0.55,0.75,0.90] - dp, /noer, $
        xtickname = replicate(' ', 60), $
        ytickname = replicate(' ', 60), $
        xthick = 3, ythick = 3, $
        charsize = 1.7, charthick = 4
  oplot, x1, y1, thick = 4
  oplot, x2, y2, thick = 4
  oplot, x3, y3, thick = 4
  polyfill, xxx, yyy, $
            /line_fill, col = 'cccccc'x, thick = 1, $
            spacing = 0.05, orien = -15
  cgtext, !X.CRANGE[0]+0.05, !Y.CRANGE[1]-0.1, /data, $
          names[0], charsize = 1.25, charthick = 5
  for ii = 0, n_elements(tr) - 1 do begin
     qui = where(regions eq tr[ii])
     oploterror, reform(vj[0,qui,0]), reform(uv[0,qui,0]), $
                 reform(errs[0,qui]), reform(errs[0,qui]), $
                 psym = 8, symsize = 1.5, /nohat 
     oplot, reform(vj[0,qui,0]), reform(uv[0,qui,0]), psym = 8, col = cols[ii]
;     if plotboth then $
;        oploterror, reform(vj[0,qui,1]), reform(uv[0,qui,1]), $
;                    reform(errs[0,qui]), reform(errs[0,qui]), $
;                    psym = 8, symsize = 0.7, /nohat, col = cols[ii], $
;                    errcol = cols[ii], errthick = 1
;     oplot, reform(vj[0,qui,1]), reform(uv[0,qui,1]), symsize = 0.7, psym = 8, col = cols[ii]
  endfor
  plot, [0.75,1.6-1d-6], [1+1d-6,2-1d-6], /xsty, /nodat, $
        pos = [0.75,0.55,0.95,0.90] - dp, /noer, $
        xtickname = replicate(' ', 60), $
        ytickname = replicate(' ', 60), ysty = 8+1, $
        xthick = 3, ythick = 3, $
        charsize = 1.7, charthick = 4
  axis, yaxis = 1, yran = !Y.CRANGE, /ysty, $
        ythick = 3, charthick = 4, charsize = 1.7
  oplot, x1, y1, thick = 4
  oplot, x2, y2, thick = 4
  oplot, x3, y3, thick = 4
  polyfill, xxx, yyy, $
            /line_fill, col = 'cccccc'x, thick = 1, $
            spacing = 0.05, orien = -15
  cgtext, !X.CRANGE[0]+0.05, !Y.CRANGE[1]-0.1, /data, $
          names[1], charsize = 1.25, charthick = 5
  for ii = 0, n_elements(tr) - 1 do begin
     qui = where(regions eq tr[ii])
     oploterror, reform(vj[1,qui,0]), reform(uv[1,qui,0]), $
                 reform(errs[1,qui]), reform(errs[1,qui]), $
                 psym = 8, symsize = 1.5, /nohat 
     oplot, [reform(vj[1,qui,0])], [reform(uv[1,qui,0])], psym = 8, col = cols[ii]
;     if plotboth then $
;        oploterror, reform(vj[1,qui,1]), reform(uv[1,qui,1]), $
;                    reform(errs[1,qui]), reform(errs[1,qui]), $
;                    psym = 8, symsize = 0.7, /nohat, col = cols[ii], $
;                    errcol = cols[ii], errthick = 1
;     oplot, reform(vj[1,qui,1]), reform(uv[1,qui,1]), symsize = 0.7, psym = 8, col = cols[ii]
  endfor
  plot, [0.75,1.6-1d-6], [1+1d-6,2-1d-6], /xsty, /ysty, /nodat, $
        pos = [0.55,0.20,0.75,0.55] - dp, /noer, $
;        xtickname = replicate(' ', 60), $
        ytickname = replicate(' ', 60), $
        xthick = 3, ythick = 3, $
        charsize = 1.7, charthick = 4
  oplot, x1, y1, thick = 4
  oplot, x2, y2, thick = 4
  oplot, x3, y3, thick = 4
  polyfill, xxx, yyy, $
            /line_fill, col = 'cccccc'x, thick = 1, $
            spacing = 0.05, orien = -15
  ox = (mvj[q1[-1]] - reform(vj[2,tot,0]))[0] ;; small offsets
  oy = (muv[q1[-1]] - reform(uv[2,tot,0]))[0]
  print, ox, oy
  oplot, mvj[mq] - ox, muv[mq] - oy, thick = 4, linesty = 2, col = 'aaaaaa'x
  cgtext, reform(vj[2,tot,0])+0.225, reform(uv[2,tot,0]) + 0.07, /data, 'Opt. evo track', $
          charsize = 1, charthick = 4, orien = 55, align = 0.5, col = '555555'x
;  oplot, mvj[q2] - ox, muv[q2] - oy, thick = 6
;  oplot, mvj[q1] - ox, muv[q1] - oy, thick = 3, col = 'aaaaaa'x
  cgtext, !X.CRANGE[0]+0.05, !Y.CRANGE[1]-0.1, /data, $
          names[2], charsize = 1.25, charthick = 5
  for ii = 0, n_elements(tr) - 1 do begin
     qui = where(regions eq tr[ii])
     oploterror, reform(vj[2,qui,0]), reform(uv[2,qui,0]), $
                 reform(errs[2,qui]), reform(errs[2,qui]), $
                 psym = 8, symsize = 1.5, /nohat 
     oplot, reform(vj[2,qui,0]), reform(uv[2,qui,0]), psym = 8, col = cols[ii]
;     if plotboth then $
;        oploterror, reform(vj[2,qui,1]), reform(uv[2,qui,1]), $
;                    reform(errs[2,qui]), reform(errs[2,qui]), $
;                    psym = 8, symsize = 0.7, /nohat, col = cols[ii], $
;                    errcol = cols[ii], errthick = 1
;     oplot, reform(vj[2,qui,1]), reform(uv[2,qui,1]), symsize = 0.7, psym = 8, col = cols[ii]
  endfor
  legend, /top, /right, box = 1, /clear, $
          ['Opt', 'Inn', 'Mid', 'Out'], $
          col = cols, psym = 8, charsize = 1.2, $
          charthick = 4, pspacing = 0.5
  
  plot, [0.75,1.6-1d-6], [1+1d-6,2-1d-6], /xsty, /nodat, $
        pos = [0.75,0.20,0.95,0.55] - dp, /noer, $
;        xtickname = replicate(' ', 60), $
        ytickname = replicate(' ', 60), $
        ysty = 8+1, $
        xthick = 3, ythick = 3, $
        charsize = 1.7, charthick = 4
  axis, yaxis = 1, yran = !Y.CRANGE, /ysty, $
        ythick = 3, charsize = 1.7, charthick = 4
  oplot, x1, y1, thick = 4
  oplot, x2, y2, thick = 4
  oplot, x3, y3, thick = 4
  polyfill, xxx, yyy, $
            /line_fill, col = 'cccccc'x, thick = 1, $
            spacing = 0.05, orien = -15
  cgtext, !X.CRANGE[0]+0.05, !Y.CRANGE[1]-0.1, /data, $
          names[3], charsize = 1.25, charthick = 5
  for ii = 0, n_elements(tr) - 1 do begin
     qui = where(regions eq tr[ii])
     oploterror, reform(vj[3,qui,0]), reform(uv[3,qui,0]), $
                 reform(errs[3,qui]), reform(errs[3,qui]), $
                 psym = 8, symsize = 1.5, /nohat
     oplot, reform(vj[3,qui,0]), reform(uv[3,qui,0]), psym = 8, col = cols[ii]
     if plotboth then $
        oploterror, reform(vj[3,qui,1]), reform(uv[3,qui,1]), $
                    reform(errs[3,qui]), reform(errs[3,qui]), $
                    psym = 8, symsize = 0.7, /nohat, col = cols[ii], $
                    errcol = cols[ii], errthick = 1
;     oplot, reform(vj[3,qui,1]), reform(uv[3,qui,1]), symsize = 0.7, psym = 8, col = cols[ii]
  endfor
  cgtext, mean(vj[3,value_locate(regions, tr),0]), $
          1.1, $
          /data, align = 0.5, 'LgNml', charsize = 1, charthick = 2
  cgtext, mean(vj[3,value_locate(regions, tr),1]), $
          1.1, $
          /data, align = 0.5, 'DelExp', charsize = 1, charthick = 2
  cgcolorbar, pos = [0.325,0.225,0.425,0.25] + [0,0.05,0,0.05], $
              charsize = 1.2, charthick = 3, $
              minrange = minavs, maxrange = maxavs, $
              title = 'A!D!18V!X!N', /discrete, $
              xtickint = 1, textthick = 3, ncolors = ncolors, /top
  device, /close
  spawn, 'gv uvj.eps &'
  
  stop
  
end
;plotuvj, 'study5_lgn_Bestfits.list', 'study5_exp_Bestfits.list', /plotboth
