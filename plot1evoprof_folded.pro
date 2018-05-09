pro plot1evoprof_folded, basicFits, evoFits, output

  obsdat = mrdfits(basicFits, 1)
  evodat = mrdfits(evoFits, 1)

;  name = strmid(basicFits, 0, 7)

  case strmid(basicFits, 0, 7) of
     '00900_1': name = 'SSF'
     '00451_2': name = 'CSF'
     '01916_2': name = 'PSB'
     '00660_2': name = 'PAS'
  endcase
  
  z    = evodat[0].REDSHIFT
  time = evodat[0].TIME
  tobs = evodat[0].TOBS
  
  obsdex = value_locate(time, tobs) ;; LEA 20180412
  edexes = value_locate(time, tobs+[-1,0,1])
  etot = where(time ge tobs-1 and time le tobs+1)
  nepochs = n_elements(edexes)
  
  regions = ['INNFL', 'INTER', 'OUTER']
  cols = [255, '00a500'x, 'ff5500'x] ;'00aa00'x, 'ffa500'x, 
  nregions = n_elements(regions)

  tdata = obsdat[where(obsdat.REGION eq 'INNFL' $
                       OR obsdat.REGION eq 'INTER' $
                       OR obsdat.REGION eq 'OUTER')]
  r = tdata.RNORM
  r = r[sort(r)]

  cgloadct, 33, ncolors = nepochs, clip = [20,200], /rev
  
;  !p.multi = [0,3,0]
  xm = 0.0
  plotsym, 0, /fill
  for jj = 0, nregions - 1 do begin

     edat = evodat[where(evodat.REGION eq regions[jj])]
     bdat = obsdat[where(obsdat.REGION eq regions[jj])]
     
     m  = alog10(edat.MGH[*,1])
     em = 0.5/alog(10) * (edat.MGH[*,2] - edat.MGH[*,0]) / 10.^m

     s  = alog10(edat.SFH[*,1])
     es = 0.5/alog(10) * (edat.SFH[*,2] - edat.SFH[*,0]) / 10.^s

     spk = max(s)
     
     xxx = [(m - em)[etot], reverse((m + em)[etot])]
     yyy = [(s - es)[etot], reverse((s + es)[etot])]

     area = alog10(bdat.AREA_PHYS)

;     print, area
     
     if regions[jj] eq 'INNFL' then begin
        mgh_inner = m
        mgh_inner_err = em
        inn_area = area
        sfh_inner = s - area
        sfh_inner_err = es        
        inn_xxx = xxx - area
        inn_yyy = yyy - area
        inn_m   = m - area
        inn_s   = s - area
     endif else if regions[jj] eq 'INTER' then begin
        mgh_inter = m
        mgh_inter_err = em
        int_area = area
        sfh_inter = s - area
        sfh_inter_err = es
        int_xxx = xxx - area
        int_yyy = yyy - area
        int_m   = m - area
        int_s   = s - area
     endif else if regions[jj] eq 'OUTER' then begin
        mgh_outer = m
        mgh_outer_err = em
        out_area = area
        sfh_outer = s - area
        sfh_outer_err = es
        out_xxx = xxx - area
        out_yyy = yyy - area
        out_m   = m - area
        out_s   = s - area
     endif
  endfor

;  stop
  
  set_plot, 'PS'
  device, filename = output, $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = 14., ysize = 4, /in

  plot, inn_m, inn_s, /nodat, $
        xtitle = 'log '+greek('Sigma')+$
        '!D!18M!X!L*!N [M!D'+sunsymbol()+'!N kpc!e-2!N]', $
        ytitle = 'log '+greek('Sigma')+$
        '!D!18SFR!X!N [M!D'+sunsymbol()+'!N yr!E-1!N kpc!E-2!N]', $
        charsize = 1.4, charthick = 5, xthick = 4, ythick = 4, $
        pos = [0.1,0.2,0.325,0.85] + [0,0,xm,0], $
        xran = [floor(min([inn_m[edexes], $
                           int_m[edexes],$
                           out_m[edexes]])), $
                ceil(max([inn_m[edexes],$
                          int_m[edexes],$
                          out_m[edexes]]))-1d-6], /xsty, $
        yran = [floor(min([inn_s[edexes],$      
                           int_s[edexes],$      
                           out_s[edexes]])), $  
                max([inn_s,$           ;[edexes]
                     int_s,$           ;[edexes]
                     out_s]) + 0.3], $      ;[edexes]
        ytick_get = yticks
  foo = max(sfh_inner, pkdex)
  polyfill, inn_m[pkdex] + mgh_inner_err[pkdex] * [-1,-1,1,1], $
            [!Y.CRANGE[0], !Y.CRANGE[1], $
             !Y.CRANGE[1], !Y.CRANGE[0]], $
            col = '0055ff'x
  oplot, replicate(inn_m[pkdex], 2), !Y.CRANGE, $
         linesty = 2, thick = 4
  polyfill, (inn_xxx > !X.CRANGE[0]) < !X.CRANGE[1], $
            (inn_yyy > !Y.CRANGE[0]) < !Y.CRANGE[1], col = cols[0], $
            /line_fill, spacing = 0.05, thick = 1, orien = 45
  oplot, inn_m[0:obsdex], inn_s[0:obsdex], thick = 7
  oplot, inn_m[0:obsdex], inn_s[0:obsdex], thick = 4, col = cols[0]
  oplot, inn_m[obsdex+1:*], inn_s[obsdex+1:*], thick = 7, linesty = 2
  oplot, inn_m[obsdex+1:*], inn_s[obsdex+1:*], thick = 4, col = cols[0], linesty = 2
  psyms = [2,8,4]
  for ii = 0, nepochs - 1 do begin
     oplot, [inn_m[edexes[ii]]], [inn_s[edexes[ii]]], psym = psyms[ii], symsize = 1.7, thick = 6
     oplot, [inn_m[edexes[ii]]], [inn_s[edexes[ii]]], psym = psyms[ii], symsize = 1.25, col = cols[0], thick = 4
  endfor
  polyfill, (int_xxx > !X.CRANGE[0]) < !X.CRANGE[1], $
            (int_yyy > !Y.CRANGE[0]) < !Y.CRANGE[1], col = cols[1], $
            /line_fill, spacing = 0.05, thick = 1, orien = -45
  oplot, int_m[0:obsdex], int_s[0:obsdex], thick = 7
  oplot, int_m[0:obsdex], int_s[0:obsdex], thick = 4, col = cols[1]
  oplot, int_m[obsdex+1:*], int_s[obsdex+1:*], thick = 7, linesty = 2
  oplot, int_m[obsdex+1:*], int_s[obsdex+1:*], thick = 4, col = cols[1], linesty = 2

  for ii = 0, nepochs - 1 do begin
     oplot, [int_m[edexes[ii]]], [int_s[edexes[ii]]], psym = psyms[ii], symsize = 1.7, thick = 6
     oplot, [int_m[edexes[ii]]], [int_s[edexes[ii]]], psym = psyms[ii], symsize = 1.25, col = cols[1], thick = 4
  endfor
  polyfill, (out_xxx > !X.CRANGE[0]) < !X.CRANGE[1], $
            (out_yyy > !Y.CRANGE[0]) < !Y.CRANGE[1], $
            col = cols[2], $
            /line_fill, spacing = 0.05, thick = 1, orien = 0
  oplot, out_m[0:obsdex], out_s[0:obsdex], thick = 7
  oplot, out_m[0:obsdex], out_s[0:obsdex], thick = 4, col = cols[2]
  oplot, out_m[obsdex+1:*], out_s[obsdex+1:*], thick = 7, linesty = 2
  oplot, out_m[obsdex+1:*], out_s[obsdex+1:*], thick = 7, col = cols[2], linesty = 2
  for ii = 0, nepochs - 1 do begin
     oplot, [out_m[edexes[ii]]], [out_s[edexes[ii]]], psym = psyms[ii], symsize = 1.7, thick = 6             
     oplot, [out_m[edexes[ii]]], [out_s[edexes[ii]]], psym = psyms[ii], symsize = 1.25, col = cols[2], thick = 4
  endfor
  for ii = 0, nepochs - 1 do begin
     if ii ne 1 then $
        col = cgcolor(string(fix(ii))) $
     else  $
        col = 0
     oplot, [inn_m[edexes[ii]], int_m[edexes[ii]], out_m[edexes[ii]]], $
            [inn_s[edexes[ii]], int_s[edexes[ii]], out_s[edexes[ii]]], $
            thick = 4, linesty = 1, col = col
  endfor
  rp = r
  if name eq 'SSF' then $
     cgtext, !X.CRANGE[0] + 0.2, yticks[-1] - 0.5 * (yticks[1] - yticks[0]), '!18'+name+'!X', $
             charsize = 1.75, charthick = 5, /data $
  else $
     cgtext, !X.CRANGE[1] - 0.2, yticks[0] + 0.5 * (yticks[1] - yticks[0]), '!18'+name+'!X', $
             charsize = 1.75, charthick = 5, /data, align = 1
  legend, pos = [!X.WINDOW[0], !Y.WINDOW[0]+0.2], /norm, box = 0, $ ;, $
          ['Inner', 'Middle', 'Outer'], $
;          ['!18<r/r!X!de!N!18>!X=0', ' '+string(rp[1], f = '(F4.2)'), ' '+string(rp[2], f = '(F4.2)')], $
          col = cols, linesty = 0, thick = 10, $
          charsize = 1.5, charthick = 4, $
          pspacing = 0.75, spacing = 1.5
  
  xxx = [time[etot], reverse(time[etot])]  
  for ii = 0, nregions - 1 do begin

     case regions[ii] of
        'INNFL': begin
           s  = sfh_inner
           es = sfh_inner_err
        end
        'INTER': begin
           s  = sfh_inter
           es = sfh_inter_err
        end
        'OUTER': begin
           s  = sfh_outer
           es = sfh_outer_err
        end
     endcase

     xw = !X.WINDOW[1] - !X.WINDOW[0]
     
     if regions[ii] eq 'INNFL' then begin
        plot, time, s, /nodat, $
              yran = !Y.CRANGE, $
              xtitle = '!18t!X [Gyr]', $
;              ytitle = 'log '+greek('Sigma')+'!D!18SFR!X!N [M!D'+sunsymbol()+'!N yr!E-1!N kpc!E-2!N]', $
              ytickname = replicate(' ', 60), $
              xran = [min(time[edexes])-0.1, max(time[edexes])+0.2], xsty = 8+1, $
              charsize = 1.5, charthick = 5, xthick = 4, ythick = 4, $
              pos = [!X.WINDOW[1],!Y.WINDOW[0],!X.WINDOW[1]+xw,!Y.WINDOW[1]], /noer, xtick_get = xticks
        axis, xaxis = 1, $
              charsize = 1.5, charthick = 5, xthick = 4, $
              xtickv = xticks, xticks = n_elements(xticks-1), $
              xtickname = string(z[value_locate(time,xticks)], f = '(F4.2)'), $
              xtitle = '!18z!X'
     endif
        
     yyy = [(s-es)[etot], reverse((s+es)[etot])]
     polyfill, (xxx > !X.CRANGE[0]) < !X.CRANGE[1], $
               (yyy > !Y.CRANGE[0]) < !Y.CRANGE[1], col = cols[ii], $
               /line_fill, thick = 1, spacing = 0.05, orien = ii * 10
     oplot, time[0:obsdex], s[0:obsdex], thick = 6
     oplot, time[0:obsdex], s[0:obsdex], col = cols[ii], thick = 3
     oplot, time[obsdex+1:*], s[obsdex+1:*], thick = 6, linesty = 2
     oplot, time[obsdex+1:*], s[obsdex+1:*], thick = 6, col = cols[ii], linesty = 2
     for jj = 0, nepochs - 1 do begin
        oplot, [time[edexes[jj]]], [s[edexes[jj]]], psym = psyms[jj], symsize = 1.7, thick = 6
        oplot, [time[edexes[jj]]], [s[edexes[jj]]], psym = psyms[jj], symsize = 1.25, col = cols[ii], thick = 3
     endfor
     
  endfor
  oplot, replicate(tobs, 2), !Y.CRANGE, linesty = 5, thick = 4
  cgtext, tobs-0.15, !Y.CRANGE[0]+0.25, /data, $
          '!18t!X!Dobs!N', orien = 90, charsize = 1.5, charthick = 4
  
  masses   = dblarr(nregions,nepochs)
  masserrs = dblarr(nregions,nepochs)
  tmass    = dblarr(nepochs)
  tmasserr = dblarr(nepochs)

  areas = [inn_area,int_area,out_area]
  
  xxx = [r, reverse(r)]
  for jj = 0, nregions - 1 do begin

     case regions[jj] of
        'INNFL': begin
           m = inn_m
           em = mgh_inner_err
        end
        'INTER':begin
           m = int_m
           em = mgh_inter_err
        end
        'OUTER':begin
           m = out_m
           em = mgh_outer_err
        end
     endcase
     
     masses[jj,*]   = m[edexes]
     masserrs[jj,*] = em[edexes]
     
  endfor

  anchms = masses[0,1]
  
  plot, [-0.1,1.75], minmax(masses - anchms), /xsty, $
;        xtitle = '!18r/r!S!X!De!N!R!Eobs!N', $
        ytitle = 'log '+greek('Sigma')+'!D!18M!X!L*!N(!18r!X,!18t!X)!18/'+$
        greek('Sigma')+'!S!DM!L!X*,inner!R!Eobs!N', $
        /nodat, $
        charsize = 1.5, charthick = 5, $
        xthick = 4, ythick = 4, $
        pos = [!X.WINDOW[1]+0.1,!Y.WINDOW[0],$
               !X.WINDOW[1]+xw/2.2+0.1,!Y.WINDOW[1]], /noer
  cgtext, !X.WINDOW[1]+0.01, 0.02, /norm, $
          '!18r/r!S!X!De!N!R!Eobs!N', charsize = 1.5, charthick = 4, align = 0.5
  for ii = 0, nepochs - 1 do begin
     anchm = anchms;[ii]
;     ii = ints[jj]

;     tmass[ii]    = alog10(total(10.^masses[*,ii]))
;     tmasserr[ii] = sqrt(total(masserrs[*,ii]^2))
     
     ;; From wolfram alpha. need sigma(f=log(x)-log(x+y+z))^2
     ;; = (df/dx)^2*sigma_x^2
     ;; = [(y+z) / x(log(10))(x+y+z)]^2 * sigma_x^2 = (sigma_x/x)^2 *
     ;; [(y+z)/(x+y+z)]^2 = sigma_ln(x)^2 * some stuff

;     em1 = (10.^masserrs[0,ii] - 1); * alog(10) 
;     em2 = (10.^masserrs[1,ii] - 1); * alog(10) 
;     em3 = (10.^masserrs[2,ii] - 1); * alog(10) 
;     em1 = sqrt((em1 * total(10.^masses[1:2,ii]))^2 / total(10.^masses[*,ii])^2)
;     em2 = sqrt((em2 * total(10.^masses[[0,2],ii]))^2 / total(10.^masses[*,ii])^2)
;     em3 = sqrt((em3 * total(10.^masses[0:1,ii]))^2 / total(10.^masses[*,ii])^2)

;     err = [em1,em2,em3]
     
     err = masserrs[*,ii]  
     
     yyy = [masses[*,ii] - err, $
            reverse(masses[*,ii] + err)] - anchm;- tmass[ii]
     case ii of
        0: early_yyy = yyy
        1: obs_yyy   = yyy
        2: late_yyy  = yyy
     endcase
  endfor
  polyfill, xxx, (early_yyy > !Y.CRANGE[0]) < !Y.CRANGE[1], col = cgcolor(string(fix(0))), $
            /line_fill, thick = 1, spacing = 0.05, orien = 45
  polyfill, xxx, (late_yyy > !Y.CRANGE[0]) < !Y.CRANGE[1], col = cgcolor(string(fix(2))), $
            /line_fill, thick = 1, spacing = 0.05, orien = -45
  polyfill, xxx, (obs_yyy > !Y.CRANGE[0]) < !Y.CRANGE[1], col = '777777'x
  oplot, r, masses[*,2] - anchm, thick = 7                                   ;tmass[2]s[2]
  oplot, r, masses[*,2] - anchm, thick = 4, col = cgcolor(string(fix(2)))    ;tmass[2]s[2]
  oplot, r, masses[*,0] - anchm, thick = 7                                   ;tmass[0]s[0]
  oplot, r, masses[*,0] - anchm, thick = 4, col = cgcolor(string(fix(0)))    ;tmass[0]s[0]
  oplot, r, masses[*,1] - anchm, thick = 7                                   ;tmass[1]s[1]

  qui = where(obsdat.REGION eq 'INNFL')
  zobs = obsdat[qui].Z[0]
  z1 = z[edexes[0]]
  z2 = z[edexes[-1]]
  tags = ['!18z!X='+string(z1, f = '(F3.1)'), $
          string(zobs, f = '(F3.1)')+' (!18z!X!Dobs!N)', $
          string(z2, f = '(F3.1)')]
  legend, /right, pos = [!X.CRANGE[1], !Y.CRANGE[1]], /data, box = 0, $
          tags, $
          col = [cgcolor('0'), 0, cgcolor(string(fix(nepochs-1)))], $
          linesty = 0, thick = 10, pspacing = 0.5, $
          charsize = 1.25, charthick = 4, spacing =1.75

  d12 = masses[*,0] - masses[*,1]
  e12 = sqrt(masserrs[*,0]^2 + masserrs[*,1]^2)
  d23 = masses[*,2] - masses[*,1]
  e23 = sqrt(masserrs[*,1]^2 + masserrs[*,2]^2) 
  
  plot, [-0.1,1.75], [min(d12)-0.1,max(d23)+0.1], /xsty, $
;        xtitle = '!18r/r!S!X!De!N!R!Eobs!N', $
;        ytitle = 'log '+greek('Sigma')+'!D!18M!X!L*!N(!18r!X,!18t!X)!18/'+$
;        greek('Sigma')+'!S!DM!L!X*,inner!R!Eobs!N', $
        yminor = 2, $
        ysty = 8+1, $
        ytickname = replicate(' ',60), $
        /nodat, $
        charsize = 1.5, charthick = 5, $
        xthick = 4, ythick = 4, $
        pos = [!X.WINDOW[1]+0.02,!Y.WINDOW[0],$
               !X.WINDOW[1]+xw/2.2+0.025,!Y.WINDOW[1]], /noer
  axis, yaxis = 1, charsize = 1.5, charthick = 5, ythick = 4, $
        ytitle = 'log '+greek('Sigma')+'!D!18M!X!L*!N(!18r!X,!18t!X)!18/'+$
        greek('Sigma')+'!S!DM!L!X*!R!Eobs!N(!18r!X)', ysty = 1, yminor = 2
  polyfill, xxx, ([-masserrs[*,1], reverse(masserrs[*,1])] > !Y.CRANGE[0]) < !Y.CRANGE[1], $
            col = '777777'x
  polyfill, xxx, ([d12 - e12, reverse(d12 + e12)] > !Y.CRANGE[0]) < !Y.CRANGE[1], $
            col = cgcolor(string(fix(0))), $
            /line_fill, thick = 1, spacing = 0.05, orien = 45
  polyfill, xxx, ([d23 - e23, reverse(d23 + e23)] > !Y.CRANGE[0]) < !Y.CRANGE[1], $
            col = cgcolor(string(fix(2))), $
            /line_fill, thick = 1, spacing = 0.05, orien = 45
  oplot, r, masses[*,2] - masses[*,1], thick = 7
  oplot, r, masses[*,2] - masses[*,1], thick = 4, col = cgcolor(string(fix(2)))
  oplot, r, masses[*,0] - masses[*,1], thick = 7
  oplot, r, masses[*,0] - masses[*,1], thick = 4, col = cgcolor(string(fix(0)))
  oplot, r, replicate(0, 3), thick = 7
  
  device, /close
;  spawn, 'gv '+output+' &'
  set_plot, 'X'
  
  print, '- - - - - - - - '
  print, time[pkdex], z[pkdex]
  print, inn_s[pkdex], spk, sfh_inner_err[pkdex]
  print, inn_m[pkdex], mgh_inner_err[pkdex]
  print, '- - - - - - - - '

  
end

pro dosamp
  plot1evoprof_folded, '00900_1_lgn_bestfit.fits', '00900_1_lgn_recons.fits', '00900_1_evostuff.eps'
  plot1evoprof_folded, '00660_2_lgn_bestfit.fits', '00660_2_lgn_recons.fits', '00660_2_evostuff.eps'
  plot1evoprof_folded, '00451_2_lgn_bestfit.fits', '00451_2_lgn_recons.fits', '00451_2_evostuff.eps'
  plot1evoprof_folded, '01916_2_lgn_bestfit.fits', '01916_2_lgn_recons.fits', '01916_2_evostuff.eps'
end
;plot1evoprof_folded, '01266_1_lgn_bestfit.fits',
;'01266_1_lgn_recons.fits', '01266_1_evostuff.eps'
;plot1evoprof_folded, '01931_1_lgn_bestfit.fits', '01931_1_lgn_recons.fits', '01931_1_evostuff.eps'

