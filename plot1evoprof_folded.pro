pro plot1evoprof_folded, basicFits, evoFits, output

  obsdat = mrdfits(basicFits, 1)
  evodat = mrdfits(evoFits, 1)

  z    = evodat[0].REDSHIFT
  time = evodat[0].TIME
  tobs = evodat[0].TOBS
  
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

  set_plot, 'PS'
  device, filename = output, $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = 10, ysize = 4, /in
  
;  !p.multi = [0,3,0]
  xm = 0.01
  plotsym, 0, /fill
  for jj = 0, nregions - 1 do begin

     edat = evodat[where(evodat.REGION eq regions[jj])]
     bdat = obsdat[where(obsdat.REGION eq regions[jj])]
     
     m  = alog10(edat.MGH[*,1])
     em = 0.5/alog(10) * (edat.MGH[*,2] - edat.MGH[*,0]) / 10.^m

     s  = alog10(edat.SFH[*,1])
     es = 0.5/alog(10) * (edat.SFH[*,2] - edat.SFH[*,0]) / 10.^s

     xxx = [(m - em)[etot], reverse((m + em)[etot])]
     yyy = [(s - es)[etot], reverse((s + es)[etot])]

     area = alog10(bdat.AREA_PHYS)

     print, area
     
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
  
  plot, inn_m, inn_s, /nodat, $
        xtitle = 'log '+greek('Sigma')+'!D!18M!X!L*!N [M!D'+sunsymbol()+'!N kpc!e-2!N]', $
        ytitle = 'log '+greek('Sigma')+'!D!18SFR!X!N [M!D'+sunsymbol()+'!N yr!E-1!N kpc!E-2!N]', $
        charsize = 1.2, charthick = 4, xthick = 3, ythick = 3, $
        pos = [0.075,0.15,0.3,0.9] + [0,0,xm,0], $
        xran = [floor(min([inn_m[edexes],int_m[edexes],out_m[edexes]])), $
                ceil(max([inn_m[edexes],int_m[edexes],out_m[edexes]]))], $
        yran = [floor(min([inn_s[edexes],int_s[edexes],out_s[edexes]])), $
                ceil(max([inn_s[edexes],int_s[edexes],out_s[edexes]]))], $
        ytick_get = yticks
  polyfill, (inn_xxx > !X.CRANGE[0]) < !X.CRANGE[1], $
            (inn_yyy > !Y.CRANGE[0]) < !Y.CRANGE[1], col = cols[0], $
            /line_fill, spacing = 0.05, thick = 1, orien = 45
  oplot, inn_m, inn_s, thick = 7
  oplot, inn_m, inn_s, thick = 4, col = cols[0]
  oplot, inn_m[edexes], inn_s[edexes], psym = 8, symsize = 1.25
  oplot, inn_m[edexes], inn_s[edexes], psym = 8, symsize = 1, col = cols[0]
  polyfill, (int_xxx > !X.CRANGE[0]) < !X.CRANGE[1], $
            (int_yyy > !Y.CRANGE[0]) < !Y.CRANGE[1], col = cols[1], $
            /line_fill, spacing = 0.05, thick = 1, orien = -45
  oplot, int_m, int_s, thick = 7
  oplot, int_m, int_s, thick = 4, col = cols[1]
  oplot, int_m[edexes], int_s[edexes], psym = 8, symsize = 1.25
  oplot, int_m[edexes], int_s[edexes], psym = 8, symsize = 1, col = cols[1]
  polyfill, (out_xxx > !X.CRANGE[0]) < !X.CRANGE[1], $
            (out_yyy > !Y.CRANGE[0]) < !Y.CRANGE[1], $
            col = cols[2], $
            /line_fill, spacing = 0.05, thick = 1, orien = 0
  oplot, out_m, out_s, thick = 7
  oplot, out_m, out_s, thick = 4, col = cols[2]
  oplot, out_m[edexes], out_s[edexes], psym = 8, symsize = 1.25
  oplot, out_m[edexes], out_s[edexes], psym = 8, symsize = 1, col = cols[2]
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
  legend, /bottom, /left, box = 0, $;pos = [!X.CRANGE[0], 0.25*!Y.CRANGE[1]], /data, $
          ['!18<r/r!X!de!N!18>!X=0', ' '+string(rp[1], f = '(F4.2)'), ' '+string(rp[2], f = '(F4.2)')], $
          col = cols, linesty = 0, thick = 6, $
          charsize = 1.1, charthick = 3, pspacing = 0.5, spacing = 1.5
  
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
     
     if regions[ii] eq 'INNFL' then begin
        plot, time, s, /nodat, $
              yran = !Y.CRANGE, $
              xtitle = '!18t!X [Gyr]', $
              ytitle = 'log '+greek('Sigma')+'!D!18SFR!X!N [M!D'+sunsymbol()+'!N yr!E-1!N kpc!E-2!N]', $
              xran = [min(time[edexes])-0.2, max(time[edexes])+0.2], xsty = 8+1, $
              charsize = 1.2, charthick = 4, xthick = 3, ythick = 3, $
              pos = [0.4,0.15,0.6,0.9] + [2*xm,0,3*xm,0], /noer, xtick_get = xticks
        axis, xaxis = 1, $
              charsize = 1.2, charthick = 4, xthick = 3, $
              xtickv = xticks, xticks = n_elements(xticks-1), $
              xtickname = string(z[value_locate(time,xticks)], f = '(F4.2)'), $
              xtitle = '!18z!X'
     endif
        
     yyy = [(s-es)[etot], reverse((s+es)[etot])]
     polyfill, (xxx > !X.CRANGE[0]) < !X.CRANGE[1], $
               (yyy > !Y.CRANGE[0]) < !Y.CRANGE[1], col = cols[ii], $
               /line_fill, thick = 1, spacing = 0.05, orien = ii * 10
     oplot, time, s, thick = 6
     oplot, time, s, col = cols[ii], thick = 3
     
  endfor
  oplot, replicate(tobs, 2), !Y.CRANGE, linesty = 5, thick = 4
  cgtext, tobs-0.15, !Y.CRANGE[0]+0.15, /data, $
          '!18t!X!Dobs!N', orien = 90, charsize = 1.1, charthick = 3
  
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
  
  plot, [-0.1,2], minmax(masses - anchms), /xsty, $
        xtitle = '!18r/r!S!X!De!N!R!Eobs!N', $
        ytitle = 'log '+greek('Sigma')+'!D!18M!X!L*!N(!18r!X,!18t!X)!18/!XM!S!D*,inner!R!Eobs!N', $
        /nodat, $
        charsize = 1.2, charthick = 4, xthick = 3, ythick = 3, $
        pos = [0.7,0.15,0.9,0.9] + [4*xm,0,5*xm,0], /noer
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
  legend, /top, /right, box = 0, $
          tags, $
          col = [cgcolor('0'), 0, cgcolor(string(fix(nepochs-1)))], $
          linesty = 0, thick = 6, pspacing = 1, $
          charsize = 1.1, charthick = 3, spacing =1.5

  device, /close
  spawn, 'gv '+output+' &'
  
end
;plot1evoprof_folded, '00900_1_lgn_bestfit.fits','00900_1_lgn_recons.fits', '00900_1_evostuff.eps'
;plot1evoprof_folded, '00660_2_lgn_bestfit.fits', '00660_2_lgn_recons.fits', '00660_2_evostuff.eps'
;plot1evoprof_folded, '00451_2_lgn_bestfit.fits',
;'00451_2_lgn_recons.fits', '00451_2_evostuff.eps'

;plot1evoprof_folded, '01266_1_lgn_bestfit.fits',
;'01266_1_lgn_recons.fits', '01266_1_evostuff.eps'
;plot1evoprof_folded, '01931_1_lgn_bestfit.fits', '01931_1_lgn_recons.fits', '01931_1_evostuff.eps'
;plot1evoprof_folded, '01916_2_lgn_bestfit.fits', '01916_2_lgn_recons.fits', '01916_2_evostuff.eps'
