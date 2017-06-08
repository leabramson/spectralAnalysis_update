pro plot1evoprof_folded, basicFits, evoFits, output

  obsdat = mrdfits(basicFits, 1)
  evodat = mrdfits(evoFits, 1)

  z    = evodat.INNFL_REDSHIFT
  time = evodat.INNFL_TIME
  tobs = evodat.INNFL_TOBS
  
  edexes = value_locate(time, tobs+[-1,0,1])
  etot = where(time ge tobs-1 and time le tobs+1)
  nepochs = n_elements(edexes)
  
  regions = ['INNFL', 'INTER', 'OUTER']
  cols = [255, '007700'x, 'ff5500'x] ;'00aa00'x, 'ffa500'x, 
  nregions = n_elements(regions)

  tdata = obsdat[where(obsdat.REGION ne 'OPTFL')]
  r = tdata.RNORM
  r = r[sort(r)]
  r = r[2:*]

  cgloadct, 33, ncolors = nepochs, clip = [20,200], /rev
  
  set_plot, 'PS'
  device, filename = output, $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = 10, ysize = 4, /in
  
;  !p.multi = [0,3,0]
  xm = 0.01
  plotsym, 0, /fill
  for jj = 0, nregions - 1 do begin
     case jj of
        0: begin
           
           m  = 0.5 * (evodat.INNFL_MGH[*,1] + evodat.INNFL_MGH[*,1])
           s  = 0.5 * (evodat.INNFL_SFH[*,1] + evodat.INNFL_SFH[*,1])
           
           em = 0.5 * (evodat.INNFL_MGH[*,2] - evodat.INNFL_MGH[*,0]) / m / alog(10)
           es = 0.5 * (evodat.INNFL_SFH[*,2] - evodat.INNFL_SFH[*,0]) / s / alog(10)

           mgh_inner = m
           mgh_inner_err = em

           sfh_inner = s
           sfh_inner_err = es
           
        end
        1: begin

           m = 0.5 * (evodat.INTUP_MGH[*,1] + evodat.INTDN_MGH[*,1])
           s = 0.5 * (evodat.INTUP_SFH[*,1] + evodat.INTDN_SFH[*,1])
           
           em1 = 0.5 * (evodat.INTUP_MGH[*,2] - evodat.INTUP_MGH[*,0]) / m / alog(10)
           em2 = 0.5 * (evodat.INTDN_MGH[*,2] - evodat.INTDN_MGH[*,0]) / m / alog(10)
           em  = sqrt(em1^2 + em2^2) / sqrt(2.)

           es1 = 0.5 * (evodat.INTUP_SFH[*,2] - evodat.INTUP_SFH[*,0]) / s / alog(10)
           es2 = 0.5 * (evodat.INTDN_SFH[*,2] - evodat.INTDN_SFH[*,0]) / s / alog(10)
           es  = sqrt(es1^2 + es2^2) / sqrt(2.)

           mgh_inter = m
           mgh_inter_err = em

           sfh_inter = s
           sfh_inter_err = es
           
        end
        2: begin

           m = 0.5 * (evodat.OUTUP_MGH[*,1] + evodat.OUTDN_MGH[*,1])
           s = 0.5 * (evodat.OUTUP_SFH[*,1] + evodat.OUTDN_SFH[*,1])
           
           em1 = 0.5 * (evodat.OUTUP_MGH[*,2] - evodat.OUTUP_MGH[*,0]) / m / alog(10)
           em2 = 0.5 * (evodat.OUTDN_MGH[*,2] - evodat.OUTDN_MGH[*,0]) / m / alog(10)
           em  = sqrt(em1^2 + em2^2) / sqrt(2.)

           es1 = 0.5 * (evodat.OUTUP_SFH[*,2] - evodat.OUTUP_SFH[*,0]) / s / alog(10)
           es2 = 0.5 * (evodat.OUTDN_SFH[*,2] - evodat.OUTDN_SFH[*,0]) / s / alog(10)
           es  = sqrt(es1^2 + es2^2) / sqrt(2.)

           mgh_outer = m
           mgh_outer_err = em

           sfh_outer = s
           sfh_outer_err = es

        end
     endcase

     m = alog10(m)
     s = alog10(s)
     
     xxx = [(m - em)[etot], reverse((m + em)[etot])]
     yyy = [(s - es)[etot], reverse((s + es)[etot])]

     case jj of
        0: begin
           inn_xxx = xxx
           inn_yyy = yyy
           inn_m   = m
           inn_s   = s
        end
        1: begin
           int_xxx = xxx
           int_yyy = yyy
           int_m   = m
           int_s   = s
        end
        2: begin
           out_xxx = xxx
           out_yyy = yyy
           out_m   = m
           out_s   = s
        end
     endcase
  endfor
  plot, inn_m, inn_s, /nodat, $
        xtitle = 'log !18M!X!D*!N [M!D'+sunsymbol()+'!N]', $
        ytitle = 'log !18SFR!X [M!D'+sunsymbol()+'!N yr!E-1!N]', $
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
     oplot, alog10(([mgh_inner[edexes[ii]], mgh_inter[edexes[ii]], mgh_outer[edexes[ii]]])), $
            alog10(([sfh_inner[edexes[ii]], sfh_inter[edexes[ii]], sfh_outer[edexes[ii]]])), $
            thick = 4, linesty = 1, col = col
  endfor
  rp = r
  legend, /bottom, /left, box = 0, $;pos = [!X.CRANGE[0], 0.25*!Y.CRANGE[1]], /data, $
          ['!18<r/r!X!de!N!18>!X=0', ' '+string(rp[1], f = '(F4.2)'), ' '+string(rp[2], f = '(F4.2)')], $
          col = cols, linesty = 0, thick = 6, $
          charsize = 1.1, charthick = 3, pspacing = 0.5, spacing = 1.5

  xxx = [time[etot], reverse(time[etot])]  
  for ii = 0, n_elements(regions) - 1 do begin
     case ii of
        0: begin
           s  = evodat.INNFL_SFH[*,1]
           es = 0.5 * (evodat.INNFL_SFH[*,2] - evodat.INNFL_SFH[*,0]) / s / alog(10)

;           q = where(obsdat.REGION eq 'INNFL')
;           so  = obsdat[q].SFR[0]
;           eso = 0.5*(obsdat[q].SFR[2] - obsdat[q].SFR[1]) / so / alog(10)           
        end
        1: begin
           s   = 0.5 * (evodat.INTUP_SFH[*,1] + evodat.INTDN_SFH[*,1])
           es1 = 0.5 * (evodat.INTUP_SFH[*,2] - evodat.INTUP_SFH[*,0]) / s / alog(10)
           es2 = 0.5 * (evodat.INTDN_SFH[*,2] - evodat.INTDN_SFH[*,0]) / s / alog(10)
           es  = sqrt(es1^2 + es2^2) / sqrt(2.)

;           q1  = where(obsdat.REGION eq 'INTUP')
;           q2  = where(obsdat.REGION eq 'INTDN')
;           so1 = obsdat[q1].SFR[0]
;           so2 = obsdat[q2].SFR[0]
;           eso1 = 0.5*(obsdat[q1].SFR[2] - obsdat[q1].SFR[1]) / so1 / alog(10)
;           eso2 = 0.5*(obsdat[q2].SFR[2] - obsdat[q2].SFR[1]) / so2 / alog(10)
;           so = 0.5 * (so1 + so2)
;           eso = sqrt(eso1^2+eso2^2) / sqrt(2.)
        end
        2: begin
           s   = 0.5 * (evodat.OUTUP_SFH[*,1] + evodat.OUTDN_SFH[*,1])
           es1 = 0.5 * (evodat.OUTUP_SFH[*,2] - evodat.OUTUP_SFH[*,0]) / s / alog(10)
           es2 = 0.5 * (evodat.OUTDN_SFH[*,2] - evodat.OUTDN_SFH[*,0]) / s / alog(10)
           es  = sqrt(es1^2 + es2^2) / sqrt(2.)

;           q1  = where(obsdat.REGION eq 'OUTUP')
;           q2  = where(obsdat.REGION eq 'OUTDN')
;           so1 = obsdat[q1].SFR[0]
;           so2 = obsdat[q2].SFR[0]
;           eso1 = 0.5*(obsdat[q1].SFR[2] - obsdat[q1].SFR[1]) / so1 / alog(10)
;           eso2 = 0.5*(obsdat[q2].SFR[2] - obsdat[q2].SFR[1]) / so2 / alog(10)
;           so = 0.5 * (so1 + so2)
;           eso = sqrt(eso1^2+eso2^2) / sqrt(2.)
        end
     endcase     
     s = alog10(s)

     if ii eq 0 then begin
        plot, time, alog10(evodat.INNFL_SFH[*,1]), /nodat, $
              yran = !Y.CRANGE, $
              xtitle = '!18t!X [Gyr]', $
              ytitle = 'log !18SFR!X [M!D'+sunsymbol()+'!N yr!E-1!N]', $
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
;     oploterror, [tobs], [alog10(so)], [eso], psym = 8, symsize = 1, $
;                 errcol = cols[ii], col = cols[ii], errthick = 4
     
  endfor
  oplot, replicate(tobs, 2), !Y.CRANGE, linesty = 5, thick = 4
  cgtext, tobs-0.15, !Y.CRANGE[0]+0.15, /data, $
          '!18t!X!Dobs!N', orien = 90, charsize = 1.1, charthick = 3
  
  masses   = fltarr(nregions,nepochs)
  masserrs = fltarr(nregions,nepochs)
  tmass    = fltarr(nepochs)
  tmasserr = fltarr(nepochs)

  xxx = [r, reverse(r)]
  for jj = 0, nregions - 1 do begin
     case jj of
        0: begin
           m  = 0.5 * (evodat.INNFL_MGH[*,1] + evodat.INNFL_MGH[*,1])
           em = 0.5 * (evodat.INNFL_MGH[*,2] - evodat.INNFL_MGH[*,0]) / m / alog(10)
        end
        1: begin
           m = 0.5 * (evodat.INTUP_MGH[*,1] + evodat.INTDN_MGH[*,1])
           em1 = 0.5 * (evodat.INTUP_MGH[*,2] - evodat.INTUP_MGH[*,0]) / m / alog(10)
           em2 = 0.5 * (evodat.INTDN_MGH[*,2] - evodat.INTDN_MGH[*,0]) / m / alog(10)
           em  = sqrt(em1^2 + em2^2) / sqrt(2.)
        end
        2: begin
           m = 0.5 * (evodat.OUTUP_MGH[*,1] + evodat.OUTDN_MGH[*,1])
           em1 = 0.5 * (evodat.OUTUP_MGH[*,2] - evodat.OUTUP_MGH[*,0]) / m / alog(10)
           em2 = 0.5 * (evodat.OUTDN_MGH[*,2] - evodat.OUTDN_MGH[*,0]) / m / alog(10)
           em  = sqrt(em1^2 + em2^2) / sqrt(2.)
        end
     endcase

     masses[jj,*]   = alog10(m[edexes])
     masserrs[jj,*] = em[edexes]
     
  endfor

  plot, [-0.1,2], [-1,0], /xsty, $
        xtitle = '!18r/r!S!X!De!N!R!Eobs!N', $
        ytitle = 'log [!18M!X!D*!N(!18r!X,!18t!X)!18/!XM!D*,tot!N(!18t!X)]', $
        /nodat, $
        charsize = 1.2, charthick = 4, xthick = 3, ythick = 3, $
        pos = [0.7,0.15,0.9,0.9] + [4*xm,0,5*xm,0], /noer
  for ii = 0, nepochs - 1 do begin
     tmass[ii]    = alog10(total(10.^masses[*,ii]))
     tmasserr[ii] = sqrt(total(masserrs[*,ii]^2))
     
     yyy = [masses[*,ii] - sqrt(tmasserr[ii]^2 + masserrs[*,ii]^2), $
            reverse(masses[*,ii] + sqrt(tmasserr[ii]^2 + masserrs[*,ii]^2))] - tmass[ii]
     case ii of
        0: early_yyy = yyy
        1: obs_yyy   = yyy
        2: late_yyy  = yyy
     endcase
  endfor
  polyfill, xxx, early_yyy, col = cgcolor(string(fix(0))), $
            /line_fill, thick = 1, spacing = 0.05, orien = 45
  polyfill, xxx, late_yyy, col = cgcolor(string(fix(2))), $
            /line_fill, thick = 1, spacing = 0.05, orien = -45
  polyfill, xxx, obs_yyy, col = '777777'x
  oplot, r, masses[*,2] - tmass[2], thick = 7
  oplot, r, masses[*,2] - tmass[2], thick = 4, col = cgcolor(string(fix(2)))
  oplot, r, masses[*,0] - tmass[0], thick = 7
  oplot, r, masses[*,0] - tmass[0], thick = 4, col = cgcolor(string(fix(0)))
  oplot, r, masses[*,1] - tmass[1], thick = 7

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
;plot1evoprof_folded, '00900_1_lgn_bestfit.fits', '00900_1_lgn_recons.fits', '00900_2_evostuff.eps'
