pro plot1evoprof, basicFits, evoFits

  obsdat = mrdfits(basicFits, 1)
  evodat = mrdfits(evoFits, 1)

  z    = evodat.INNFL_REDSHIFT
  time = evodat.INNFL_TIME
  tobs = evodat.INNFL_TOBS
  
  edexes = value_locate(z, getredshift([tobs-1,tobs,tobs+1]))
  etot = where(z ge 0.5 and z le 2)

  regions = ['INNFL', 'INTUP', 'INTDN', 'OUTUP', 'OUTDN']
  cols = [255, '00ff00'x, '00aa00'x, 'ffa500'x, 'ff0000'x]
  nregions = n_elements(regions)

  !p.multi = [0,3,0]
  
  plot, [8,12], [-1,3], $
        /nodat
  plotsym, 0, /fill
  for jj = 0, nregions - 1 do begin
     case jj of
        0: begin
           oplot, alog10(evodat.INNFL_MGH[etot,0]), alog10(evodat.INNFL_SFH[etot,0]), col = cols[jj]
           oplot, alog10(evodat.INNFL_MGH[edexes,0]), alog10(evodat.INNFL_SFH[edexes,0]), col = cols[jj], psym = 8
        end
        1: begin
           oplot, alog10(evodat.INTUP_MGH[etot,0]), alog10(evodat.INTUP_SFH[etot,0]), col = cols[jj]
           oplot, alog10(evodat.INTUP_MGH[edexes,0]), alog10(evodat.INTUP_SFH[edexes,0]), col = cols[jj], psym = 8
        end
        2: begin
           oplot, alog10(evodat.INTDN_MGH[etot,0]), alog10(evodat.INTDN_SFH[etot,0]), col = cols[jj]
           oplot, alog10(evodat.INTDN_MGH[edexes,0]), alog10(evodat.INTDN_SFH[edexes,0]), col = cols[jj], psym = 8
        end
        3: begin
           oplot, alog10(evodat.OUTUP_MGH[etot,0]), alog10(evodat.OUTUP_SFH[etot,0]), col = cols[jj]
           oplot, alog10(evodat.OUTUP_MGH[edexes,0]), alog10(evodat.OUTUP_SFH[edexes,0]), col = cols[jj], psym = 8
        end
        4: begin
           oplot, alog10(evodat.OUTDN_MGH[etot,0]), alog10(evodat.OUTDN_SFH[etot,0]), col = cols[jj]
           oplot, alog10(evodat.OUTDN_MGH[edexes,0]), alog10(evodat.OUTDN_SFH[edexes,0]), col = cols[jj], psym = 8
        end
     endcase
          
  endfor

  massback = fltarr(nregions,n_elements(edexes),3)
  for ii = 0, n_elements(edexes) - 1 do $
     massback[*,ii,*] = [evodat.OUTDN_MGH[edexes[ii],*], evodat.INTDN_MGH[edexes[ii],*], $
                         evodat.INNFL_MGH[edexes[ii],*], $
                         evodat.INTUP_MGH[edexes[ii],*], evodat.OUTUP_MGH[edexes[ii],*]]
  
  plot, [-2,2], [-1,0], $
        xtitle = 'R/R!De!N', ytitle = 'log M/M!Dinner!N', $
        /nodat
  tdata = obsdat[where(obsdat.REGION ne 'OPTFL')]
  tdata = tdata[sort(tdata.RNORM)]

  prof = tdata.LMASS[0] - alog10(total(10.^tdata.LMASS[0]))
  err = 0.5 * (tdata.LMASS[2] - tdata.LMASS[1])
  err[where(tdata.REGION ne 'INNFL')] = $
     sqrt(err[where(tdata.REGION ne 'INNFL')]^2 + err[where(tdata.REGION eq 'INNFL')]^2)
  
  oploterror, tdata.RNORM, prof, err
  oplot, tdata.RNORM, alog10(massback[*,0,0] / total(massback[*,0,0])), col = 255
  oplot, tdata.RNORM, alog10(massback[*,-1,0] / total(massback[*,-1,0])), col = 'ffa500'x

  plot, time, alog10(evodat.INNFL_SFH[*,0]), /nodat, $
        yran = [0,3], $
        xtitle = '!18t!X [Gyr]', $
        ytitle = 'log !18SFR!X [M!D'+sunsymbol()+'!N yr!E-1!N]', $
        xran = minmax(time[edexes]), /xsty
  xxx = [time, reverse(time)]  
  for ii = 0, n_elements(regions) - 2 do begin
     case ii of
        0: sfhs = evodat.INNFL_SFH
        1: sfhs = evodat.INTUP_SFH
        2: sfhs = evodat.INTDN_SFH
        3: sfhs = evodat.OUTUP_SFH
        4: sfhs = evodat.OUTDN_SFH
     endcase     
     yyy = alog10([SFHS[*,0], reverse(SFHS[*,2])])
     polyfill, (xxx > !X.CRANGE[0]) < !X.CRANGE[1], yyy > !Y.CRANGE[0], col = cols[ii], $
               /line_fill, thick = 1, spacing = 0.05, orien = ii * 10
     oplot, time, alog10(SFHS[*,1]), thick = 4
     oplot, time, alog10(SFHS[*,1]), col = cols[ii]
  endfor
  oplot, replicate(tobs, 2), !Y.CRANGE, linesty = 5
  
  stop
   
  
end
