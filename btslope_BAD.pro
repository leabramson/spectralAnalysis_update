pro btslope, inlist, output

  readcol, inlist, sumfiles, f= 'A'
  readcol, repstr(inlist, 'Bestfit', 'Recon'), evofiles, f = 'A'   
  nfiles = n_elements(sumfiles)

  bt    = fltarr(nfiles)
  bterr = fltarr(nfiles)

  dmdt  = fltarr(3, nfiles)
  edmdt = fltarr(3, nfiles)
  
  dt = 0.02
  rad = findgen(101) * dt
  rs = fltarr(3, nfiles)
  for ii = 0, nfiles - 1 do begin

     ;; Calculate 1D B/T
     sumdata = mrdfits(sumfiles[ii], 1)
     sumdata = sumdata[where(sumdata.REGION ne 'OPTFL')]

     minn  = sumdata[0].LMASS[0]
     eminn = 0.5 * (sumdata[0].LMASS[2] - sumdata[0].LMASS[1])

     mint   = mean(sumdata[1:2].LMASS)
     emint1 = 0.5 * (sumdata[1].LMASS[2] - sumdata[1].LMASS[1])
     emint2 = 0.5 * (sumdata[2].LMASS[2] - sumdata[2].LMASS[1])
     emint  = sqrt(emint1^2 + emint2^2)
     
     mout   = mean(sumdata[3:4].LMASS)
     emout1 = 0.5 * (sumdata[3].LMASS[2] - sumdata[3].LMASS[1])
     emout2 = 0.5 * (sumdata[4].LMASS[2] - sumdata[4].LMASS[1])
     emout  = sqrt(emout1^2 + emout2^2)

     m  = [minn, mint, mout]
     em = [eminn, emint, emout]
     tr = sumdata[[0,2,4]].RNORM

     rs[*,ii] = tr
     
     t = interpol(m, tr, rad)     

     bt[ii]    = total(10.^t[where(rad le 0.5)]) / total(10.^t)
     bterr[ii] = (10.^sqrt(total(em^2)) - 1) * bt[ii]

     ;; Calculate the mass growth histories and shit
     evodata = mrdfits(evofiles[ii], 1)
     q = where(evodata.REGION eq 'INNFL')
     time    = evodata[q].TIME
     z       = evodata[q].REDSHIFT
     tobs    = evodata[q].TOBS
     epochs  = value_locate(time, tobs+[-1,0])
     
     mghInner0 = evodata[q].MGH[epochs[0],1]
     mghInner1 = evodata[q].MGH[epochs[1],1]
     mghInnErr0 = 0.5 * (evodata[q].MGH[epochs[0],2] - evodata[q].MGH[epochs[0],0])
     mghInnErr1 = 0.5 * (evodata[q].MGH[epochs[1],2] - evodata[q].MGH[epochs[1],0])

     mghInter0 = 0.5 * (evodata.INTUP_MGH[epochs[0],1] + evodata.INTDN_MGH[epochs[0],1])
     mghInter1 = 0.5 * (evodata.INTUP_MGH[epochs[1],1] + evodata.INTDN_MGH[epochs[1],1])
     mghIntErr0 = sqrt($
                  (0.5 * (evodata.INTUP_MGH[epochs[0],2] - evodata.INTUP_MGH[epochs[0],0]))^2 + $
                  (0.5 * (evodata.INTDN_MGH[epochs[0],2] - evodata.INTDN_MGH[epochs[0],0]))^2)
     mghIntErr1 = sqrt($
                  (0.5 * (evodata.INTUP_MGH[epochs[1],2] - evodata.INTUP_MGH[epochs[1],0]))^2 + $
                  (0.5 * (evodata.INTDN_MGH[epochs[1],2] - evodata.INTDN_MGH[epochs[1],0]))^2)

     mghOuter0 = 0.5 * (evodata.OUTUP_MGH[epochs[0],1] + evodata.OUTDN_MGH[epochs[0],1])
     mghOuter1 = 0.5 * (evodata.OUTUP_MGH[epochs[1],1] + evodata.OUTDN_MGH[epochs[1],1])
     mghOutErr0 = sqrt($
                  (0.5 * (evodata.OUTUP_MGH[epochs[0],2] - evodata.OUTUP_MGH[epochs[0],0]))^2 + $
                  (0.5 * (evodata.OUTDN_MGH[epochs[0],2] - evodata.OUTDN_MGH[epochs[0],0]))^2)
     mghOutErr1 = sqrt($
                  (0.5 * (evodata.OUTUP_MGH[epochs[1],2] - evodata.OUTUP_MGH[epochs[1],0]))^2 + $
                  (0.5 * (evodata.OUTDN_MGH[epochs[1],2] - evodata.OUTDN_MGH[epochs[1],0]))^2)

     dmdti  = alog10(mghInner1/mghInner0)
     edmdti = sqrt((1./alog(10) * mghInnErr1/mghInner1)^2 + $
                   (1./alog(10) * mghInnErr0/mghInner0)^2)

     dmdtint  = alog10(mghInter1/mghInter0)
     edmdtint = sqrt((1./alog(10) * mghIntErr1/mghInter1)^2 + $
                     (1./alog(10) * mghIntErr0/mghInter0)^2)

     dmdtout  = alog10(mghOuter1/mghOuter0)
     edmdtout = sqrt((1./alog(10) * mghOutErr1/mghOuter1)^2 + $
                     (1./alog(10) * mghOutErr0/mghOuter0)^2)

     dmdt[*,ii]  = [dmdti, dmdtint, dmdtout]
     edmdt[*,ii] = [edmdti, edmdtint, edmdtout]
     
  endfor

  rp = mean(rs, dim = 2)

  set_plot, 'PS'
  device, filename = output, $
          /col, /encap, /decomp, bits_per_pixel = 8, $
          xsize = 7, ysize= 5, /in
  plotsym, 0, /fill
  plot, bt, dmdt[0,*], /nodat, $
        xtitle = '!18B/T!X (!18r/r!X!De!N<2)', $
        ytitle = 'log d!18M!X!D*!N/d!18t!X [Gyr!e-1!N; 1 Gyr ago to !18t!X!Dobs!N]', $
        xran = [min(bt)-0.1,max(bt)+0.1], yran = [-0.3,1.5], $
        charsize = 1.25, charthick = 4, xthick = 3, ythick = 3
  cols = [255, '00aa00'x, 'ff5500'x]
  for ii = n_elements(cols) - 1, 0, -1 do begin
     oploterror, bt, dmdt[ii,*], bterr, edmdt[ii,*], $
                 psym = 8, $
                 symsize = 1.25, errthick = 3, errcol = '777777'x
     oplot, bt, dmdt[ii,*], $
                 psym = 8, col = cols[ii], $;errcol = cols[ii], $
                 symsize = 1
  endfor
  legend, /top, /left, box = 0, $
          ['!18<r/re>!X=0', ' '+string(rp[1], f = '(F4.2)'), ' '+string(rp[2], f = '(F4.2)')], $
          col = cols, psym= 8, charsize = 1.1, charthick = 3, pspacing = 0.5, spacing = 1.5
  device, /close
  set_plot, 'X'
  spawn, 'gv '+output
  
end
; btslope, 'allLgnBestfits.list', 'btslope.eps'
; btslope, 'study5_lgn_Bestfits.list', 'btslope5.eps'

pro t0tau, inlist, output

  readcol, inlist, sumfiles, f= 'A'
  readcol, repstr(inlist, 'Bestfit', 'Recon'), evofiles, f = 'A'   
  nfiles = n_elements(sumfiles)

  t0    = fltarr(nfiles, 3)
  t0err = fltarr(nfiles, 3)

  tau    = fltarr(nfiles, 3)
  tauerr = fltarr(nfiles, 3)

  z = fltarr(nfiles)
  for ii = 0, nfiles - 1 do begin

     data = mrdfits(sumfiles[ii], 1)

     for jj = 0, 2 do begin
        case jj of
           0: begin
              q = where(data.REGION eq 'INNFL')              
              innt0  = data[q].T0[0]
              innt0e = 0.5 * (data[q].T0[2] - data[q].T0[1])
              inntau  = data[q].TAU[0]
              inntaue = 0.5 * (data[q].TAU[2] - data[q].TAU[1])
              z[ii] = data[q].Z[0]
           end
           1: begin
              q1 = where(data.REGION eq 'INTUP')
              q2 = where(data.REGION eq 'INTDN')
              intt0  = 0.5 * (data[q1].T0[0] + data[q2].T0[0])              
              intt0e = sqrt((0.5 * (data[q1].T0[2] - data[q1].T0[1]))^2 + (0.5 * (data[q2].T0[2] - data[q2].T0[1]))^2)
              inttau  = 0.5 * (data[q1].TAU[0] + data[q2].TAU[0])
              inttaue = sqrt((0.5 * (data[q1].TAU[2] - data[q1].TAU[1]))^2 + (0.5 * (data[q2].TAU[2] - data[q2].TAU[1]))^2)
           end
           2: begin
              q1 = where(data.REGION eq 'OUTUP')
              q2 = where(data.REGION eq 'OUTDN')
              outt0  = 0.5 * (data[q1].T0[0] + data[q2].T0[0])              
              outt0e = sqrt((0.5 * (data[q1].T0[2] - data[q1].T0[1]))^2 + (0.5 * (data[q2].T0[2] - data[q2].T0[1]))^2)
              outtau  = 0.5 * (data[q1].TAU[0] + data[q2].TAU[0])
              outtaue = sqrt((0.5 * (data[q1].TAU[2] - data[q1].TAU[1]))^2 + (0.5 * (data[q2].TAU[2] - data[q2].TAU[1]))^2)
           end
        endcase        
     endfor

     t0[ii,*]     = [innt0,intt0,outt0]
     t0err[ii,*] = [innt0e,intt0e,outt0e]

     tau[ii,*]     = [inntau,inttau,outtau]
     tauerr[ii,*] = [inntaue,inttaue,outtaue]
     
  endfor

  plotsym, 0, /fill
  cgloadct, 33
  !p.multi = [0,ceil(nfiles/2), 2]
  tpeak = alog(getage(2))
  tdat = mrdfits('../mikeGenBasic.fits', 1)
  for ii = 0, nfiles- 1 do begin
     plot, [alog(t0[ii,0]/1d9)], [tau[ii,0]], psym = 1, /nodat, $
           xtitle = 'ln !18T!X!D0!N [Gyr]', ytitle = 'ln '+greek('tau')+' [Gyr]', $
           /iso, xran = [0.5,2.5], yran = [0,2]
     cgtext, !X.CRANGE[0], !Y.CRANGE[1] * 0.9, /data, $
             '!18z!X='+string(z[ii], f = '(F4.2)'), charsize = 1
     oplot, thinit(tdat.T0,4), thinit(tdat.TAU,4), psym = 1, col = '777777'x, thick = 1
     tt0 = findgen((!X.CRANGE[1] - !X.CRANGE[0])/0.05+1)*0.05 + !X.CRANGE[0]
     ttau = sqrt(tt0 - tpeak)
     tp2 = alog(getage(z[ii]))
     ttau2 = sqrt(tt0 - tp2)
     oplot, tt0, ttau
     oplot, tt0, ttau2, linesty = 5 
     oploterror, [alog(t0[ii,0]/1d9)], [tau[ii,0]], $
                 0.5/alog(10) * t0err[ii,0]/t0[ii,0], tauerr[ii,0], $
                 psym = 8, errcol = 255, errthick = 1, symsize = 1.1
     oplot, [alog(t0[ii,0]/1d9)], [tau[ii,0]], col = 255, psym = 8
     oploterror, [alog(t0[ii,1]/1d9)], [tau[ii,1]], $
                 0.5/alog(10) * t0err[ii,1]/t0[ii,1], tauerr[ii,1], $
                 psym = 8, errcol = '00ff00'x, errthick = 1, symsize = 1.1
     oplot, [alog(t0[ii,1]/1d9)], [tau[ii,1]], psym = 8, col = '00ff00'x
     oploterror, [alog(t0[ii,2]/1d9)], [tau[ii,2]], $
                 0.5/alog(10) * t0err[ii,2]/t0[ii,2], tauerr[ii,2], $
                 psym = 8, errcol = 'ff0000'x, errthick = 1, symsize = 1.1
     oplot, [alog(t0[ii,2]/1d9)], [tau[ii,2]], psym = 8, col = 'ff0000'x
  endfor
    
end
; t0tau, 'allLgnBestfits.list', 't0tauslope.eps'
; t0tau, 'study3_lgn_Bestfits.list', 't0tauslope3.eps'
