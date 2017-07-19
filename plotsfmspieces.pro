pro plotsfmspieces, indata, output

  readcol, indata, dirs, f = 'A'
  nobj = n_elements(dirs)

  regions  = ['INNER', 'INTER', 'OUTER', 'OPTFL']
  nregions = n_elements(regions)
  
  sfr   = fltarr(nregions, nobj, 3)
  mstel = fltarr(nregions, nobj, 3)

  !P.MULTI = [0,2,ceil(nobj/2.)]
if 0 then begin
  window, 2, xsize = 800, ysize = 1000, retain = 2
  
  plotsym, 0, /fill
  psyms = [1,2,4,5,6,8]
  for ii = 0, nobj - 1 do begin
     analyzetofits, dirs[ii], 'tmp.fits'     
     tdata = mrdfits('tmp.fits', 1)
     if ii eq 0 then $
        regions = tdata.REGION
     ha = 2.8 * tdata.HA_FLUX
     m  = tdata.LMASS 

     plot, findgen(10), findgen(4) - 2, /nodat, $
           xran = [9.5,11.5], yran = [-1,2.5], $
           xtitle = 'log Mstel', ytitle = 'log SFR'
     oplot, !X.CRANGE, !Y.CRANGE, linesty = 1, thick = 1
     for jj = 0, 5 do begin
        sfr[jj,ii,*] = alog10([gethasfr(ha[0,jj],tdata[3].Z[0]), $
                               gethasfr(ha[1,jj],tdata[3].Z[0]), $
                               gethasfr(ha[2,jj],tdata[3].Z[0])])

        case regions[jj] of
           'INNFL': col = 255      
           'INTER': col = '00a500'x
           'OUTER': col = 'ffa500'x
;           'OUTUP': col = 'ffa500'x
;           'OUTDN': col = 'ff0000'x
;           'OPTFL': col = '777777'x
        endcase

        if sfr[jj,ii,0] gt !Y.CRANGE[0] then $
           oploterror, [m[0,jj]], [sfr[jj,ii,0]], 0.5 * (sfr[jj,ii,2] - sfr[jj,ii,1]), $
                       col = long(col), errcol = long(col), psym = 8, symsize = 2 $
        else $
           oplot, m[0,jj] + [-0.3,0.3], !Y.CRANGE[0] + [0.1,0.1], col = long(col)
        legend, /bottom, /right, box = 0, $
                [repstr(dirs[ii], '_pyspecfitPhotResults', ''), $
                 'z='+string(tdata[3].Z[0], f = '(F5.3)')]
        
     endfor
     
  endfor
  plot, findgen(10), findgen(4) - 2, /nodat, $
        xran = [9.5,11.5], yran = [-1,2], $
        xtickname = replicate(' ',60), ytickname = replicate(' ',60), $
        xsty = 4, ysty = 4
  legend, /bottom, /right, $
          ['R = 0', 'R=0.5 Re', '', 'R=Re', '', 'Total'], $
          col = [255, '00ff00'x, '00aa00'x, 'ffa500'x, 'ff0000'x, '777777'x], $
          psym = 8, pspacing = 1, box = 0

  write_jpeg, 'testResSFMS.jpeg', tvrd(true = 1), true = 1
endif

  age = sfr
  
  window, 2, xsize = 800, ysize = 1000, retain = 2
  
  plotsym, 0, /fill
  psyms = [1,2,4,5,6,8]
  for ii = 0, nobj - 1 do begin
     analyzetofits, dirs[ii], 'tmp.fits'     
     tdata = mrdfits('tmp.fits', 1)
;     if ii eq 0 then $
;        regions = tdata.REGION
     tage = tdata.AGE

     plot, findgen(10), findgen(4) - 2, /nodat, $
           xran = [9.5,11.5], yran = [6.8,alog10(getage(1)) + 9], /ysty, $
           xtitle = 'log Mstel', ytitle = 'log Age'
     
     for jj = 0, nregions - 1 do begin
        age[jj,ii,*] = [tage[0,jj], $
                        tage[1,jj], $
                        tage[2,jj]]
        m  = tdata.LMASS
        
        case regions[jj] of
           'INNFL': col = 255      
           'INTUP': col = '00ff00'x
           'INTDN': col = '00aa00'x
           'OUTUP': col = 'ffa500'x
           'OUTDN': col = 'ff0000'x
           'OPTFL': col = '777777'x
        endcase
        
        if alog10(age[jj,ii,0]) gt !Y.CRANGE[0] then $
           oploterror, [m[0,jj]], [alog10(age[jj,ii,0])], 0.5/alog(10) * (age[jj,ii,2] - age[jj,ii,1]) / age[jj,ii,0], $
                       col = long(col), errcol = long(col), psym = 8, symsize = 2 $
        else $
           oplot, m[0,jj] + [-0.3,0.3], !Y.CRANGE[0] + [0.1,0.1], col = long(col)
        legend, /bottom, /right, box = 0, $
                [repstr(dirs[ii], '_pyspecfitPhotResults', ''), $
                 'z='+string(tdata[3].Z[0], f = '(F5.3)')]
        
     endfor
     
  endfor
  plot, findgen(10), findgen(4) - 2, /nodat, $
        xran = !X.crange, yran = !Y.CRANGE, $
        xtickname = replicate(' ',60), ytickname = replicate(' ',60), $
        xsty = 4, ysty = 4
  legend, /bottom, /right, $
          ['R = 0', 'R=0.5 Re', '', 'R=Re', '', 'Total'], $
          col = [255, '00ff00'x, '00aa00'x, 'ffa500'x, 'ff0000'x, '777777'x], $
          psym = 8, pspacing = 1, box = 0

  write_jpeg, 'testResAge.jpeg', tvrd(true = 1), true = 1
end
;plotsfmspieces, 'neatObjs.list'
