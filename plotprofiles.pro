pro plotprofiles
  
  readcol, 'allLgnBestfits.list', files, f = 'A'
  readcol, 'allExpBestfits.list', files2, f = 'A'
  nfiles = n_elements(files)

  cgloadct, 33, ncolors = nfiles, clip = [10,240]

  !P.multi = [0,nfiles/2+1,2]
  window, 0, xsize = 900, ysize = 500
  ar = 500./900
  for ii = 0, nfiles - 1 do begin
     plot, [-2,2], [-1,0], $
           xtitle = 'R/R!De!N', ytitle = 'log M/M!Dinner!N', $
           /nodat
     for jj = 0, 1 do begin
        case jj of
           0: file = files[ii]
           1: file = files2[ii]
        endcase

        data = mrdfits(file, 1)
        opt = where(data.REGION ne 'OPTFL')
        tdata = data[opt]
        
        tdata = tdata[sort(tdata.RNORM)]

        prof = tdata.LMASS[0] - alog10(total(10.^tdata.LMASS[0]))
        err = 0.5 * (tdata.LMASS[2] - tdata.LMASS[1])
        err[where(tdata.REGION ne 'INNFL')] = $
           sqrt(err[where(tdata.REGION ne 'INNFL')]^2 + err[where(tdata.REGION eq 'INNFL')]^2)

        xxx = [tdata[2].RNORM, tdata[3].RNORM, tdata[4].RNORM]
        xxx = [xxx, reverse(xxx)]

        ;fold it
        int = 0.5 * (prof[1] + prof[3])
        out = 0.5 * (prof[0] + prof[4])
        interr = sqrt((err[1]^2 + err[3]^2)/2)
        outerr = sqrt((err[4]^2 + err[4]^2)/2)
        prof2 = [prof[2], int, out]
        err2  = [err[2], interr, outerr]
        yyy = [prof2-err2, reverse(prof2+err2)]
        case jj of
           0: begin
              polyfill, xxx, yyy, col = cgcolor(string(fix(ii)))
             oploterror, tdata.RNORM, prof, err, errcol = cgcolor(string(fix(ii))), col = cgcolor(string(fix(ii)))
           end
           1: begin
              polyfill, xxx, yyy, /line_fill, spacing = 0.05, thick = 1, orien = 45
              oploterror, tdata.RNORM, prof, err, errcol = cgcolor(string(fix(ii))), col = cgcolor(string(fix(ii))), linesty = 2
           end
        endcase
        legend, /bottom, /right, box = 0, $
                strmid(files[ii],0,7), charsize = 1

        im = read_tiff(tdata[0].TIFF_IMAGE)
        s = size(im, /dim)
        cgimage, im, pos = [!X.WINDOW[0], !Y.WINDOW[1]-0.05/ar, $
                            !X.WINDOW[0]+0.05, !Y.WINDOW[1]], /over
        
     endfor
        
  endfor

  k = get_kbrd(1)
  !P.multi = [0,nfiles/2+1,2]
  window, 0, xsize = 900, ysize = 500

  for ii = 0, nfiles - 1 do begin
     plot, [0,2], [-12,-8], $
           xtitle = 'R/R!De!N', ytitle = 'log sSFR [yr!E-1!N]', $
           /nodat
     for jj = 0, 1 do begin
        case jj of
           0: file = files[ii]
           1: file = files2[ii]
        endcase

        data = mrdfits(file, 1)
        opt = where(data.REGION ne 'OPTFL')
        tdata = data[opt]
        
        tdata = tdata[sort(tdata.RNORM)]

        print, strmid(files[ii],0,7), tdata[2].SFR[0]
        
        prof = alog10(tdata.SFR[0]) - tdata.LMASS[0]
;        prof -= prof[2]
        err1 = 0.5 * (tdata.LMASS[2] - tdata.LMASS[1])
        err2 = 0.5/alog(10) * (tdata.SFR[2] - tdata.SFR[1])/tdata.SFR[0]

        ul = where(tdata.SFR[0] lt 2 * 10.^err2, nul, compl = gd)
        prof[ul] = alog10(2 * 10.^err2[ul]) - tdata[ul].LMASS[0]
        err = sqrt(err1^2 + err2^2)
;        err[where(tdata.REGION ne 'INNFL')] = $
;           sqrt(err[where(tdata.REGION ne 'INNFL')]^2 + err[where(tdata.REGION eq 'INNFL')]^2)

        rf = [tdata[2].RNORM, tdata[3].RNORM, tdata[4].RNORM]
        xxx = [rf, reverse(rf)]

        ;fold it
        int = 0.5 * (prof[1] + prof[3])
        out = 0.5 * (prof[0] + prof[4])
        interr = sqrt((err[1]^2 + err[3]^2)/2)
        outerr = sqrt((err[4]^2 + err[4]^2)/2)
        prof2 = [prof[2], int, out]
        err2  = [err[2], interr, outerr]
        yyy = [prof2-err2, reverse(prof2+err2)] > !Y.CRANGE[0]

        case jj of
           0: $;if nul eq 0 then $
              polyfill, xxx, yyy, col = cgcolor(string(fix(ii)))
           1: $;if nul eq 0 then $
              polyfill, xxx, yyy, /line_fill, spacing = 0.05, thick = 1, orien = 45
        endcase
        oplot, !X.crange, replicate(alog10(4./getage(tdata[0].Z[0])) - 9, 2), linesty = 5, col = '777777'x
        legend, /bottom, /right, box = 0, $
                strmid(files[ii],0,7), charsize = 1

        im = read_tiff(tdata[0].TIFF_IMAGE)
        s = size(im, /dim)
        cgimage, im, pos = [!X.WINDOW[0], !Y.WINDOW[0], !X.WINDOW[0]+0.025, !Y.WINDOW[0]+0.025/ar], /over
        
     endfor
        
  endfor

  k = get_kbrd(1)
  
  !P.multi = [0,nfiles/2+1,2]
  window, 0, xsize = 900, ysize = 500

  for ii = 0, nfiles - 1 do begin
     plot, [-2,2], [0,max(data.HA_EW) * 2], $
           xtitle = 'R/R!De!N', ytitle = 'EW(H'+greek('alpha')+' ['+textoidl('\AA')+']', $
           /nodat
     for jj = 0, 1 do begin
        case jj of
           0: file = files[ii]
           1: file = files2[ii]
        endcase

        data = mrdfits(file, 1)
        opt = where(data.REGION ne 'OPTFL')
        tdata = data[opt]
        
        tdata = tdata[sort(tdata.RNORM)]        
        prof = tdata.HA_EW[0]
        err  = 0.5 * (tdata.HA_EW[2] - tdata.HA_EW[1])

        rf = [tdata[2].RNORM, tdata[3].RNORM, tdata[4].RNORM]
        xxx = [rf, reverse(rf)]

        ;fold it
        int = 0.5 * (prof[1] + prof[3])
        out = 0.5 * (prof[0] + prof[4])
        interr = sqrt((err[1]^2 + err[3]^2)/2)
        outerr = sqrt((err[4]^2 + err[4]^2)/2)
        prof2 = [prof[2], int, out]
        err2  = [err[2], interr, outerr]
        yyy = [prof2-err2, reverse(prof2+err2)] > !Y.CRANGE[0]

        case jj of
           0: begin
              oploterror, tdata.RNORM, prof, err, errcol = cgcolor(string(fix(ii))), col = cgcolor(string(fix(ii)))
              polyfill, xxx, yyy, col = cgcolor(string(fix(ii))), /line_fill, spacing = 0.05, thick = 1, orien = -45
           end
           1: begin
              oploterror, tdata.RNORM, prof, err, errcol = cgcolor(string(fix(ii))), col = cgcolor(string(fix(ii))), linesty = 2
              polyfill, xxx, yyy, /line_fill, spacing = 0.05, thick = 1, orien = 45
           end
        endcase
        legend, /bottom, /right, box = 0, $
                strmid(files[ii],0,7), charsize = 1

        im = read_tiff(tdata[0].TIFF_IMAGE)
        s = size(im, /dim)
        cgimage, im, pos = [!X.WINDOW[0], !Y.WINDOW[0], !X.WINDOW[0]+0.025, !Y.WINDOW[0]+0.025/ar], /over
        
     endfor
        
  endfor

  stop
  
  if 0 then begin
  readcol, 'allExpBestfits.list', files, f = 'A'
  nfiles = n_elements(files)

  cgloadct, 33, ncolors = nfiles, clip = [10,240]

  plot, [-2,2], [-0.7,0], $
        xtitle = 'R/R!De!N', ytitle = 'log M/M!Dinner!N', $
        /nodat
  for ii = 0, nfiles - 1 do begin
     data = mrdfits(files[ii], 1)
     opt = where(data.REGION ne 'OPTFL')
     tdata = data[opt]

     tdata = tdata[sort(tdata.RNORM)]
     
     xxx = [tdata.RNORM, reverse(tdata.RNORM)]
     prof = tdata.LMASS[0] - tdata[where(tdata.REGION eq 'INNFL')].LMASS[0]
     err = 0.5 * (tdata.LMASS[2] - tdata.LMASS[1])
     err[where(tdata.REGION ne 'INNFL')] = $
        sqrt(err[where(tdata.REGION ne 'INNFL')]^2 + err[where(tdata.REGION eq 'INNFL')]^2)
     yyy = [prof-err, reverse(prof+err)]
;     polyfill, xxx, yyy, col = cgcolor(string(fix(ii)))
     oploterror, tdata.RNORM, prof, err, errcol = cgcolor(string(fix(ii))), col = cgcolor(string(fix(ii)))
     
  endfor
  endif
  stop
  
end
