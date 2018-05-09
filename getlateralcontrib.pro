;;
;; Calculate how much flux in the INNER
;; extraction region is coming from the
;; INNER, INTER, OUTER, and total regions
;; rotated by 90 degrees. Attempts to
;; constrain what fraction of the INNER
;; signal is not r~0 but from larger r
;; regions
;;

function getlateralcontrib, infile

  data = mrdfits(infile, 1)

  ;; Reconstruct the INNER LSF.
  ;; Use this as the full spatial profile
  ;; in this region.

  lsfpar = data.LSF_INNER_PARS
  s = size(data.DIRECT_IM, /dim)
  x = 2 * s[0]
  dx = 1./3.
  trun = x/dx + 1
  run = findgen(trun) * dx - s[0]

  u       = run / lsfpar[2] 
  profile = lsfpar[0] / (u^2 + 1)^lsfpar[3] 
  profile /= total(profile) 

  ;; cumulative profile
  cprof = total(profile, /cum)

  ;; get regional fractions
  print, [data.NPIX_INNER[0], data.NPIX_INTER_UP[0], data.NPIX_OUTER_UP[0]]
  innSize = data.NPIX_INNER[0] / dx
  intSize = data.NPIX_INTER_UP[0] / dx 
  outSize = data.NPIX_OUTER_UP[0] / dx

  inner = n_elements(run)/2 + [-innSize/2., innSize/2.]; - 1/dx
  inter = inner[1] + [0,intSize]
  outer = inter[1] + [0,outSize]
  onstamp = where(abs(run) le s[0]/2)
  
  hit = onstamp[value_locate(cprof[onstamp], [0.25,0.75])]
  
  plot, run, cprof, xran = [-15,15]
  oplot, replicate(run[n_elements(run)/2], 2), !Y.CRANGE
  oplot, !X.CRANGE, [0.5,0.5]
  oplot, run[onstamp], cprof[onstamp], col = '777777'x, thick = 3
  oplot, run[inner[0]:inner[1]], cprof[inner[0]:inner[1]], $
         col = '0000ff'x, thick = 4
  oplot, run[inter[0]:inter[1]], cprof[inter[0]:inter[1]], $
         col = '00ff00'x, thick = 4
  oplot, run[outer[0]:outer[1]], cprof[outer[0]:outer[1]], $
         col = 'ff0000'x, thick = 4
;  oplot, replicate(run[hit[0]+1], 2), !Y.CRANGE
;  oplot, replicate(run[hit[1]+1], 2), !Y.CRANGE
;  oplot, 2*replicate(run[hit[0]+1], 2), !Y.CRANGE, linesty = 5
;  oplot, 2*replicate(run[hit[1]+1], 2), !Y.CRANGE, linesty = 5
;  oplot, 3*replicate(run[hit[0]+1], 2), !Y.CRANGE, linesty = 2
;  oplot, 3*replicate(run[hit[1]+1], 2), !Y.CRANGE, linesty = 2
;  oplot, 4*replicate(run[hit[0]+1], 2), !Y.CRANGE, linesty = 1
;  oplot, 4*replicate(run[hit[1]+1], 2), !Y.CRANGE, linesty = 1
  stmpFrac = cprof[onstamp[-1]] - cprof[onstamp[0]]
  innFrac = cprof[inner[-1]] - cprof[inner[0]]
  intFrac = (cprof[inter[-1]] - cprof[inter[0]]) * 2
  outFrac = (cprof[outer[-1]] - cprof[outer[0]]) * 2
  totFrac = total([innFrac,intFrac,outFrac])

;  stop
  
  output = {RUN: run, PROF: profile, CPROF: cprof, $
            INNER: inner, INTER: inter, OUTER: outer, $
            STMPSIZE: s[0], STMPFRAC: stmpFrac, TOTFRAC: totFrac, $
            INNFRAC: innfrac, INTFRAC: intfrac, OUTFRAC: outFrac}
  
  RETURN, output
end

;;

pro plotstuff

  files = ['../MACS2129/00451_2_B.fits', $
           '../MACS0744/00660_2_B.fits', $
           '../MACS1149/00900_1_B.fits', $
           '../MACS1423/01916_2_B.fits'] 
  mfiles = ['../MACS2129/00451_2_B_allbandImages.fits', $
            '../MACS0744/00660_2_B_allbandImages.fits',$
            '../MACS1149/00900_1_B_allbandImages.fits',$
            '../MACS1423/01916_2_B_allbandImages.fits']
  clusters = 'MACS'+['2129','0744','1149','2129']
  names = ['CSF', 'PAS', 'SSF', 'PSB']
  galstys = [0,0,0,0];[0,1,2,3]
  
  nfiles = n_elements(files)

  ;; LSF, 1D INNER, 2D INNER
  innFs = fltarr(nfiles,3)
  intFs = fltarr(nfiles,3)
  outFs = fltarr(nfiles,3)
  
  totFs = fltarr(nfiles,3)
  stmpFs = fltarr(nfiles)

  fluxes = fltarr(4, nfiles)
  cfluxes = fluxes
  
  snrs = []
  for ii = 0, nfiles - 1 do begin
     foo = getlateralcontrib(files[ii])
     foo2 = getstamplateralcontrib(master = files[ii], $
                                   multiband = mfiles[ii], $
                                   cluster = clusters[ii]);, shift = 1)

     innFs[ii,0] = foo.INNFRAC
     intFs[ii,0] = foo.INTFRAC
     outFs[ii,0] = foo.OUTFRAC

     innFs[ii,1:*] = [foo2.PROF_INN_INN, foo2.INN_INN] / [foo2.PROF_TOT, foo2.INN_FLUX]
     intFs[ii,1:*] = [foo2.PROF_INT_INN, foo2.INT_INN] / [foo2.PROF_TOT, foo2.INN_FLUX]
     outFs[ii,1:*] = [foo2.PROF_OUT_INN, foo2.OUT_INN] / [foo2.PROF_TOT, foo2.INN_FLUX]
          
     totFs[ii,0] = foo.TOTFRAC
     totFs[ii,1:*] = [total([innFs[ii,1],intFs[ii,1],outFs[ii,1]]), $
                      total([innFs[ii,2],intFs[ii,2],outFs[ii,2]])]

;     fluxes[*,ii] = [foo2.INN_INN, foo2.INT_INN, foo2.OUT_INN, foo2.INN_FLUX]
;     cfluxes[*,ii] = [foo2.PROF_INN_INN, foo2.PROF_INT_INN, foo2.PROF_OUT_INN, foo2.PROF_TOT]
     
     stmpFs[ii] = foo.STMPFRAC

     
;     if ii eq 2 then stop
     
     foo = mrdfits(files[ii], 1, /silent)
     foor = mrdfits(repstr(files[ii], "_B", "_R"), 1, /silent)
    
     u = where(~foo.MASK_INNER AND $
               foo.LAMBDA ge 8400 AND foo.LAMBDA le 1.12d4)
     ur = where(~foor.MASK_INNER AND $
                foor.LAMBDA ge 1.1d4 AND foor.LAMBDA le 1.65d4)

     sn = median([foo.F_INNER[u] / sqrt(foo.VAR_INNER[u]), $
                  foor.F_INNER[ur] / sqrt(foor.VAR_INNER[ur])])
     snrs = [snrs, sn]
  endfor

  rs = findgen(6) + 1
;  as = transpose([[innFs], [intFs], [outFs], [totFs], [stmpFs]])

  tcols = [255, '00a500'x, 'ff0000'x, '777777'x, 0]
  plotsym, 0, /fill

  reg = [' ','INNER', 'INTER', 'OUTER', $
         'Tot.', 'OnStamp',' ']

  hit = (where(names eq 'SSF'))[0]
  print, ''
  print, 'SSF'
  print, [innFs[hit,0], intFs[hit,0], outFs[hit,0]]
  print, innFs[hit,0] / [innFs[hit,0], intFs[hit,0], outFs[hit,0]]
  print, [innFs[hit,2], intFs[hit,2], outFs[hit,2]]
  print, innFs[hit,2] / [innFs[hit,2], intFs[hit,2], outFs[hit,2]]
  hit = (where(names eq 'PSB'))[0]
  print, ''
  print, 'PSB'
  print, [innFs[hit,0], intFs[hit,0], outFs[hit,0]]
  print, innFs[hit,0] / [innFs[hit,0], intFs[hit,0], outFs[hit,0]]
  print, [innFs[hit,2], intFs[hit,2], outFs[hit,2]]
  print, innFs[hit,2] / [innFs[hit,2], intFs[hit,2], outFs[hit,2]]
;  stop
  
;  print, names
;  print, snrs

;  names += ',!18<S/N>!X!DINNER!N='+string(snrs, f = '(f4.1)')

  set_plot, 'PS'
  device, filename = 'INNERfluxFrac.eps', $
          /col, /encap, /decomp, bits_per_pix = 8, $
          xsize = 5, ysize = 4, /in
  plot, rs, rs, psym = 2, /nodat, $
        xran = [0,6], yran = [0,1.05], /ysty, $
;        xtitle = 'Zone', $
        ytitle = 'INNER flux fractional contribution', $
        xminor = 1, xtickname = reg, $
        xticks = n_elements(rs), $
        xthick = 4, ythick = 4, $
        charthick = 4, charsize = 1.25, $
        pos = [0.15,0.15,0.95,0.95]

  ;; Cartoon
  r = 0.3 & rscale = 5
  xc = 0.9 & yc = 0.875
  xoff = 1.2
  theta = findgen(16+1) * 2.*!pi/16.
  cgcolorfill, xc + r * cos(theta), $
               yc + r/rscale * sin(theta), col = 'ff5500'x, $
               /line_fill, orien = 45, thick = 1, spacing = 0.05
  cgcolorfill, xc + 0.5 * r * cos(theta), $
               yc + 0.5 * r/rscale * sin(theta), col = '00ffff'x, $
               /line_fill, orien = -45, thick = 2, spacing = 0.025
;  oplot, r * cos(theta) + xc, r/rscale * sin(theta) + yc, $
;         thick = 4, col = '777777'x
  polyfill, xc + 0.3 * [-1,-1,1,1], yc + r/rscale/7. * [-1,1,1,-1], col = 255
  oplot, xc + 0.3 * [-1,1], yc + r/rscale/2.5 * [1,1], col = '00a500'x, thick = 2
  oplot, xc + 0.3 * [-1,1], yc - r/rscale/2.5 * [1,1], col = '00a500'x, thick = 2
  oplot, xc + 0.3 * [-1,1], yc + r/rscale/1.5 * [1,1], col = 'ff0000'x, thick = 2
  oplot, xc + 0.3 * [-1,1], yc - r/rscale/1.5 * [1,1], col = 'ff0000'x, thick = 2
  arrow, xc + 0.4, yc, xc + 0.8, yc, hsize = -0.2, thick = 6, /data, col = '777777'x

  cgcolorfill, xc + xoff+ r * cos(theta), $
               yc + r/rscale * sin(theta), col = 'ff5500'x, $
               /line_fill, orien = 45, thick = 1, spacing = 0.05
  cgcolorfill, xc + xoff + 0.5 * r * cos(theta), $
               yc + 0.5 * r/rscale * sin(theta), col = '00ffff'x, $
               /line_fill, orien = -45, thick = 2, spacing = 0.025
;  oplot, r * cos(theta) + xoff + xc, r/rscale * sin(theta) + yc, $
;         thick = 4, col = '777777'x
  polyfill, xoff + xc + 0.3 * [-1,-1,1,1], yc + r/rscale/7. * [-1,1,1,-1], col = 'aaaaaa'x
  polyfill, xoff + xc + r/7 * [-1,-1,1,1], yc + r/rscale/7. * [-1,1,1,-1], col = 255
  polyfill, xoff + xc + r/[-2.5,-2.5,-7,-7], yc + r/rscale/7. * [-1,1,1,-1], col = '00a500'x
  polyfill, xoff + xc + r/[2.5,2.5,7,7], yc + r/rscale/7. * [-1,1,1,-1], col = '00a500'x
  polyfill, xoff + xc + r/[-1.5,-1.5,-2.5,-2.5], yc + r/rscale/7. * [-1,1,1,-1], col = 'ff0000'x
  polyfill, xoff + xc + r/[1.5,1.5,2.5,2.5], yc + r/rscale/7. * [-1,1,1,-1], col = 'ff0000'x
  oplot, xoff + xc + r/7.  * [1,1], yc + r/rscale * [-1,1], col = 255, thick = 2
  oplot, xoff + xc - r/7.  * [1,1], yc + r/rscale * [-1,1],col = 255, thick = 2
  oplot, xoff + xc + r/2.5 * [1,1], yc + r/rscale * [-1,1], col = '00a500'x, thick = 2
  oplot, xoff + xc - r/2.5 * [1,1], yc + r/rscale * [-1,1], col = '00a500'x, thick = 2
  oplot, xoff + xc + r/1.5 * [1,1], yc + r/rscale * [-1,1], col = 'ff0000'x, thick = 2
  oplot, xoff + xc - r/1.5 * [1,1], yc + r/rscale * [-1,1], col = 'ff0000'x, thick = 2

  ;; Useful  thresholds
  oplot, !X.CRANGE, replicate(0.5, 2), thick = 10, $
         col = '777777'x, linesty = 5
  oplot, !X.CRANGE, replicate(1.0, 2), thick = 10, $
         col = '777777'x, linesty = 5

  ;; Division between regions and global fractions
  oplot, replicate(3.5,2), !Y.CRANGE, thick = 10

  cgloadct, 33, /rev, ncolors = 9, clip = [10,240]
  cols = cgcolor(['4','0','5','2'])

  roffs = [-0.2,-0.1,0.1,0.2] / 2.
  
  ;; Actual points & labels
  for ii = 0, nfiles - 1 do begin
     cgtext, rs[0]-0.3, mean(innFs[ii,[0,2]], dim = 2), /data, $
             string(snrs[ii], f = '(f4.1)'), align = 1, $
             charsize = 1, charthick = 4
     oplot, rs[0:2] + roffs[ii], [innFs[ii,0], intFs[ii,0], outFs[ii,0]], $
            thick = 10, linesty = 2                                                   ;, linesty = galstys[ii], thick = 8
     oplot, rs[0:2] + roffs[ii], [innFs[ii,0], intFs[ii,0], outFs[ii,0]], $
            col = cols[ii], thick = 5, linesty = 2                                                   ;, linesty = galstys[ii], thick = 8
     oplot, rs[0:2] + roffs[ii], [innFs[ii,2], intFs[ii,2], outFs[ii,2]], $
            thick = 10  
     oplot, rs[0:2] + roffs[ii], [innFs[ii,2], intFs[ii,2], outFs[ii,2]], $
            col = cols[ii], thick = 5      ;, linesty = galstys[ii], thick = 8

     oplot, rs[-3:-2] + roffs[ii], [totFs[ii,0],stmpFs[ii]], thick = 10, linesty = 2
     oplot, rs[-3:-2] + roffs[ii], [totFs[ii,2],stmpFs[ii]], thick = 10
     oplot, rs[-3:-2] + roffs[ii], [totFs[ii,0],stmpFs[ii]], thick = 5, col = cols[ii], linesty = 2
     oplot, rs[-3:-2] + roffs[ii], [totFs[ii,2],stmpFs[ii]], thick = 5, col = cols[ii]
     for jj = 0, 2 do begin
        oplot, [rs[jj]] + roffs[ii], [([innFs[ii,2], intFs[ii,2], outFs[ii,2]])[jj]], $
               psym = 8, symsize = 2
        oplot, [rs[jj]] + roffs[ii], [([innFs[ii,2], intFs[ii,2], outFs[ii,2]])[jj]], $
               psym = 8, symsize = 1.5, col = cols[ii]
     endfor
     oplot, [rs[3]] + roffs[ii], [totFs[ii,2]], $
            psym = 8, symsize = 2
     oplot, [rs[3]] + roffs[ii], [totFs[ii,2]], $
            psym = 8, symsize = 1.5, col = cols[ii]
     oplot, [rs[4]] + roffs[ii], [stmpFs[ii]], $
            psym = 8, symsize = 2
     oplot, [rs[4]] + roffs[ii], [stmpFs[ii]], $
            psym = 8, symsize = 1.5, col = cols[ii]
  endfor
  cgtext, rs[0]-0.3, min(innFs[*,*])-0.1, /data, $
          '!18<S/N>!X pix!E-1!N', align = 0.5, $
          charsize = 1, charthick = 4

  ;; Legend
  legend, /bottom, /right, box = 0, $
          [names, 'LSF', 'Image'], $
          linesty = [galstys,2,0], pspacing = 1.5, $
          psym = [0,0,0,0,0,-8], $
          charsize = 1.15, charthick = 4, thick = 6, col = [cols, 0, 0]

  device, /close
  set_plot, 'X'
  spawn, 'gv INNERfluxFrac.eps &'
  

end
