pro mosaic, inlist, output

  readcol, inlist, files, f= 'A'
  nfiles = n_elements(files)

;  multiplot, /init
;  multiplot, [ceil(nfiles/2), 2], $
;             xgap = 0, ygap = 0
  xm = 0.05
  ym = 0.05
  if nfiles gt 5 then begin
     nx = 5.
     ny = 2.
  endif else begin
     nx = nfiles
     ny = 1.
  endelse
  xw = (1 - 2*xm) / nx
  xs = xm + findgen(nx) * xw

  set_plot, 'PS'
  if nfiles gt 5 then $
     xsize = 8.5 $
  else $
     xsize = 8.5 * nfiles/5.
  ysize = xsize/2.
  
  device, filename = output, $
          xsize = xsize, ysize = ysize, /in, $
          /col, /encap, /decomp, bits_per_pix = 8
  ar = ysize/xsize
  marg = [0.025,0.025/ar, 0.025, 0.025/ar]
  for ii = 0, nfiles - 1 do begin
     tdata  = mrdfits(files[ii], 1)
     im     = read_tiff(tdata[0].TIFF_IMAGE)
     tim = im
     for jj = 0, 2 do $
        tim[jj,*,*] = rotate(reform(im[jj,*,*]), 7)
     im = tim
     
     foo    = mrdfits(tdata[0].FITSFILE, 1)
     foo2   = mrdfits(foo.FILENAME, 'sci', head, /silent)
     pscale = sxpar(head, 'CD2_2')

     hlrObs = foo.RE ;; pix
     zpick  = tdata[0].Z[0]     
     kpc    = 1. / zang(1.0, zpick, /silent) ;; kpc / asec
     kpcPix = kpc * pscale                   ;; kpc / pix
     hlrKpc = hlrObs * kpcPix

     s  = size(im, /dim)
     dx = s[1]
     dy = s[2]

     width = median(foo.outerup - foo.outerdn)
     outer = median(foo.outerup - foo.interUp)
     inter = median(foo.interup - foo.innerUp)
     inner = median(foo.innerup - foo.innerDn)

     yctr = s[2]/2
     x = replicate(55,2)

     stats = '!18z!X=!N'+string(zpick, f = '(F4.2)')+$
             '!C!18r!X!De!N='+string(hlrObs*pscale, f = '(F3.1)')+'"='+$
             string(hlrkpc, f = '(F3.1)')+' kpc'

     if nfiles gt 5 then begin
        if ii le nfiles/2-1 then begin
           ys = ym+xw/ar
           jj = ii
        endif else begin
           ys = ym
           jj = ii-nfiles/2
        endelse
     endif else begin
        ys = ym
        jj = ii
     endelse
     
     plot, findgen(dx), findgen(dy), $
           xtickname = replicate(' ', 60), $
           ytickname = replicate(' ', 60), /xsty, /ysty, $
           pos = [xs[jj],ys,xs[jj]+xw,ys+xw/ar]+marg, /noer
     cgimage, im, /over
     oplot, replicate(10,2), 5+[0, hlrObs], thick = 8, col = 'ffffff'x
     oplot, replicate(15,2), 5+[0, 5/kpcPix], thick = 8, col = 'ffffff'x
     oplot, x, yctr+width*[-0.5,0.5], thick = 8, col = 'ffffff'x
     oplot, x, yctr-width/2.+[0,outer], thick = 8, col = 'ff0000'x
     oplot, x, yctr+width/2.-[0,outer], thick = 8, col = 'ffa500'x
     oplot, x, yctr-width/2.+outer+[0,inter], thick = 8, col = '00aa00'x
     oplot, x, yctr+width/2.-outer-[0,inter], thick = 8, col = '00ff00'x
     oplot, x, yctr-width/2.+outer+inter+[0,inner], thick = 8, col = 255
     if ii eq 0 then begin
        cgtext, 6, 2+hlrObs/2., $
                '!18r!X!De!N', align = 0, charsize = 0.9, /data, orien = 90, $
                charthick = 3, col = 'ffffff'x
        cgtext, 22, 3, '5 kpc', align = 0, charsize = 0.9, /data, orien = 90, $
                charthick = 3, col = 'ffffff'x
     endif
     cgtext, !X.CRANGE[0]+2, !Y.CRANGE[1]-7, /data, $
             stats, charsize = 0.9, charthick = 3, col = 'ffffff'x
     cgtext, !X.CRANGE[1]-1, !Y.CRANGE[0]+5, /data, $
             strmid(files[ii],0,7), charsize = 0.8, charthick = 3, col = 'ffffff'x, align = 1
     plot, findgen(dx) * pscale, findgen(dy) * pscale, $
           xtickname = replicate(' ', 60), $
           ytickname = replicate(' ', 60), $
           /xsty, /ysty, xminor = 2, yminor = 2, /nodat, $
           pos = [xs[jj],ys,xs[jj]+xw,ys+xw/ar]+marg, /noer, col = 'ffffff'x, $
           xthick = 3, ythick = 3

     if ii eq 0 then begin
        xxs = !X.WINDOW[0]
        yys = !Y.WINDOW[1]
     endif
     if ii eq nfiles - 1 then begin
        xxst = !X.WINDOW[1]
        yyst = !Y.WINDOW[0]
     endif
;     multiplot
     
  endfor

  cgtext, xm+nx/2.*xw+0.025, xm/2, greek('delta')+'!18x!X!N ["]', $
          charsize = 1.25, charthick = 4, align = 0.5, /norm
  if nfiles gt 5 then begin
     cgtext, xm/2, xm+(xw+0.025)/ar, greek('delta')+'!18y!X!N ["]', $
             charsize = 1.25, charthick = 4, align = 0.5, /norm, orien = 90
  endif else begin
     cgtext, xm/2, mean(!Y.WINDOW), greek('delta')+'!18y!X!N ["]', $
             charsize = 1.25, charthick = 4, align = 0.5, /norm, orien = 90
  endelse
     
  device, /close
  set_plot, 'X'
  spawn, 'gv '+output+' &'
;  multiplot, restore
;  cleanplot

end
;mosaic, 'allExpBestfits.list', 'mosaic.eps'
;mosaic, 'study3_lgn_Bestfits.list', 'mosaic3.eps'
;mosaic, 'study5_lgn_Bestfits.list', 'mosaic4.eps'
