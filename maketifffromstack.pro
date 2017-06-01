pro maketifffromstack, imstack, outname, $
                       min = min, $
                       max = max, $
                       scale = scale, $
                       ds9 = ds9, $
                       RGB = rgb
  
  if NOT keyword_set(min) then min = replicate(0, 3)
  if NOT keyword_set(max) then max = replicate(1, 3)
  if N_ELEMENTS(min) eq 1 then min = replicate(min, 3)
  if N_ELEMENTS(max) eq 1 then max = replicate(max, 3)
  if NOT keyword_set(scale) then scale = 'log'
  if keyword_set(DS9) then begin
     spawn, 'ds9 &'
     wait, 1
  endif
     
  data = mrdfits(imstack, 1, /silent)

  if NOT keyword_set(rgb) then begin

     red = where(data.FILTER le 160, nred)
     grn = where(data.FILTER ge 606, ngrn)
     blu = where(data.FILTER gt 300 AND data.FILTER lt 606, nblu)

     if nred gt 1 then $
        r = median(data[red].IMAGE, dim = 3, /even) $
     else $
        r = data[red].IMAGE
     if ngrn gt 1 then $
        g = median(data[grn].IMAGE, dim = 3)$
     else $
        g = data[grn].IMAGE
     if nblu gt 1 then $
        b = median(data[blu].IMAGE, dim = 3) $
     else $
        b = data[blu].IMAGE
     
  endif else begin

     red = where(data.FILTER eq rgb[0])
     grn = where(data.FILTER eq rgb[1])
     blu = where(data.FILTER eq rgb[2])
          
     r = data[red].IMAGE
     g = data[grn].IMAGE
     b = data[blu].IMAGE
     
  endelse
  
  mwrfits, r, 'tmpRed.fits', /create
  mwrfits, g, 'tmpGrn.fits', /create
  mwrfits, b, 'tmpBlu.fits', /create
  
  spawn, 'xpaset -p ds9 rgb'
  spawn, 'xpaset -p ds9 rgb red'
  spawn, 'xpaset -p ds9 file tmpRed.fits'
  spawn, 'xpaset -p ds9 scale limits '+string(min[0])+' '+string(max[0])
  spawn, 'xpaset -p ds9 scale '+scale
  spawn, 'xpaset -p ds9 rgb green'
  spawn, 'xpaset -p ds9 file tmpGrn.fits'
  spawn, 'xpaset -p ds9 scale limits '+string(min[1])+' '+string(max[1])
  spawn, 'xpaset -p ds9 scale '+scale
  spawn, 'xpaset -p ds9 rgb blue'
  spawn, 'xpaset -p ds9 file tmpBlu.fits'
  spawn, 'xpaset -p ds9 scale limits '+string(min[2])+' '+string(max[2])
  spawn, 'xpaset -p ds9 scale '+scale
  spawn, 'xpaset -p ds9 export tiff '+outname

  im = read_tiff(outname)
  s = size(im, /dim)

  plot, findgen(s[1]), findgen(s[2]), /iso, /nodat, $
        title = imstack, xtitle = 'pix', ytitle = 'pix'
  cgimage, im, /over, stretch = 2
  
end
