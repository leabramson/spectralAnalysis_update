pro maketifffromstack, imstack, outname, $
                       min = min, $
                       max = max, $
                       scale = scale, $
                       ds9 = ds9, $
                       RGB = rgb

  if NOT keyword_set(min) then min = 0
  if NOT keyword_set(max) then max = 1
  if NOT keyword_set(scale) then scale = 'log'

  if keyword_set(DS9) then begin
     spawn, 'ds9 &'
     wait, 1
  endif
     
  data = mrdfits(imstack, 1, /silent)

  if NOT keyword_set(rgb) then begin
     red = where(data.FILTER le 160)
     grn = where(data.FILTER ge 606)
     blu = where(data.FILTER gt 300 AND data.FILTER lt 606)
  endif else begin
     red = where(data.FILTER eq rgb[0])
     red = where(data.FILTER eq rgb[1])
     red = where(data.FILTER eq rgb[2])
  endelse
     
  r = median(data[red].IMAGE, dim = 3);, /nan)
  g = median(data[grn].IMAGE, dim = 3);, /nan)
  b = median(data[blu].IMAGE, dim = 3);, /nan)
  
  mwrfits, r, 'tmpRed.fits', /create
  mwrfits, g, 'tmpGrn.fits', /create
  mwrfits, b, 'tmpBlu.fits', /create
  
  spawn, 'xpaset -p ds9 rgb'
  spawn, 'xpaset -p ds9 rgb red'
  spawn, 'xpaset -p ds9 file tmpRed.fits'
  spawn, 'xpaset -p ds9 scale limits '+string(min)+' '+string(max)
  spawn, 'xpaset -p ds9 scale '+scale
  spawn, 'xpaset -p ds9 rgb green'
  spawn, 'xpaset -p ds9 file tmpGrn.fits'
  spawn, 'xpaset -p ds9 scale limits '+string(min)+' '+string(max)
  spawn, 'xpaset -p ds9 scale '+scale
  spawn, 'xpaset -p ds9 rgb blue'
  spawn, 'xpaset -p ds9 file tmpBlu.fits'
  spawn, 'xpaset -p ds9 scale limits '+string(min)+' '+string(max)
  spawn, 'xpaset -p ds9 scale '+scale
  spawn, 'xpaset -p ds9 export tiff '+outname

  im = read_tiff(outname)
  s = size(im, /dim)

  plot, findgen(s[1]), findgen(s[2]), /iso, /nodat, $
        title = imstack, xtitle = 'pix', ytitle = 'pix'
  cgimage, im, /over, stretch = 2
  
end
