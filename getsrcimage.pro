;; Point to Takahiro's convolced cluster images given a real
;; cluster input name
function getClusterImData, clusterName

  tPrefix = '/data/code/mtaka/HFF/Public/IMAGING/'

  clusterName = strupcase(clusterName)
  
  ; [01,02,03,04,05,06,07,08,09,10,99] =
  ; [A2744,M0416,M0717,M1149,A0370,A1063,M0744,M1423,M2129,R1347,XDF]
  case clusterName of
     'A2744'   : imdir = tPrefix+'01/' 
     'ABEL2744': imdir = tPrefix+'01/' 
     'M0416'   : imdir = tPrefix+'02/'
     'MACS0416': imdir = tPrefix+'02/'
     'M0717'   : imdir = tPrefix+'03/'
     'MACS0717': imdir = tPrefix+'03/'
     'M1149'   : imdir = tPrefix+'04/'
     'MACS1149': imdir = tPrefix+'04/'
     'A0370'   : imdir = tPrefix+'05/'
     'ABEL0370': imdir = tPrefix+'05/'
     'A1063'   : imdir = tPrefix+'06/'
     'RXJC2248': imdir = tPrefix+'06/'
     'M0744'   : imdir = tPrefix+'07/'
     'MACS0744': imdir = tPrefix+'07/'
     'M1423'   : imdir = tPrefix+'08/'
     'MACS1423': imdir = tPrefix+'08/'
     'M2129'   : imdir = tPrefix+'09/'
     'MACS2129': imdir = tPrefix+'09/'
     'R1347'   : imdir = tPrefix+'10/'
     'RXJC1347': imdir = tPrefix+'10/'
  endcase

  spawn, 'ls '+imdir+'cls_f???w_v*_???.fits > tmp.data'
  readcol, 'tmp.data', files, f = 'A', /silent, /quick
  nfiles = n_elements(files)

  filters = strarr(nfiles)
  for ii = 0, nfiles - 1 do $
     filters[ii] = strmid(files[ii], strpos(files[ii], 'cls_')+5, 3)

  filters = filters[sort(filters)]
  uf      = UNIQ(filters)
  nuf     = n_elements(uf)
  filters = filters[uf]

  ims = strarr(nuf)
  rms = strarr(nuf)
  for ii = 0, nuf - 1 do begin
     rms[ii] = files[uf[ii]]
     ims[ii] = repstr(files[uf[ii]], 'rms', 'drz')
  endfor
     
  savedata = {IMAGES  : ims,$
              RMS_MAPS: rms, $
              FILTERS : filters}

  RETURN, savedata
end

;; Correlate the SOBs...
;;
function getAng, p, $
                 imToMap = imToMap, $                 
                 imToMatch  = imToMatch, $
                 Err = Err

  a     = p[0]
  dx    = p[1]
  dy    = p[2]
  theta = p[3]
  
  tim = shift(rot(imToMap * a, theta, /cubic, /interp), dx, dy)

;  stop

  resids = (imToMatch - tim)/Err
  resids = reform(reform(resids, n_elements(imToMatch), 1))
  
  RETURN, resids
end

;; Get cutouts of an arbitrary size for a source given an RA, DEC,
;; and image name.
;; OK -- I think what you have to do is READ IN ALL THE STAMPS, make a
;;       stack, and then use *that* to find image offsets, BECAUSE THE
;;       GLASS DIRECT IMAGES COME FROM SOME FUCKING IR STACK.
pro getsrcimage, clusterName, ra, dec, $
                 MASTER = master, $
                 SIZE   = size, $
                 ROTATE = rotate, $
                 PA_OUT = pa_out, $
                 DS9    = ds9, $
                 OUTPUT = output

  if keyword_set(DS9) then begin
     spawn, 'ds9 &'
     wait, 1
  endif
  if NOT keyword_set(ROTATE) AND keyword_set(PA_OUT) $
  then rotate = 1 else $
     if NOT keyword_set(ROTATE) then rotate = 0
  if NOT keyword_set(OUTPUT) then $
     output = repstr(master, '.fits', '_allbandImages.fits')
  
  files  = getClusterImData(clusterName)
  images = files.IMAGES
  rmses  = files.RMS_MAPS
  filts  = files.FILTERS

  ;; Figure out what the grism saw
  ;; Ok -- we're gonna have to cross-correlate because shit's bananas
  mdata   = mrdfits(master, 1)  
  mstamp   = mdata.DIRECT_IM
  mstamp  *= mdata.EXPTIME
  mstmperr = mrdfits(mdata.FILENAME, 'dwht') > sqrt(abs(mstamp))
  mstmperr *=  mdata.EXPTIME
  foo      = mrdfits(mdata.FILENAME, 'sci', head)
  pixscale = sxpar(head, 'CD2_2')
  s = size(mstamp, /dim)
  width  = s[0] * pixscale ;; This is what we gotta cut out of the multiband data
  height = s[1] * pixscale
  if n_elements(ra)  eq 0 then ra  = mdata.RA
  if n_elements(dec) eq 0 then dec = mdata.DEC
  
  get_lun, lun
  openw, lun, 'foo.reg'
  printf, lun, 'fk5'
  printf, lun, f = '(%"box(%f,%f,%f,%f,0) # color = red")', $
          ra, dec, width/3600., height/3600.
  close, lun
  
  nfilters = n_elements(filts)
;  !p.multi = [0,nfilters,0]
  window, 0, xsize = 1500, ysize = 200
  multiplot, [nfilters,1]
  for ii = 0, nfilters - 1 do begin

     im     = readfits(images[ii], head, /silent)
     pixscl = abs(sxpar(head, 'CD1_1')) * 3600
     
     spawn, 'xpaset -p ds9 file '+images[ii]
     wait, 1
     spawn, 'xpaset -p ds9 regions deleteall'
     spawn, 'xpaset -p ds9 regions file foo.reg'
     spawn, 'xpaset -p ds9 pan to '+string(ra)+' '+string(dec)+' wcs fk5'
     spawn, 'xpaset -p ds9 scale limits 0 1'
     spawn, 'xpaset -p ds9 scale log'
     spawn, 'xpaget ds9 regions -coord image -format pros > foo2.dat'
     readcol, 'foo2.dat', x, y, f = 'X,X,D,D'
     x = x[0] & y = y[0]
     
     if ii eq 0 then begin
        sims = dblarr(s[0], s[1], nfilters) ;; It's gonna get sampled at the master pixscale
        srms = dblarr(s[0], s[1], nfilters)
        ets  = fltarr(nfilters)
     endif     

     if (1 - pixscale/pixscl) ge 0.1 then begin
        sideLenX = ceil(width / pixscl)                                     ;; asec / asec / pix
        sideLenY = ceil(height / pixscl)                                    ;; asec / asec / pix
        tim      = im[x-sidelenX/2:x+sidelenX/2-1, y-sidelenY/2:y+sidelenY/2-1] ;; Proper area stamp
        tim2     = congrid(tim, s[0], s[1], /center, /interp, cubic = -0.5)     ;; Proper pixscale stamp
        im       = temporary(readfits(rmses[ii]))
        etim     = im[x-sidelenX/2:x+sidelenX/2-1, y-sidelenY/2:y+sidelenY/2-1]
        etim2    = congrid(etim, s[0], s[1], /center, /interp, cubic = -0.5)
     endif else begin
        sideLenX = s[0]
        sideLenY = s[1]
        tim2     = im[x-sidelenX/2:x+sidelenX/2-1, y-sidelenY/2:y+sidelenY/2-1] ;; Proper area stamp
        im       = temporary(readfits(rmses[ii]))
        etim2    = im[x-sidelenX/2:x+sidelenX/2-1, y-sidelenY/2:y+sidelenY/2-1]
     endelse     

     sims[*,*,ii] = tim2
     srms[*,*,ii] = etim2
     ets[ii] = sxpar(head, 'EXPTIME')
     
     plot, findgen(s[0]), findgen(s[1]), /nodat, /iso
     cgimage, tim2, /over, stretch = 2
     multiplot
     
  endfor
  
  irbs = where(filts le 160, compl = obands)
  irstack = mean(sims[*,*,irbs], dim = 3)   ;; This *roughly* corresponds to the GLASS alignment image
  opstack = mean(sims[*,*,obands], dim = 3)

  tim2 = irstack; + opstack)                ;reform(sims[*,*, where(filts eq 140)])
  
  ;; Find the damn object
  foo = max(mstamp, hit1)
  foo = max(tim2, hit2)
  ctr1 = array_indices(mstamp, hit1)
  ctr2 = array_indices(tim2, hit2)
  dx = ctr1[0] - ctr2[0]
  dy = ctr1[1] - ctr2[1]

  scale = total(mstamp[20:40,20:40]) / total(tim2[20:40,20:40])
  tim2 *= scale
  resids = []
  for ii = 0, 359 do $
     resids = [resids, $
               total(((mstamp - rot(shift(tim2, dx,dy), ii, /cubic, /interp))/mstmperr)^2)]

  foo = min(resids, rang)
  thetas = indgen(360)
  rang = thetas[rang]

  tim3 = rot(tim2, rang, /cubic, /interp)
; 
;  p0 = double([1, 0, 0, rang])
;  parinfo  = {LIMITS: [0,0], LIMITED:1}
;  parinfo = replicate(parinfo, n_elements(p0))
;  parinfo[0].LIMITS  = [0,100]
;  parinfo[1].LIMITS  = [-5,5]
;  parinfo[2].LIMITS  = [-5,5]
;  parinfo[3].LIMITS  = [0,360]
;  functargs = {imToMap: tim2, imToMatch: mstamp, Err: mstmperr}    
;  scalepars = mpfit('getAng', p0, functargs = functargs, $
;                    /fastnorm, status = stat, xtol = 1d-20, ftol = 1d-20, gtol = 1d-20)
  
;  !p.multi = [0,4,0]
;  plot, findgen(s[0]) * pixscale, findgen(s[1]) * pixscale, /iso, $
;        /xsty, /ysty
;  cgimage, mstamp, /over, stretch = 2
;  plot, findgen(s[0]) * pixscl, findgen(s[1]) * pixscl, /iso, $
;        /xsty, /ysty
;  cgimage, tim, /over, stretch = 2
;  plot, findgen(s[0]) * pixscl, findgen(s[1]) * pixscl, /iso, $
;        /xsty, /ysty
;  cgimage, tim2, /over, stretch = 2
;  plot, findgen(s[0]) * pixscl, findgen(s[1]) * pixscl, /iso, $
;        /xsty, /ysty
;  cgimage, shift(rot(tim2, scalepars[3], /cubic, /interp), scalepars[1], scalepars[2]), /over, stretch = 2
;  !p.multi = 0
  
  ;; dump 'em
  newIms = sims
  newErs = srms     
  window, 1, xsize = 1500, ysize = 400
;  !p.multi = [0, nfilters+2, 0]
  multiplot, [nfilters+2, 1]
  plot, findgen(s[0]), findgen(s[1]), /iso, /nodat
  cgimage, mstamp, /over, stretch = 2
  multiplot
  for ii = 0, nfilters - 1 do begin
     newIms[*,*,ii] = rot(sims[*,*,ii], rang, /cubic, /interp)
     newErs[*,*,ii] = rot(srms[*,*,ii], rang, /cubic, /interp)
     plot, findgen(s[0]), findgen(s[1]), /iso, /nodat
     cgimage, newIms[*,*,ii], /over, stretch = 2
     multiplot
  endfor
  plot, findgen(s[0]), findgen(s[1]), /iso, /nodat
  cgimage, tim3, /over, stretch = 2
  cleanplot
  !p.multi = 0

  savedata = {IMAGE: dblarr(s[0],s[1]), $
              ERR:   dblarr(s[0],s[1]), $
              FILTER: 0, $
              EXPTIME: 0.}
  savedata = replicate(savedata, nfilters)
  for ii = 0, nfilters - 1 do begin
     savedata[ii].IMAGE = reform(newIms[*,*,ii])
     savedata[ii].ERR = reform(newErs[*,*,ii])
     savedata[ii].FILTER = filts[ii]
     savedata[ii].EXPTIME = ets[ii]
  endfor
  mwrfits, savedata, output, /create
  
end
;getSRCimage, 'ABEL2744', master = 'ABEL2744/00793_1_B.fits'
;getSRCimage, 'MACS1149', master = 'MACS1149/00753_1_B.fits'
;getSRCimage, 'MACS0744', master = 'MACS0744/00660_2_B.fits'
;getSRCimage, 'MACS1149', master = 'MACS1149/00900_1_B.fits'

;;
;;
