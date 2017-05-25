function makeCutOut, imarr, x, y, $
                     SIZE = size, $
                     ROTANG = rotang
  
  RETURN, output
end

;;
;;

pro prepForCutOut, image, filter, $
                   TARGDIR = targdir, $
                   DS9 = ds9

  if keyword_SET(DS9) then begin
     spawn, 'ds9 &'
     wait, 1
  endif
     
  ;; Load all the extracted spectra fits files
  ;; into a list
  spawn, 'ls '+targdir+'/*_B.fits > toget.list'
  readcol, 'toget.list', files, f = 'A'
  nfiles = n_elements(files)
  info = {ID: string(0), $
          RA: 0.d, $
          DEC: 0.d, $
          PA: 0., $
          Z: 0., $
          MAG: 0., $
          RE: 0., $
          PIXSCALE: 0.}
  info = replicate(info, nfiles)
  for ii = 0, nfiles - 1 do begin
     tdata = mrdfits(files[ii], 1)
     foo   = mrdfits(tdata.FILENAME, 'sci', thead, /silent)
     info[ii].ID  = strmid(files[ii], strpos(files[ii], '/')+1, strpos(files[ii], '_'))
     info[ii].RA  = tdata.RA
     info[ii].DEC = tdata.DEC
     info[ii].MAG = tdata.MAG
     info[ii].PA  = tdata.PA
     info[ii].Z   = tdata.Z
     info[ii].RE  = tdata.RE
     info[ii].PIXSCALE = sxpar(thead, 'CD2_2')
  endfor
     
  ;; Prep the image
  spawn, 'xpaset -p ds9 file '+image
  wait, 1
  im  = readfits(image, imhead)
  head = mrdfits(image, 1)
;  posAng = sxpar(head, 'PA_V3')
  if filter eq 'B' or filter eq 'V' or filter eq 'R' or filter eq 'I' then begin
     rs = 1
     posAng = head[0].PA_V3
  endif else begin
     rs = -1
     posAng = sxpar(imhead, 'PA_V3')
  endelse
  
  dPA    = info.PA - posAng
  dontRot = where(abs(dPA) lt 1., ndont)
  if ndont gt 0 then dPA[dontRot] = 0

  print, info.PA
  print, posAng
  print, dPA

;  stop
  
  for ii = 0, nfiles - 1 do begin
     close, 1
     openw, 1, 'foo.reg'
     printf, 1, 'fk5'
     printf, 1, 'circle(', info[ii].RA, ',', info[ii].DEC, ',', string(2./3600)+') # color = red'
     close, 1
     spawn, 'xpaset -p ds9 regions deleteall'
     spawn, 'xpaset -p ds9 regions file foo.reg'
     spawn, 'xpaget ds9 regions -coord image -format pros > foo2.dat'
     readcol, 'foo2.dat', x, y, f = 'X,X,D,D'

     x = x[0]
     y = y[0]

     if sxpar(imhead, 'CD1_1') lt 0 then revx = 1 else revx = 0
     if sxpar(imhead, 'CD2_2') lt 0 then revy = 1 else revy = 0
     dx = abs(3600. * sxpar(imhead, 'CD1_1'))
     
     foo = mrdfits(files[ii], 1, /silent)
     ts  = size(foo.DIRECT_IM, /dim)
     dp  = ts[0] / 2
     dp *= info[ii].PIXSCALE/dx

     cutOut = im[x-dp:x+dp,y-dp:y+dp]
;     if revx then cutOut = reverse(cutOut, 1)
;     if revy then cutOut = reverse(cutOut, 2)
     cutOut = congrid(cutOut, ts[0], ts[1])
     if rs eq 1 then $
        cutOut = rot(cutOut, dpa[ii], /interp, /cubic) $
     else begin
        if dpa[ii] gt 0 then $
           cutOut = rot(cutOut, 180 - dpa[ii], /interp, /cubic) $
        else $
           cutOut = rot(cutOut, 360 - dpa[ii], /interp, /cubic)
     endelse

     !p.multi = [0,3,0]
     plot, findgen(60), findgen(60), /nodat, /iso
     cgimage, foo.DIRECT_IM, /over, stretch = 7
     plot, findgen(60), findgen(60), /nodat, /iso
     cgimage, im[x-dp:x+dp,y-dp:y+dp], /over, stretch = 7
     plot, findgen(60), findgen(60), /nodat, /iso
     cgimage, cutout, /over, stretch = 7

     k = get_kbrd(1)
     
;     stop

     outname = targdir+'/'+$
               repstr(info[ii].ID, '.fits', '')+$
               '_'+string(info[ii].PA, f = '(I03)')+'_'+filter+'.fits'
     
     mwrfits, cutOut, outname, /create
  endfor
  

;  stop  
  
end

pro doCluster, cluster, targdir = targdir, survey

  if survey eq 'clash' then begin
     prepForCutout, '~/LocalData/Abramson/GLASSimages/hlsp_clash_hst_acs-65mas_'+cluster+'_f814w_v1_drz.fits', $
                    'I', targdir = targdir
     prepForCutout, '~/LocalData/Abramson/GLASSimages/hlsp_clash_hst_wfc3ir-65mas_'+cluster+'_f125w_v1_drz.fits', $
                    'J', targdir = targdir
     prepForCutout, '~/LocalData/Abramson/GLASSimages/hlsp_clash_hst_wfc3ir-65mas_'+cluster+'_f160w_v1_drz.fits', $
                    'H', targdir = targdir
  endif else if survey eq 'HFF' then begin
     prepForCutout, '~/LocalData/Abramson/GLASSimages/hlsp_frontier_hst_wfc3-30mas_'+cluster+'_f160w_v1.0-epoch1_drz.fits', $
                    'H', targdir = targdir
     prepForCutout, '~/LocalData/Abramson/GLASSimages/hlsp_frontier_hst_acs-30mas_'+cluster+'_f814w_v1.0-epoch1_drz.fits', $
                    'I', targdir = targdir
     prepForCutout, '~/LocalData/Abramson/GLASSimages/hlsp_frontier_hst_wfc3-30mas_'+cluster+'_f125w_v1.0-epoch1_drz.fits', $
                    'J', targdir = targdir

  endif
  
end
