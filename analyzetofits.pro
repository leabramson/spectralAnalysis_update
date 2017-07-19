;; Take all of the .analyze files for a single source and convert them
;; into a more convenient, monolithic FITS table

pro analyzetofits, inprefix, output, $
                   FITSFILE = fitsfile, $
                   NOTMASKED = notmasked, $
                   SFH = sfh

  if NOT keyword_set(NOTMASKED) then $
     notmasked = 0 $
  else $
     notmasked = 1
  
  if notmasked then $
     spawn, 'ls '+inprefix+'/?????.analyze > tmp.txt' $
  else $
     spawn, 'ls '+inprefix+'/?????.masked.analyze > tmp.txt'

  if NOT keyword_set(SFH) then sfh = 'DelayedExponential'
  
  readcol, 'tmp.txt', files, f = 'A', /silent
  nregions = n_elements(files)
  
  regions = strarr(nregions)
  for ii = 0, nregions - 1 do begin
     treg = repstr(files[ii], inprefix+'/', '')
     if notmasked then $
        regions[ii] = repstr(treg, '.analyze', '') $
     else $
        regions[ii] = repstr(treg, '.masked.analyze', '')
  endfor
  
  if sfh eq 'DelayedExponential' then begin

     ; SII6718_ew   z      SII6733_ew
     ; OIII_ew      H_ew   A_V
     ; age          ltau   logmstar   logssfr 
     ; sfr   Ha_flux   Hb_flux
     ; OIII_flux   SII6718_flux   SII6733_flux 
     
     params = ['z','logmstar','sfr', $
               'H_ew','Ha_flux','Hb_flux',$
               'OIII_ew','OIII_flux', $
               'age', 'ltau', 'A_V', $
               'evidence', 'chi2phot', 'chi2spec']
     
     savedata = {REGION: string(0), $
                 Z: fltarr(3), $
                 LMASS: fltarr(3), $
                 SFR: fltarr(3), $
                 HA_EW: fltarr(3), $
                 HA_FLUX:fltarr(3), $
                 HB_FLUX: fltarr(3) , $
                 OIII_EW: fltarr(3), $
                 OIII_FLUX: fltarr(3), $
                 AGE: fltarr(3), $
                 TAU: fltarr(3), $
                 A_V: fltarr(3), $
                 EVIDENCE: 0., $
                 CHI2PHOT: [0.,0.], $
                 CHI2SPEC: [0.,0.], $
                 RNORM: 0., $
                 RE_OBS: 0., $
                 RE_PHYS: 0., $
                 AREA_OBS: 0., $
                 AREA_PHYS: 0., $
                 FITSFILE: fitsfile, $
                 TIFF_IMAGE: repstr(fitsfile, '_B.fits', '_rgb.tiff')}
     
  endif else if sfh eq 'LogNormal' then begin

     params = ['logmstar','sfr', $
               'freeHa_ew', 'freeHb_ew', $
               'freeHa_flux','freeHb_flux',$
               'OIII_ew','OIII_flux', $
               'T0', 'tau', 'A_V', $
               'evidence', 'chi2phot', 'chi2spec']
     
     savedata = {REGION: string(0), $
                 Z: fltarr(3), $
                 LMASS: fltarr(3), $
                 SFR: fltarr(3), $
                 HA_EW: fltarr(3), $
                 HB_EW: fltarr(3), $
                 HA_FLUX:fltarr(3), $
                 HB_FLUX: fltarr(3) , $
                 OIII_EW: fltarr(3), $
                 OIII_FLUX: fltarr(3), $
                 T0: fltarr(3), $
                 TAU: fltarr(3), $
                 A_V: fltarr(3), $
                 EVIDENCE: 0., $
                 CHI2PHOT: [0.,0.], $
                 CHI2SPEC: [0.,0.], $
                 RNORM: 0., $
                 RE_OBS: 0., $
                 RE_PHYS: 0., $
                 AREA_OBS: 0., $
                 AREA_PHYS: 0., $
                 FITSFILE: fitsfile, $
                 TIFF_IMAGE: repstr(fitsfile, '_B.fits', '_rgb.tiff')}
     
  endif

  savedata = replicate(savedata, nregions)
  for ii = 0, nregions - 1 do $
     savedata[ii].REGION = regions[ii]

  if sfh eq 'LogNormal' then begin
     inprefix2 = repstr(inprefix, 'LogNormal', 'Phot') ;; Lognormals don't have a redshift

     if notmasked then $
        tfile = inprefix2+'/OPTFL.analyze' $
     else $
        tfile = inprefix2+'/OPTFL.masked.analyze'

     spawn, 'cat '+tfile+' | grep z > tmp.txt'
     readcol, 'tmp.txt', f, flo, fhi, f = 'X,F,F,F', /silent
     savedata[*].Z         = [f[0], flo[0], fhi[0]]
  endif
  
  master = where(regions eq 'OPTFL') ;; base everything off the optimum extraction
  print, f = '(%"\n\n Reading data ... \n\n")'
  for jj = 0, n_elements(params) - 1 do begin
     spawn, 'cat '+files[master]+' | grep '+params[jj]+' > tmp.txt'
     readcol, 'tmp.txt', f, flo, fhi, f = 'X,F,F,F', /silent
     case params[jj] of
        'z'          : savedata[master].Z         = [f[0], flo[0], fhi[0]]
        'logmstar'   : savedata[master].LMASS     = [f[0], flo[0], fhi[0]]
        'sfr'        : savedata[master].SFR       = [f[1], flo[1], fhi[1]] ;; gets logssfr too...
        'H_ew'       : savedata[master].HA_EW     = [f[0], flo[0], fhi[0]]
        'freeHa_ew'  : savedata[master].HA_EW     = [f[0], flo[0], fhi[0]]
        'freeHb_ew'  : savedata[master].HB_EW     = [f[0], flo[0], fhi[0]]
        'Ha_flux'    : savedata[master].HA_FLUX   = [f[0], flo[0], fhi[0]]
        'Hb_flux'    : savedata[master].HB_FLUX   = [f[0], flo[0], fhi[0]]
        'freeHa_flux': savedata[master].HA_FLUX   = [f[0], flo[0], fhi[0]]
        'freeHb_flux': savedata[master].HB_FLUX   = [f[0], flo[0], fhi[0]]
        'OIII_ew'    : savedata[master].OIII_EW   = [f[0], flo[0], fhi[0]]
        'OIII_flux'  : savedata[master].OIII_FLUX = [f[0], flo[0], fhi[0]]
        'age'        : savedata[master].AGE       = [f[0], flo[0], fhi[0]]
        'ltau'       : savedata[master].TAU       = 10.^[f[0], flo[0], fhi[0]]
        'A_V'        : savedata[master].A_V       = [f[0], flo[0], fhi[0]]
        'T0'         : savedata[master].T0        = [f[0], flo[0], fhi[0]]
        'tau'        : savedata[master].TAU       = [f[0], flo[0], fhi[0]]
        'evidence'   : begin
           readcol, 'tmp.txt', f, f = 'X,F', /silent
           savedata[master].EVIDENCE       = f[0]
        end
        'chi2phot' : begin
           readcol, 'tmp.txt', chi2, n, f = 'X,X,F,F', delimiter = '/=', /silent
           savedata[master].CHI2PHOT       = [chi2,n]
        end
        'chi2spec' : begin
           readcol, 'tmp.txt', chi2, n, f = 'X,X,F,F', delimiter = '/=', /silent
           savedata[master].CHI2SPEC       = [chi2,n]
        end
     endcase
  endfor

  radii = getradii(fitsfile)
  areas = getarea(fitsfile, savedata[master].Z[0])
  
  for ii = 0, nregions - 1 do begin
    
     if regions[ii] ne 'OPTFL' then begin

        hit = where(radii.REGION eq regions[ii])
        savedata[ii].RNORM   = radii[hit].RNORM
        savedata[ii].RE_OBS  = radii[hit].RE_OBS
        savedata[ii].RE_PHYS = radii[hit].RE_OBS / $
                               zang(1., savedata[where(savedata.REGION eq 'OPTFL')].Z[0], /silent)

        hit = where(areas.REGION eq regions[ii])
        savedata[ii].AREA_OBS  = areas[hit].AREA_OBS
        savedata[ii].AREA_PHYS = areas[hit].AREA_PHYS
        
        for jj = 0, n_elements(params) - 1 do begin
           spawn, 'cat '+files[ii]+' | grep '+params[jj]+' > tmp.txt'
           readcol, 'tmp.txt', f, flo, fhi, f = 'X,F,F,F', /silent
           case params[jj] of
              'z'        : begin
                 if n_elements(f) gt 0 then $
                    savedata[ii].Z         = [f[0], flo[0], fhi[0]] $
                 else $
                    savedata[ii].Z         = savedata[master].Z
              end
              'logmstar' : savedata[ii].LMASS     = [f[0], flo[0], fhi[0]]
              'sfr'      : savedata[ii].SFR       = [f[1], flo[1], fhi[1]] ;; gets logssfr too...
              'freeHa_ew': savedata[ii].HA_EW     = [f[0], flo[0], fhi[0]]
              'freeHb_ew': savedata[ii].HB_EW     = [f[0], flo[0], fhi[0]]
              'H_ew'     : savedata[ii].HA_EW     = [f[0], flo[0], fhi[0]]
              'Ha_flux'  : savedata[ii].HA_FLUX   = [f[0], flo[0], fhi[0]]
              'Hb_flux'  : savedata[ii].HB_FLUX   = [f[0], flo[0], fhi[0]]
              'freeHa_flux': savedata[ii].HA_FLUX = [f[0], flo[0], fhi[0]]
              'freeHb_flux': savedata[ii].HB_FLUX = [f[0], flo[0], fhi[0]]
              'OIII_ew'  : savedata[ii].OIII_EW   = [f[0], flo[0], fhi[0]]
              'OIII_flux': savedata[ii].OIII_FLUX = [f[0], flo[0], fhi[0]]
              'age'      : savedata[ii].AGE       = [f[0], flo[0], fhi[0]]
              'ltau'     : savedata[ii].TAU       = 10.^[f[0], flo[0], fhi[0]]
              'A_V'      : savedata[ii].A_V       = [f[0], flo[0], fhi[0]]
              'T0'       : savedata[ii].T0        = [f[0], flo[0], fhi[0]]
              'tau'      : savedata[ii].TAU       = [f[0], flo[0], fhi[0]]
              'evidence' : begin
                 readcol, 'tmp.txt', f, f = 'X,F', /silent
                 savedata[master].EVIDENCE       = f[0]
              end
              'chi2phot' : begin
                 readcol, 'tmp.txt', chi2, n, f = 'X,X,F,F', delimiter = '/=', /silent
                 savedata[master].CHI2PHOT       = [chi2,n]
              end
              'chi2spec' : begin
                 readcol, 'tmp.txt', chi2, n, f = 'X,X,F,F', delimiter = '/=', /silent
                 savedata[master].CHI2SPEC       = [chi2,n]
              end
           endcase
        endfor
     endif
  endfor
  
  mwrfits, savedata, output, /create
  
end
