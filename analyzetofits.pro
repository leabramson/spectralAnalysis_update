;; Take all of the .analyze files for a single source and convert them
;; into a more convenient, monolitich FITS tabls

pro analyzetofits, inprefix, output, $
                   NOTMASKED = notmasked

  if NOT keyword_set(NOTMASKED) then notmasked = 0 else notmasked = 1

  if notmasked then $
     spawn, 'ls '+inprefix+'/?????.analyze > tmp.txt' $
  else $
     spawn, 'ls '+inprefix+'/?????.masked.analyze > tmp.txt'
  
  readcol, 'tmp.txt', files, f = 'A'
  nregions = n_elements(files)
  
  regions = strarr(nregions)
  for ii = 0, nregions - 1 do begin
     treg = repstr(files[ii], inprefix+'/', '')
     if notmasked then $
        regions[ii] = repstr(treg, '.analyze', '') $
     else $
        regions[ii] = repstr(treg, '.masked.analyze', '')
  endfor

;  SII6718_ew   z   SII6733_ew   OIII_ew   H_ew   A_V   age   ltau   logmstar   logssfr 
;  sfr   Ha_flux   Hb_flux   OIII_flux   SII6718_flux   SII6733_flux 

  params = ['z','logmstar','sfr','H_ew','Ha_flux','Hb_flux','OIII_ew','OIII_flux', $
           'age', 'ltau']

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
              TAU: fltarr(3)}
  savedata = replicate(savedata, nregions)
  for ii = 0, nregions - 1 do $
     savedata[ii].REGION = regions[ii]

  master = where(regions eq 'OPTFL') ;; base everything off the optimum extraction
  for jj = 0, n_elements(params) - 1 do begin
     spawn, 'cat '+files[master]+' | grep '+params[jj]+' > tmp.txt'
     readcol, 'tmp.txt', f, flo, fhi, f = 'X,F,F,F'
     case params[jj] of
        'z'        : savedata[master].Z         = [f[0], flo[0], fhi[0]]
        'logmstar' : savedata[master].LMASS     = [f[0], flo[0], fhi[0]]
        'sfr'      : savedata[master].SFR       = [f[1], flo[1], fhi[1]] ;; gets logssfr too...
        'H_ew'     : savedata[master].HA_EW     = [f[0], flo[0], fhi[0]]
        'Ha_flux'  : savedata[master].HA_FLUX   = [f[0], flo[0], fhi[0]]
        'Hb_flux'  : savedata[master].HB_FLUX   = [f[0], flo[0], fhi[0]]
        'OIII_ew'  : savedata[master].OIII_EW   = [f[0], flo[0], fhi[0]]
        'OIII_flux': savedata[master].OIII_FLUX = [f[0], flo[0], fhi[0]]
        'age'      : savedata[master].AGE       = [f[0], flo[0], fhi[0]]
        'ltau'     : savedata[master].TAU       = 10.^[f[0], flo[0], fhi[0]]
     endcase
  endfor

  for ii = 0, nregions - 1 do begin
     if regions[ii] ne 'OPTFL' then begin
        for jj = 0, n_elements(params) - 1 do begin
           spawn, 'cat '+files[ii]+' | grep '+params[jj]+' > tmp.txt'
           readcol, 'tmp.txt', f, flo, fhi, f = 'X,F,F,F'
           case params[jj] of
              'z'        : begin
                 if n_elements(f) gt 0 then $
                    savedata[ii].Z         = [f[0], flo[0], fhi[0]] $
                 else $
                    savedata[ii].Z         = savedata[master].Z
              end
              'logmstar' : savedata[ii].LMASS     = [f[0], flo[0], fhi[0]]
              'sfr'      : savedata[ii].SFR       = [f[1], flo[1], fhi[1]] ;; gets logssfr too...
              'H_ew'     : savedata[ii].HA_EW     = [f[0], flo[0], fhi[0]]
              'Ha_flux'  : savedata[ii].HA_FLUX   = [f[0], flo[0], fhi[0]]
              'Hb_flux'  : savedata[ii].HB_FLUX   = [f[0], flo[0], fhi[0]]
              'OIII_ew'  : savedata[ii].OIII_EW   = [f[0], flo[0], fhi[0]]
              'OIII_flux': savedata[ii].OIII_FLUX = [f[0], flo[0], fhi[0]]
              'age'      : savedata[ii].AGE       = [f[0], flo[0], fhi[0]]
              'ltau'     : savedata[ii].TAU       = 10.^[f[0], flo[0], fhi[0]]
           endcase
        endfor
     endif
  endfor

  mwrfits, savedata, output, /create
  
end
