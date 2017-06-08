;; Plot the observed integrated SFMS and then that fixed to a common
;; reference redshift
pro plotsfms

  readcol, 'allLgnBestfits.list', files, f = 'A'
  nfiles = n_elements(files)
  
  masses  = fltarr(nfiles, 3)
  sfrs    = fltarr(nfiles, 3)
  sfrs_ha = fltarr(nfiles, 3)
  zs      = fltarr(nfiles, 3)

  for ii = 0, nfiles - 1 do begin
     data = mrdfits(files[ii], 1)
     opt = where(data.REGION eq 'OPTFL')
     data = data[opt]

     masses[ii,*]  = data.LMASS
     sfrs[ii,*]    = data.SFR
     zs[ii,*]      = data.Z
     sfrs_ha[ii,*] = data.HA_FLUX * 1d-18 * dluminosity(zs[ii,0], /cm)^2 * 4 *!pi / 1.26d41
  endfor
  
;  z_common = median(zs[*,0]);
;
;  readcol, 'allRecons.list', files, f = 'A'
;  nfiles = n_elements(files)
;  cmasses  = fltarr(nfiles, 3)
;  csfrs    = fltarr(nfiles, 3)
;  for ii = 0, nfiles - 1 do begin
;     data = mrdfits(files[ii], 1)
;     hit = value_locate(data.OPTFL_REDSHIFT, z_common)
;          
;     cmasses[ii,*] = alog10(data.OPTFL_MGH[hit,*])
;     csfrs[ii,*]    = data.OPTFL_SFH[hit,*]
;  endfor
    
  plotsym, 0, /fill
;  !p.multi = [0,2,0]
  plot, masses[*,0], alog10(sfrs[*,0]), /nodat, $
        xran = [10,11.5], yran = [-1,2], $
        xtitle = 'log !18M!X!D*!N [M!D'+sunsymbol()+'!N]', $
        ytitle = 'log !18SFR!X [M!D'+sunsymbol()+'!N yr!E-1!N]', /iso
  oploterror, masses[*,0], alog10(sfrs[*,0]), $
              0.5 * (masses[*,2] - masses[*,1]), $
              0.5/alog(10) * (sfrs[*,2] - sfrs[*,1])/sfrs[*,0], $
              psym = 8, col = '777777'x, errcol = '777777'x

  stop

  !p.multi =[0,3,nfiles/3+1]
  regions = ['INNFL', 'INTUP', 'INTDN', 'OUTUP', 'OUTDN', 'OPTFL']
  cols = [255, '00ff00'x, '007700'x, 'ffa500'x, 'ff0000'x, '777777'x]
  for ii = 0, nfiles - 1 do begin
     data = mrdfits(files[ii], 1)     
     qui = where(data.REGION eq 'OPTFL')
     tsfrl = data.HA_FLUX[0] * 1d-18 * dluminosity(zs[ii,0], /cm)^2 * 4 *!pi / 1.26d41
     tssfrl = alog10(tsfrl) - data.LMASS[0]
     plot, [data[qui].LMASS[0]], alog10([data[qui].SFR[0]]) - data[qui].LMASS, $
;           xran = minmax(data.LMASS), yran = alog10(minmax(data.SFR[0]/10.^data.LMASS[0])), $
;           xran = minmax(data.LMASS), yran = minmax(tssfrl), $
           xran = [9,11.5], yran = [-12,-8.5], $
           xtitle = 'log !18M!X!D*!N [M!D'+sunsymbol()+'!N]', $
           ytitle = 'log !18sSFR!X [yr!E-1!N]', $
           /nodat
     for jj = 0, n_elements(regions) - 1 do begin
        qui = where(data.REGION eq REGIONS[JJ])
        tdata = data[qui]
        mass = tdata.LMASS[0]
        sfr = tdata.SFR[0]
        ssfr    = alog10(sfr) - mass
        mass_err = 0.5 * (tdata.LMASS[2] - tdata.LMASS[1])
        sfr_err  = 0.5/alog(10) * (tdata.SFR[2] - tdata.SFR[1])/tdata.SFR[0]
        ssfr_err = sqrt(mass_err^2 + sfr_err^2)
        hasfr = tdata.HA_FLUX * 1d-18 * dluminosity(zs[ii,0], /cm)^2 * 4 *!pi / 1.26d41
        hasfr_err = 0.5/alog(10) * (hasfr[2] - hasfr[1]) / hasfr[0]
        hassfr_err = sqrt(hasfr_err^2 + mass_err^2)
        hassfr = alog10(hasfr) - mass        
        oploterror, [mass], [hassfr], mass_err, hassfr_err, $
                    psym = 8, col = cols[jj], errcol = cols[jj]
        oplot, [mass], [ssfr], psym = 8, symsize = 1.3
        oploterror, [mass], [ssfr], mass_err, ssfr_err, $
                    psym = 8, col = cols[jj], errcol = cols[jj]
        legend, /bottom, /right, box = 0, $
                repstr(files[ii], '_lgn_bestfit.fits', ''), $
                charsize = 1
     endfor
  endfor

;  sfrs    = fltarr(n_elements(regions), nfiles, 3)
;  sfrs_Ha = fltarr(n_elements(regions), nfiles, 3)
;  for ii = 0, nfiles - 1 do begin
;
;     data = mrdfits(files[ii], 1)     
;     for jj = 0, n_elements(regions) - 1 do begin
;        tdata = data[where(data.REGION eq regions[jj])]
;        sfrs[jj,ii,*] = tdata.SFR
;        sfrs_Ha[jj,ii,*] = tdata.HA_FLUX * 1d-18 * dluminosity(zs[ii,0], /cm)^2 * 4 *!pi / 1.26d41
;     endfor
;     
;  endfor
;
;  !p.multi = 0
;  
;  plot, alog10(sfrs), alog10(sfrs_Ha), $
;        /iso, /nodat, xran = [-4,4], yran = [-4,4], /xsty, /ysty
;  one_one
;  for jj = 0, n_elements(regions) - 1 do $
;     oploterror, alog10(sfrs[jj,*,0]), alog10(sfrs_ha[jj,*,0]), $
;                 0.5/alog(10) * (sfrs[jj,*,2] - sfrs[jj,*,1])/sfrs[jj,*,0], $
;                 0.5/alog(10) * (sfrs_ha[jj,*,2] - sfrs_ha[jj,*,1])/sfrs_Ha[jj,*,0], $
;                 col = cols[jj], errcol = cols[jj], psym = 8
  
;  oploterror, masses[*,0], alog10(sfrs_ha[*,0]), $
;              0.5 * (masses[*,2] - masses[*,1]), $
;              0.5/alog(10) * (sfrs_ha[*,2] - sfrs_ha[*,1])/sfrs_ha[*,0], $
;              psym = 8, col = '0000ff'x, errcol = '0000ff'x

  
  
  stop
  
end
