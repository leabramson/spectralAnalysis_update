pro getHSTimaging, cluster, SURVEY = survey

  opf = ['f814w']
  irf = ['f125w', 'f160w']

  print, ''
  print, 'Getting '+survey+' imaging for bands : '+opf[0:*]+' | '+irf[0:*]
  print, 'for cluster : '+cluster
  print, ''
  
  if survey eq 'HFF' then begin

     for ii = 0, n_elements(irf) - 1 do begin
        print, ' >>> Getting '+irf[ii]+' data ... '
        spawn, 'wget https://archive.stsci.edu/pub/hlsp/frontier/'+cluster+'/images/hst/v1.0-epoch1/'+$
               'hlsp_frontier_hst_wfc3-60mas_'+cluster+'_'+irf[ii]+'_v1.0-epoch1_drz.fits'
     endfor

     print, '' 

     for ii = 0, n_elements(opf) - 1 do begin
        print, ' >>> Getting '+opf[ii]+' data ... '
        spawn, 'wget https://archive.stsci.edu/pub/hlsp/frontier/'+cluster+'/images/hst/v1.0-epoch1/'+$
               'hlsp_frontier_hst_acs-60mas_'+cluster+'_'+opf[ii]+'_v1.0-epoch1_drz.fits'
     endfor

     print, ''
     
  endif else if survey eq 'clash' then begin

     for ii = 0, n_elements(irf) - 1 do begin
        print, ' >>> Getting '+irf[ii]+' data ... '
        spawn, 'wget https://archive.stsci.edu/missions/hlsp/clash/'+cluster+'/data/hst/scale_65mas/'+$
               'hlsp_clash_hst_wfc3ir-65mas_'+cluster+'_'+irf[ii]+'_v1_drz.fits'
     endfor

     print, ''

     for ii = 0, n_elements(opf) - 1 do begin
        print, ' >>> Getting '+opf[ii]+' data ... '
        spawn, 'wget https://archive.stsci.edu/missions/hlsp/clash/'+cluster+'/data/hst/scale_65mas/'+$
               'hlsp_clash_hst_acs-65mas_'+cluster+'_'+opf[ii]+'_v1_drz.fits'
     endfor

     print, ''
     
  endif

  print, 'DATA OBTAINED'
  print, ''
  
end

pro cycle


;  getHSTimaging, 'macs1423', survey = 'clash'
  getHSTimaging, 'macs1149', survey = 'HFF'
  getHSTimaging, 'macs2129', survey = 'clash'
  getHSTimaging, 'rxj1347' , survey = 'clash'
  
end
