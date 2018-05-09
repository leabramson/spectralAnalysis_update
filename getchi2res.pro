pro getchi2res, lgnlist, explist

  readcol, lgnlist, lgnchains, f ='A'
  readcol, explist, expchains, f ='A'

  test = mrdfits(lgnchains[0], 1)

  regions = test.REGION
  chisPs = fltarr(n_elements(lgnchains), n_elements(regions), 2)
  chisSs = fltarr(n_elements(lgnchains), n_elements(regions), 2)

  for ii = 0, n_elements(lgnchains) - 1 do begin
     ldata = mrdfits(lgnchains[ii], 1)
     edata = mrdfits(expchains[ii], 1)

     chisPs[ii,*,0] = median(ldata.CHI_PHOT, dim = 1)
     chisSs[ii,*,0] = median(ldata.CHI_SPEC, dim = 1)

     chisPs[ii,*,1] = median(edata.CHI_PHOT, dim = 1)
     chisSs[ii,*,1] = median(edata.CHI_SPEC, dim = 1)
  endfor

  plot, chisSs[0,*,0] - chisSs[0,*,1], xtickname = regions, psym = 1, $
        /nodat, yran = [-60,70], xran = [-0.1,4.1], /xsty
  cols = ['ffa500'x, '0000ff'x, 'ff0000'x, '0055ff'x]
  oplot, !X.CRANGE, [1,1]
  for ii = 0, n_elements(lgnchains) - 1 do $
     oplot, chisSs[ii,*,0] - chisSs[ii,*,1], psym = 1, col = cols[ii]

  print, f = '(%"\n %s %s %s %s")', $
         regions
  for ii = 0, n_elements(lgnchains) - 1 do begin
     print, ''
     print, transpose(chisPs[ii,*,0])
     print, transpose(chisPs[ii,*,1])
     print, transpose(chisSs[ii,*,0])
     print, transpose(chisSs[ii,*,1])
     print, ''
     print, transpose(chisPs[ii,*,0] - chisPs[ii,*,1]) / chisqr_cvf(0.32, 17) 
     print, transpose(chisSs[ii,*,0] - chisSs[ii,*,1]) / chisqr_cvf(0.32, 300) 
     print, ''                                         
  endfor                                               

  stop
  
end
;getchi2res, 'study5_lgn_Chains.list', 'study5_exp_Chains.list'

