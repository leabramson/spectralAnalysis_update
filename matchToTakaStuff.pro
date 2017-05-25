pro matchToTakaStuff, cluster, fitsList, output

  addDir = '~/LocalData/Abramson/GLASS_official_products/'
  
  cluster = strupcase(cluster)
  addDir = addDir+cluster+'/Catalogs/'

  radecToID = addDir+'TAKA_IDtoRADEC.cat'
  fastStuff = addDir+'TAKA_FAST.output'
  morphoStuff = addDir+'TAKA_GALFIT.output'
;  case cluster of
;     'ABEL0370': 
;     'ABEL2744':
;     'MACS0416':
;     'MACS0717':
;     'MACS0744':
;     'MACS1149':
;     'MACS1423':
;     'MACS2129':
;     'RXJC1347':
;     'RXJC2248':
;  endcase
  
  readcol, fitsList, files, f = 'A'
  nfiles = n_elements(files)

  readcol, raDecToID, $
           id1, ra, dec, $
           f = 'L,D,D'
  readcol, fastStuff, $
           id2, z, ltau, metal, lage, Av, lmass, lsfr, lssfr, la2t, chi2, $
           f = 'L,F,F,F,F,F,F,F,F,F,F'
  readcol, morphoStuff, $
           id3, x, y, gmag, ae, n, q, pa, $
           f = 'L,F,F,F,F,F,F,F'

  close, 1
  openw, 1, output, width = 1024
  printf, 1, '#ID ID_TAKA RA DEC MAG GFMAG_TAKA Z_GLASS Z_TAKA MASS_TAKA SFR_TAKA N Re Q'
  for ii = 0, nfiles - 1 do begin
     data = mrdfits(files[ii], 1)

     mag     = data.MAG
     z_GLASS = data.Z
     tra     = data.RA
     tdec    = data.DEC
     tid     = strmid(files(ii), 0, 5)

     spherematch, ra, dec, $
                  tra, tdec, $
                  1.d/3600, $
                  hit, sinds, len
     if n_elements(hit) eq 1 then begin
        hit  = hit[0]
        hit2 = (where(id2 eq id1[hit], n2))[0]
        hit3 = (where(id3 eq id1[hit], n3))[0]

        if n2 eq 1 AND n3 eq 1 then $
           printf, 1, $
                   tid, ' ', tra, tdec, mag, gmag[hit3], $
                   z_GLASS, z[hit], lmass[hit2], lsfr[hit2], $
                   n[hit3], ae[hit3], q[hit3] $
        else if n2 eq 1 AND n3 eq 0 then $
           printf, 1, $
                   tid, ' ', tra, tdec, mag, -99, $
                   z_GLASS, z[hit], lmass[hit2], lsfr[hit2], $
                   -99, -99, -99
     endif else stop

  endfor
  close, 1

end

pro run1, cluster

  check = file_info('matchToTakaStuff.pro')
  if NOT check.EXISTS then $
     spawn, 'ln -s ../matchToTakaStuff.pro .'
  spawn, "cat pyspecWorthy_sourceInfo.list | grep '_B.fits' > tmp.list"
  matchToTakaStuff, cluster, 'tmp.list', 'fitworthyWithTakaData.list'
  spawn, "cat fitworthyWithTakaData.list"
  
end
