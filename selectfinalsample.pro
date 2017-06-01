pro selectfinalsample, output

  spawn, $
     'ls ????????/?????_?_pyspecfitPhotResults/OPTFL.masked.analyze > tmpRes.txt'
  readcol, 'tmpRes.txt', files, f = 'A'
  nfiles = n_elements(files)

  hit = []
  zs  = []
  plot, findgen(60), findgen(60), /iso, /nodat
  for ii = 0, nfiles - 1 do begin
     spawn, 'cat '+files[ii]+' | grep z > tmpZ.txt'
     readcol, 'tmpZ.txt', z, f = 'X,F'
     print, f = '(%"%s : %f")', $
            files[ii], z
     im = read_tiff(repstr(files[ii], '_pyspecfitPhotResults/OPTFL.masked.analyze', '_rgb.tiff'))
     cgimage, im, /over, stretch = 5
     galaxy = ''
     read, galaxy, prompt = 'Is this a galaxy? [Y]/n : '
     galaxy = strupcase(galaxy)
     if z gt 1 AND z lt 2 AND galaxy ne 'N' then begin
        hit = [hit, ii]
        zs = [zs, z]
     endif
  endfor

  print, f = '(%"\n\n FILES TO KEEP : \n\n")'
  print, files[hit]

  close, 1
  openw, 1, output, width = 64
  printf, 1, '# '+systime()
  printf, 1, '# FINAL FILE SELECTION; 1 < z < 2'
  for ii = 0, n_elements(hit) - 1 do $
     printf, 1, f = '(%"%s %f")', $
             repstr(files[hit[ii]], '_pyspecfitPhotResults/OPTFL.masked.analyze', ''), zs[ii]
  close, 1
  
end

;;selectfinalsample, '20170530_finalSampleForAnalysis.list'
;;
;;

pro cullworthy

  readcol, '20170530_finalSampleForAnalysis.list', $
           prefix, z, f= 'A,F'
  nfiles = n_elements(prefix)
  for ii = 0, nfiles - 1 do begin
     cluster = strmid(prefix[ii], 0,8)
     spawn, 'rm '+cluster+'/forLognormalFitting.list'
  endfor
  for ii = 0, nfiles - 1 do begin
     cluster = strmid(prefix[ii], 0,8)
     source  = strmid(prefix[ii], 9,7)
     check = file_info(cluster+'/forLognormalFitting.list')
     if NOT check.EXISTS then $
        spawn, 'cat '+cluster+'/pyspecWorthy_sourceInfo.list | grep '+source+' > '+cluster+'/forLognormalFitting.list' $
     else $
        spawn, 'cat '+cluster+'/pyspecWorthy_sourceInfo.list | grep '+source+' >> '+cluster+'/forLognormalFitting.list' 
  endfor
  
end
