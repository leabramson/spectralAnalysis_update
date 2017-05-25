;; Do one set of 5 extractions for a sigle object at a single PA
;; i.e., everything for "00001_1" in a given FIELD (e.g., "RXJC2248")

pro runpyspeconeset, field, setPrefix

  checklist = field+'/fileBreakdownFor1Dextraction.list'
  check = file_info(checklist)
  if NOT check.EXISTS then $
     organizeoneonefield, field

  readcol, checklist, $
           fname, ID, PA, grism, pz, z, comment = '#', $
           f = 'A,A,A,A,F,F'

  readyToGo = strarr(n_elements(ID))
  for ii = 0, n_elements(ID) - 1 do $
     readyToGo = id[ii]+'_'+pa[ii]

  hits = where(readyToGo eq setPrefix, nhits)

  if nhits eq 2 then begin

     ;; Find a "redshift"
     tpz = pz[hits[0]]
     tz  = z[hits[0]]

     if tpz lt 0 then begin
        zstart = tz
        zsran = 0.5
     endif else if tz lt 0 then begin
        zstart = tpz
        zsran = 0.5
     endif else if tz gt 0 AND tpz gt 0 then begin
        zstart = 0.5 * (tpz + tz)
        zsran = (sqrt(0.5) * abs(tpz - tz)) < 0.5
     endif

     pyspecfitonetrace, id[hits[0]], pa[hits[0]], $
                        'INNFL.masked', zstart, zsran, $
                        ncores = 10
     outname = id[hits[0]]+'_'+pa[hits[0]]+'_pyspecfitResults'
     check = file_info(outname)
     if NOT check.exists then $
        spawn, 'mkdir '+outname
     spawn, 'mv INNFL.masked* '+outname+'/'
     spawn, 'cat '+outname+'/INNFL.masked.analyze | grep "z" > tz.txt'
     readcol, 'tz.txt', zfit, f = 'X,F'
     
     print, ''
     print, '!!! INNER EXTRACTION FIT RESULTS !!!'
     print, ''
     spawn, 'cat '+outname+'/INNFL.masked.analyze'
     print, ''
     print, ''
     print, 'Input redshift : ', zfit
     plot1res, outname+'/INNFL.masked.masked.bestspec0', outname+'/INNFL.masked.masked.bestspec1', zfit
     wait, 2.0
     
     pyspecfitonetrace, id[hits[0]], pa[hits[0]], $
                        'INTUP.masked.masked', zfit, 0.1, $
                        ncores = 10
     spawn, 'mv INTUP.masked* '+outname+'/'
     spawn, 'cat '+ outname+'/INNFL.masked.analyze | grep "z"'
     spawn, 'cat '+ outname+'/INTUP.masked.analyze | grep "z"'
     print, 'Input redshift : ', zfit
     wait, 2.0
     
     pyspecfitonetrace, id[hits[0]], pa[hits[0]], $
                        'INTDN.masked', zfit, 0.1, $
                        ncores = 10
     spawn, 'mv INTDN.masked* '+outname+'/'
     spawn, 'cat '+ outname+'/INNFL.masked.analyze | grep "z"'
     spawn, 'cat '+ outname+'/INTUP.masked.analyze | grep "z"'
     spawn, 'cat '+ outname+'/INTDN.masked.analyze | grep "z"'
     print, 'Input redshift : ', zfit
     wait, 2.0
     
     pyspecfitonetrace, id[hits[0]], pa[hits[0]], $
                        'OUTUP.masked', zfit, 0.1, $
                        ncores = 10
     spawn, 'mv OUTUP.masked* '+outname+'/'
     spawn, 'cat '+ outname+'/INNFL.masked.analyze | grep "z"'
     spawn, 'cat '+ outname+'/INTUP.masked.analyze | grep "z"'
     spawn, 'cat '+ outname+'/INTDN.masked.analyze | grep "z"'
     spawn, 'cat '+ outname+'/OUTUP.masked.analyze | grep "z"'
     print, 'Input redshift : ', zfit
     wait, 2.0
     
     pyspecfitonetrace, id[hits[0]], pa[hits[0]], $
                        'OUTDN.masked', zfit, 0.1, $
                        ncores = 10
     spawn, 'mv OUTDN.masked* '+outname+'/'
     spawn, 'cat '+ outname+'/INNFL.masked.analyze | grep "z"'
     spawn, 'cat '+ outname+'/INTUP.masked.analyze | grep "z"'
     spawn, 'cat '+ outname+'/INTDN.masked.analyze | grep "z"'
     spawn, 'cat '+ outname+'/OUTUP.masked.analyze | grep "z"'
     spawn, 'cat '+ outname+'/OUTDN.masked.analyze | grep "z"'
     print, 'Input redshift : ', zfit
  

     print, ''
     print, ' >>> ALL DONE! <<< '
     print, ''
       
  endif else begin

     print, ''
     print, ' >>> OBJECT LACKS COVERAGE IN BOTH GRISMS AT THIS PA | NO FITTING POSSIBLE <<< '
     print, ''

  endelse

end
