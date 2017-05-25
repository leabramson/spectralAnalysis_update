pro pyspecfitfitonefield

  check = file_info('fitspec.py')
  if NOT check.exists then $
     spawn, 'ln -s ../../fitspec.py .'
  
;  readcol, 'fileBreakdownFor1Dextraction.list', $
;           ID, PA, GRISM, ZP, Z, $
;           f = 'X,A,A,A,F,F'

  sid = id[sort(id)]
  sz  = z[sort(ID)]
  spz = zp[sort(ID)]

  uid = sid[uniq(sid)]
  uz  = sz[uniq(sid)]
  upz = spz[uniq(sid)]
  ngals = n_elements(uid)

  for ii = 0, ngals - 1 do begin

     ;; Find a "redshift"
     tpz = upz[ii]
     tz  = uz[ii]

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
        
     cd, uid[ii]+'_1Dextractions'
     test_probe1 = uid[ii]+'_1_B_INNFL.pyspec'
     test_probe2 = uid[ii]+'_2_B_INNFL.pyspec'
     readcol, test_probe1, flx1, msk1, f = 'X,F,X,I'
     readcol, test_probe2, flx2, msk2, f = 'X,F,X,I'

     flag1 = total(flx1(where(msk1)))
     flag2 = total(flx2(where(msk2)))

;     stop
     
     if flag1 gt total(msk1) then begin
        pyspecfitonetrace, uid[ii], 1, 'INNFL', zstart, zsran, $
                           ncores = 10
        outname = uid[ii]+'_1_pyspecfitResults'
        check = file_info(outname)
        if NOT check.exists then $
           spawn, 'mkdir '+outname
        spawn, 'mv INNFL* '+outname+'/'
        spawn, 'cat '+outname+'/INNFL.analyze | grep "z" > tz.txt'
        readcol, 'tz.txt', zfit, f = 'X,F'

        print, ''
        print, '!!! INNER EXTRACTION FIT RESULTS !!!'
        print, ''
        spawn, 'cat '+outname+'/INNFL.analyze'
        print, ''
        print, ''
        print, 'Input redshift : ', zfit
        plot1res, outname+'/INNFL.bestspec0', outname+'/INNFL.bestspec1', zfit
        wait, 2.0
        
        pyspecfitonetrace, uid[ii], 1, 'INTUP', zfit, 0.1, $
                           ncores = 10
        spawn, 'mv INTUP* '+outname+'/'
        spawn, 'cat '+ outname+'/INNFL.analyze | grep "z"'
        spawn, 'cat '+ outname+'/INTUP.analyze | grep "z"'
        print, 'Input redshift : ', zfit
        wait, 2.0
        
        pyspecfitonetrace, uid[ii], 1, 'INTDN', zfit, 0.1, $
                           ncores = 10
        spawn, 'mv INTDN* '+outname+'/'
        spawn, 'cat '+ outname+'/INNFL.analyze | grep "z"'
        spawn, 'cat '+ outname+'/INTUP.analyze | grep "z"'
        spawn, 'cat '+ outname+'/INTDN.analyze | grep "z"'
        print, 'Input redshift : ', zfit
        wait, 2.0
        
        pyspecfitonetrace, uid[ii], 1, 'OUTUP', zfit, 0.1, $
                           ncores = 10
        spawn, 'mv OUTUP* '+outname+'/'
        spawn, 'cat '+ outname+'/INNFL.analyze | grep "z"'
        spawn, 'cat '+ outname+'/INTUP.analyze | grep "z"'
        spawn, 'cat '+ outname+'/INTDN.analyze | grep "z"'
        spawn, 'cat '+ outname+'/OUTUP.analyze | grep "z"'
        print, 'Input redshift : ', zfit
        wait, 2.0
        
        pyspecfitonetrace, uid[ii], 1, 'OUTDN', zfit, 0.1, $
                           ncores = 10
        spawn, 'mv OUTDN* '+outname+'/'
                spawn, 'cat '+ outname+'/INNFL.analyze | grep "z"'
        spawn, 'cat '+ outname+'/INTUP.analyze | grep "z"'
        spawn, 'cat '+ outname+'/INTDN.analyze | grep "z"'
        spawn, 'cat '+ outname+'/OUTUP.analyze | grep "z"'
        spawn, 'cat '+ outname+'/OUTDN.analyze | grep "z"'
        print, 'Input redshift : ', zfit
        wait, 2.0
        
     endif

     if flag2 gt total(msk2) then begin
        pyspecfitonetrace, uid[ii], 2, 'INNFL', zstart, zsran, $
                           ncores = 10
        outname = uid[ii]+'_2_pyspecfitResults'
        check = file_info(outname)
        if NOT check.exists then $
           spawn, 'mkdir '+outname
        spawn, 'mv INNFL* '+outname+'/'
        spawn, 'cat '+outname+'/INNFL.analyze | grep "z" > tz.txt'
        readcol, 'tz.txt', zfit, f = 'X,F'

        print, ''
        print, '!!! INNER EXTRACTION FIT RESULTS !!!'
        print, ''
        spawn, 'cat '+ outname+'/INNFL.analyze'
        print, ''
        print, ''
        print, 'Input redshift : ', zfit        
        plot1res, outname+'/INNFL.bestspec0', outname+'/INNFL.bestspec1', zfit
        wait, 2.0
        
        pyspecfitonetrace, uid[ii], 2, 'INTUP', zfit, 0.1, $
                           ncores = 10
        spawn, 'mv INTUP* '+outname+'/'
        spawn, 'cat '+ outname+'/INNFL.analyze | grep "z"'
        spawn, 'cat '+ outname+'/INTUP.analyze | grep "z"'
        print, 'Input redshift : ', zfit
        wait, 2.0
        
        pyspecfitonetrace, uid[ii], 2, 'INTDN', zfit, 0.1, $
                           ncores = 10
        spawn, 'mv INTDN* '+outname
        spawn, 'cat '+ outname+'/INNFL.analyze | grep "z"'
        spawn, 'cat '+ outname+'/INTUP.analyze | grep "z"'
        spawn, 'cat '+ outname+'/INTDN.analyze | grep "z"'
        print, 'Input redshift : ', zfit
        wait, 2.0
        
        pyspecfitonetrace, uid[ii], 2, 'OUTUP', zfit, 0.1, $
                           ncores = 10
        spawn, 'mv OUTUP* '+outname+'/'
        spawn, 'cat '+ outname+'/INNFL.analyze | grep "z"'
        spawn, 'cat '+ outname+'/INTUP.analyze | grep "z"'
        spawn, 'cat '+ outname+'/INTDN.analyze | grep "z"'
        spawn, 'cat '+ outname+'/OUTUP.analyze | grep "z"'
        print, 'Input redshift : ', zfit
        wait, 2.0
        
        pyspecfitonetrace, uid[ii], 2, 'OUTDN', zfit, 0.1, $
                           ncores = 10
        spawn, 'mv OUTDN* '+outname+'/'
        spawn, 'cat '+ outname+'/INNFL.analyze | grep "z"'
        spawn, 'cat '+ outname+'/INTUP.analyze | grep "z"'
        spawn, 'cat '+ outname+'/INTDN.analyze | grep "z"'
        spawn, 'cat '+ outname+'/OUTUP.analyze | grep "z"'
        spawn, 'cat '+ outname+'/OUTDN.analyze | grep "z"'
        print, 'Input redshift : ', zfit
        wait, 2.0
        
     endif
     cd, '..'
     
  endfor

end
