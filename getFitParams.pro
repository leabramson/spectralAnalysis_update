;; Extract one set of parameters for a list of objects

pro getFitParams, dirlist, $
                  PARAM = param, $
                  OUTPUT = output

  ;; Param can be 'age', 'HpsEW', 'OIII' as of May 5, 2016
  
  if NOT keyword_set(PARAM)  then param  = 'age'
  if NOT keyword_set(OUTPUT) then output = param+'_pyspec_results.fits'

  readcol, dirlist, dirs, f = 'A'
  ndirs = n_elements(dirs)

  savedata = {FILE: string(0), $
              PARAMETER: string(0), $
              VALUES: fltarr(3,2), $
              HI_ERR: fltarr(3,2),  $
              LO_ERR: fltarr(3,2)}
  savedata = replicate(savedata, ndirs)

  for ii = 0, ndirs - 1 do begin

     spawn, 'cat '+dirs[ii]+'*1.analyze | grep '+param+' > tmpParams1.list'
     spawn, 'cat '+dirs[ii]+'*2.analyze | grep '+param+' > tmpParams2.list'
     readcol, 'tmpParams1.list', paramFit1, paramLo1, paramHi1, f = 'X,F,F,F'
     readcol, 'tmpParams2.list', paramFit2, paramLo2, paramHi2, f = 'X,F,F,F'

     if param eq 'age' then begin
        p1Ehi = 0.434 * (paramHi1 - paramFit1) / paramFit1
        p1Elo = 0.434 * (paramLo1 - paramFit1) / paramFit1 * (-1)
        p1    = alog10(paramFit1)
        p2Ehi = 0.434 * (paramHi2 - paramFit2) / paramFit2
        p2Elo = 0.434 * (paramLo2 - paramFit2) / paramFit2 * (-1)
        p2    = alog10(paramFit2)
     endif else begin
        p1Ehi = paramHi1 - paramFit1
        p1Elo = paramFit1 - paramLo1
        p1    = paramFit1
        p2Ehi = paramHi2 - paramFit2
        p2Elo = paramFit2 - paramLo2
        p2    = paramFit2
     endelse
        
;     params[ii,*,0,0] = p1
;     params[ii,*,1,0] = p1Elo
;     params[ii,*,2,0] = p1Ehi
;     params[ii,*,0,1] = p2
;     params[ii,*,1,1] = p2Elo
;     params[ii,*,2,1] = p2Ehi

     savedata[ii].FILE        = dirs[ii]
     savedata[ii].PARAMETER   = param
     savedata[ii].VALUES[*,0] = p1
     savedata[ii].VALUES[*,1] = p2
     savedata[ii].HI_ERR[*,0] = p1Ehi
     savedata[ii].HI_ERR[*,1] = p2Ehi
     savedata[ii].LO_ERR[*,0] = p1Elo
     savedata[ii].LO_ERR[*,1] = p2Elo

  endfor

  mwrfits, savedata, output, /create
  
end
