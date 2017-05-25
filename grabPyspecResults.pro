function grabPyspecResults, filename;, $
;                            OIII = oiii

  if keyword_set(OIII) then oiii = 1 else oiii = 0

  spawn, 'cat '+filename+' | grep "z " > tmpz.dat'
  spawn, 'cat '+filename+' | grep "age" > tmpAge.dat'
  spawn, 'cat '+filename+' | grep "metal" > tmpMetal.dat'
  spawn, 'cat '+filename+' | grep "HpsEW" > tmpEW.dat'
  if oiii then $
     spawn, 'cat '+filename+' | grep "OIIIpsEW" > tmpEW_OIII.dat'
  
  readcol, 'tmpz.dat'    , zFit    , zMin    , zMax    , f = 'X,F,F,F'
  readcol, 'tmpAge.dat'  , ageFit  , ageMin  , ageMax  , f = 'X,F,F,F'
  readcol, 'tmpMetal.dat', metalFit, metalMin, metalMax, f = 'X,F,F,F'
  readcol, 'tmpEW.dat'   , EWFit   , EWMin   , EWMax   , f = 'X,F,F,F'
;  if oiii then $
  readcol, 'tmpEW_OIII.dat', EWFit_OIII, EWMin_OIII, EWMax_OIII, f = 'X,F,F,F'

;  if NOT oiii then $
;     savedata = {ZFIT: zFit[0], ZLO: zMin[0], ZHI: zMax[0], $
;                 AGEFIT: ageFit[0], AGELO: ageMin[0], AGEHI: ageMax[0], $
;                 EWFIT:EWFit[0], EWLO: EWmin[0], EWHI: EWmax[0]} $
;  else $
  savedata = {ZFIT: zFit[0], ZLO: zMin[0], ZHI: zMax[0], $
              AGEFIT: ageFit[0], AGELO: ageMin[0], AGEHI: ageMax[0], $
              METALFIT: metalFit[0], METALLO: metalMin[0], METALHI: metalMax[0], $
              EWFIT:EWFit[0], EWLO: EWmin[0], EWHI: EWmax[0], $
              EWFIT_OIII:EWFit_OIII[0], EWLO_OIII: EWmin_OIII[0], EWHI_OIII: EWmax_OIII[0]}
  
  RETURN, savedata
end

  
