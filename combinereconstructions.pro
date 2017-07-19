pro combinereconstructions, pattern, $
                            OUTPUT = output

  regions = ['INNFL', 'INTER', 'OUTER', $
             'INTUP', 'INTDN', 'OUTUP', 'OUTDN', $
             'OPTFL']
  nregions = n_elements(regions)

  test = mrdfits(repstr(pattern, 'XXXXX', regions[-1]), 1, /silent)
  niter = n_elements(test.TIME)
  
  savedata = {REGION:   string(0), $
              TIME:     test.TIME, $
              REDSHIFT: test.REDSHIFT, $
              TOBS:     test.TOBS, $
              SFH:      findgen(niter,3), $
              MGH:      findgen(niter,3), $
              SSFH:     findgen(niter,3)}
  savedata = replicate(savedata, nregions)
  savedata[-1].REGION = regions[-1]
  savedata[-1].SFH    = test.SFH
  savedata[-1].MGH    = test.MGH
  savedata[-1].SSFH   = test.SSFH
  
  for ii = 0, nregions - 2 do begin
     test = mrdfits(repstr(pattern, 'XXXXX', regions[ii]), 1, /silent)
     savedata[ii].REGION = regions[ii]     
     savedata[ii].SFH    = test.SFH
     savedata[ii].MGH    = test.MGH
     savedata[ii].SSFH   = test.SSFH
  endfor

  mwrfits, savedata, output, /create

end
