function eazylookup, bandpass

  readcol, 'eazyFilterDef.dat', ef, rf, $
           f = 'I,I', comment = '#', /quick, /silent

  eazyfilter = ef[where(rf eq bandpass)]
  
  RETURN, eazyfilter
end
