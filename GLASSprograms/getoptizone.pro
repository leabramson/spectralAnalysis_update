function getoptizone, image, varmap, $
                      CLINEUP = clineUp, $
                      CLINEDN = clineDn

  s = size(image, /dim)
  length = s[0]
  height = s[1]

  ntogo  = height - max(clineUp)
  results = fltarr(ntogo)

  ty = findgen(height)
  for ii = 0, ntogo - 1 do begin

     apstopHi = clineUp + ii
     apstopLo = clineDn - ii

     tmp = fltarr(length)
     for jj = 0, length - 1 do begin
        tap = where((ty ge clineUp[jj] AND ty lt apstopHi[jj]) $
                    OR $
                    (ty le clineDn[jj] AND ty lt apstopLo[jj]), ntap)  

        tmp[jj] = total(image[jj,tap]) / sqrt(total(varmap[jj,tap]))
     endfor

     results[ii] = total(tmp)
     
  endfor

  foo = max(results, optiAp)
  
;  plot, findgen(ntogo), results, $
;        psym = 1, /ynoz
;  oplot, [optiap], [results[optiAp]], $
;         psym = 1, col = 255, symsize = 2
;  stop
  
  RETURN, optiAp
end
