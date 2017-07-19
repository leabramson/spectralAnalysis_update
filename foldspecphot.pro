pro foldspecphot, prefix, region, $
                  NOTMASKED = notmasked

  if NOT keyword_set(notmasked) then $
     notmasked = 0 $
  else $
     notmasked = 1
  
  region = strupcase(region)

  spawn, 'rm '+prefix+'_?_'+region+'*pyspec'
  spawn, 'rm '+prefix+'_?_'+region+'*lsf'
  spawn, 'rm '+prefix+'_'+region+'*sed'
  
  case region of
     'INTER': begin
        sprefix = prefix+'_X_INT??.masked'
        pprefix = prefix+'_INT??.sed'
     end
     'OUTER': begin
        sprefix = prefix+'_X_OUT??.masked'
        pprefix = prefix+'_OUT??.sed'
     end
  endcase

  bprefix = repstr(sprefix, '_X_', '_B_')
  rprefix = repstr(sprefix, '_X_', '_R_')

  spawn, 'ls '+bprefix+'.pyspec > tmpBSpectra.list'
  spawn, 'ls '+rprefix+'.pyspec > tmpRSpectra.list'

  spawn, 'ls '+bprefix+'.lsf > tmpBLSF.list'
  spawn, 'ls '+rprefix+'.lsf > tmpRLSF.list'

  spawn, 'ls '+pprefix+' > tmpPhot.list'

  readcol, 'tmpBSpectra.list', bspecs, f = 'A'
  readcol, 'tmpRSpectra.list', rspecs, f = 'A'
  readcol, 'tmpBLSF.list', blsfs, f = 'A'
  readcol, 'tmpRLSF.list', rlsfs, f = 'A'
  readcol, 'tmpPhot.list', seds, f = 'A'

  ;; Do photometry
  readcol, seds[0], l1, f1, e1, f = 'F,F,F'
  readcol, seds[1], l2, f2, e2, f = 'F,F,F'

  l = [l1, l2]
  l = l[sort(l)]
  l = l[UNIQ(l)]
  pf = fltarr(n_elements(l))
  pe = fltarr(n_elements(l))
  for ii = 0, n_elements(l) - 1 do begin     
     u1 = where(l1 eq l[ii], n1)
     u2 = where(l2 eq l[ii], n2)
     if n1 AND n2 then begin
        pf[ii] = 0.5 * (f1[u1] + f2[u2])
        pe[ii] = sqrt(e1[u1]^2 + e2[u2]^2) / sqrt(2.)
     endif else if n1 AND NOT n2 then begin
        pf[ii] = f1[u1]
        pe[ii] = e1[u1]
     endif else if n2 AND NOT n1 then begin
        pf[ii] = f2[u2]
        pe[ii] = e2[u2]
     endif
  endfor
  
  close, 1
  openw, 1, prefix+'_'+region+'.sed'
  for ii = 0, n_elements(l) - 1 do $
     printf, 1, fix(l[ii]), pf[ii], pe[ii]
  close, 1

  ii = 0
  while ii lt 2 do begin
     case ii of
        0: begin
           sfiles = bspecs
           lfiles = blsfs
           outfilt = '_B_'
        end
        1: begin
           sfiles = rspecs
           lfiles = rlsfs
           outfilt = '_R_'
        end
     endcase
     
     readcol, sfiles[0], l1, f1, e1, m1, f = 'F,F,F,F'
     readcol, sfiles[1], l2, f2, e2, m2, f = 'F,F,F,F'

     s = fltarr(n_elements(l1))
     e = fltarr(n_elements(l1))
     m = fltarr(n_elements(l1))
     
     for jj = 0, n_elements(l1) - 1 do begin
        if m1[jj] and m2[jj] then begin
           s[jj] = 0.5 * (f1[jj] + f2[jj])
           e[jj] = sqrt(e1[jj]^2 + e2[jj]^2) / sqrt(2.)
           m[jj] = 1
        endif else if m1[jj] AND NOT m2[jj] then begin
           s[jj] = f1[jj]
           e[jj] = e1[jj]
           m[jj] = 1
        endif else if m2[jj] AND NOT m1[jj] then begin
           s[jj] = f2[jj]
           e[jj] = e2[jj]
           m[jj] = 1
        endif else if NOT m1[jj] AND NOT m2[jj] then begin
           m[jj] = 0
        endif
     endfor
     
     readcol, lfiles[0], x1, k1, f = 'F,F,'
     readcol, lfiles[1], x2, k2, f = 'F,F,'

     k = 0.5 * (k1 + k2)

     close, 1
     close, 2
     openw, 1, prefix+outfilt+region+'.masked.pyspec'
     openw, 2, prefix+outfilt+region+'.masked.lsf'
     for jj = 0, n_elements(l1) - 1 do $
        printf, 1, l1[jj], s[jj], e[jj], m[jj]
     close, 1
     for jj = 0, n_elements(x1) - 1 do $
        printf, 2, x1[jj], k[jj]
     close, 2

     ii++
     
  endwhile
  
end
;foldspecphot, '00451_2', 'INTER'
