pro join2cats, masterIn, otherCatIn, output, $
               STRIPTAG = striptag, $
               POINTING = pointing, $
               DPOSNAME = dposname
  
  master = mrdfits(masterIn, 1)
  n1    = n_elements(master)
  mra  = master.X_WORLD
  mdec = master.Y_WORLD

  slave  = mrdfits(otherCatIn, 1)
  n2 = n_elements(slave)

  newtags = tag_names(slave)
  ntags = n_elements(newtags)
  if total(newtags eq 'ALPHA_J2000') eq 0 then begin
     sra  = slave.RA
     sdec = slave.DEC
  endif else begin
     sra  = slave.ALPHA_J2000  ;; Kuang puts his photozs this way
     sdec = slave.DELTA_J2000  
  endelse
     
  spherematch, mra, mdec, sra, sdec, 0.2d/3600, $
               minds, sinds, len
  if n_elements(striptag) gt 0 then begin
     check = []    
     for ii = 0, n_elements(striptag) - 1 do $
        check = [check, where(newtags eq striptag[ii])]
;        slave = mod_struct(slave, striptag[ii], /delete)
  endif else check = !VALUES.F_NAN
;  newtags = tag_names(slave)
;  ntags = n_elements(newtags)

  tagsdone = 0
  for ii = 0, ntags - 1 do begin

     if ~total(check eq ii) then begin
     
        tmpdata = slave.(ii)
        type    = size(tmpdata, /tname)
        case type of
           'STRING' : tmpVal = 'string(-99)'
           'FLOAT'  : tmpVal = '-99.'
           'LONG'   : tmpVal = '-99L'
           'INT'    : tmpVal = '-99'
           'DOUBLE' : tmpVal = '-99.d'
        endcase
        
        if tagsdone eq 0 then begin
           com = 'tmpstruct = {'+newtags[ii]+': '+string(tmpVal)+'}'
           foo = execute(com)
           com = 'tmpstruct = replicate(tmpstruct, n1)'
           foo = execute(com)           
        endif else begin
           com = 'tmpstruct1 = {'+newtags[ii]+': '+string(tmpVal)+'}'
           foo = execute(com)
           com = 'tmpstruct1 = replicate(tmpstruct1, n1)'
           foo = execute(com)
           tmpstruct = struct_addtags(tmpstruct, tmpstruct1)
        endelse
        
        tmpstruct[minds].(tagsdone) = tmpdata[sinds]

        tagsdone++
        
     endif

  endfor

  ;; Add a column for source offsets
  if KEYWORD_SET(DPOSNAME) then begin
     com = 'tmpstruct1 = {'+dposname+': -99.d}'
     foo = execute(com)
     com = 'tmpstruct1 = replicate(tmpstruct1, n1)'
     foo = execute(com)
     tmpstruct = struct_addtags(tmpstruct, tmpstruct1)
     tmpstruct[minds].(tagsdone) = len[sinds] * 3600
  endif
  
  master = struct_addtags(master, tmpstruct)

  if keyword_set(POINTING) then begin
     tmp = {POINTING: pointing}
     tmp = replicate(tmp, n_elements(master))
     master = struct_addtags(master, tmp)
  endif

  mwrfits, master, output, /create

end
