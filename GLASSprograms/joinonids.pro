function joinonids, master, slave, $
                    MASTER_ID_TAG = master_ID_TAG, $
                    SLAVE_ID_TAG = slave_id_tag, $
                    STRIPTAG = striptag, $
                    POINTING = pointing

  if NOT keyword_set(MASTER_ID_TAG) then $
     master_id_tag = 'ID_TAKA'

  if NOT keyword_set(slave_ID_TAG) then $
     slave_id_tag = 'ID_TAKA'
  
  n1 = n_elements(master)
  n2 = n_elements(slave)

  mtags = tag_names(master)
  newtags = tag_names(slave)
  ntags = n_elements(newtags)

  midtag = (where(mtags eq master_id_tag))[0]
  sidtag = (where(newtags eq slave_id_tag))[0]

  sinds = []
  minds = []
  for ii = 0, n_elements(master) - 1 do begin
     hit = where(slave.(sidtag) eq master[ii].(midtag), nhit)
     if nhit eq 1 then begin
        hit = hit[0]
        minds = [minds, ii]
        sinds = [sinds, hit]
     endif
  endfor

  if n_elements(striptag) gt 0 then begin
     check = []    
     for ii = 0, n_elements(striptag) - 1 do $
        check = [check, where(newtags eq striptag[ii])]
  endif else check = !VALUES.F_NAN
  
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

  joined = struct_addtags(master, tmpstruct)

  if keyword_set(POINTING) then begin
     tmp = {POINTING: pointing}
     tmp = replicate(tmp, n_elements(master))
     joined = struct_addtags(joined, tmp)
  endif
  
  RETURN, joined
end

