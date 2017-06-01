pro quantInspect, inlist, outlist, $
                  CONTAM_LEVEL = contam_level, $
                  CONTAM_THRESH = contam_thresh, $
                  ZONE = zone

  ;; Quantitatively cull a list of objects based on their
  ;; contamination levels.

  if NOT keyword_set(CONTAM_LEVEL) then $ ;; contamination model is .ge. 0.1 * target 2D spec.
     contam_level = 0.1 
  if NOT keyword_set(CONTAM_THRESH) then $ ;; discard if this frac. of pixels in ZONE is above this number
     contam_thresh = 0.2
  if NOT keyword_set(ZONE) then zone = 'OUT'
  
  print, ''
  print, 'Reading file list ... '

  readcol, inlist, $
           field, $
           idNo, $
           ra, $
           dec, $
           z, $
           zq, $
           pz, $
           mag, $
           f = 'A,I,D,D,F,F,F,F', /silent
  nfiles = n_elements(idNo)
  
  print, '... '+string(nfiles, f = '(I04)')+' files found.'
  print, ''
  print, 'Analyzing contamination levels ... '
  print, ''

  close, 1
  openw, 1, outlist, width = 256
  printf, 1, '#FIELD GLASS_ID RA DEC Z ZQ KZ_PHOT MAG PA(S)'

  close, 2
  openw, 2, 'FITS_FILES_ONLY_'+outlist
  
  counter = 0
  keep = [string(0)]
  window, 0, xsize = 1200, ysize = 1200, retain = 2
  while counter le nfiles - 1 do begin

     dirname = strmid(field[counter], 0, 8) ; directory of files
     objno   = string(idNo[counter], f = '(I05)')
     spawn, 'ls '+dirname+'/'+objno+'_?_?.fits > tmp.list'
     readcol, 'tmp.list', tfiles, f = 'A', /silent
     ntfiles = n_elements(tfiles)

     ii = 0
     while ii le ntfiles - 1 do begin

        print, f = '(%" >>> Now analyzing : %s (B) and %s (R)\n")', $
               tfiles[ii], tfiles[ii+1]        
        
;        bdata = mrdfits(tfiles[ii], 1, /silent)
;        bsci    = mrdfits(bdata.FILENAME, 'sci', /silent)
;        bcontam = mrdfits(bdata.FILENAME, 'contam', /silent)
        
;        rdata = mrdfits(tfiles[ii+1], 1, /silent)
;        rsci    = mrdfits(rdata.FILENAME, 'sci', /silent)
;        rcontam = mrdfits(rdata.FILENAME, 'contam', /silent)

        bc = quantcontam(tfiles[ii]  , CFRACT = CONTAM_LEVEL)
        rc = quantcontam(tfiles[ii+1], CFRACT = CONTAM_LEVEL)

        case zone of
           'TOT': begin
              if bc.TOTAL le contam_thresh $
                 AND $
                 rc.TOTAL le contam_thresh $
              then $
                 keepflag = 1 $
              else $
                 keepflag = 0
           end
           'OPT': begin
              if bc.OPT_ZONE le contam_thresh $
                 AND $
                 rc.OPT_ZONE le contam_thresh $
              then $
                 keepflag = 1 $
              else $
                 keepflag = 0
           end
           'OUT': begin
              if bc.OUT_ZONE le contam_thresh $
                 AND $
                 rc.OUT_ZONE le contam_thresh $
              then $
                 keepflag = 1 $
              else $
                 keepflag = 0
           end
           'INT': begin
              if bc.INT_ZONE le contam_thresh $
                 AND $
                 rc.INT_ZONE le contam_thresh $
              then $
                 keepflag = 1 $
              else $
                 keepflag = 0
           end
           'INN': begin
              if bc.INN_ZONE le contam_thresh $
                 AND $
                 rc.INN_ZONE le contam_thresh $
              then $
                 keepflag = 1 $
              else $
                 keepflag = 0
           end
        endcase

        if keepflag then begin
           if bc.INN_NE le contam_thresh $
              AND $
              rc.INN_NE le contam_thresh $
           then $
              keepflag = 1 $
           else $
              keepflag = 0
        endif
        
        if keepflag then begin

           print, f = '(%"Retaining\n")'
           
           printf, 2, tfiles[ii], ' ', tfiles[ii+1]
           keep = [keep, dirname+'_'+string(objno, f = '(I05)')]
                     
           printf, 1, field[counter]+' ', idNo[counter], $
                   ra[counter], dec[counter], $
                   z[counter], zq[counter], $
                   pz[counter], mag[counter], $
                   strmid(repstr(tfiles[ii], '.fits', ''), 2, 1, /rev)
                      
        endif else print, f = '(%"Discarding\n")'
        
        ii+=2

     endwhile

     counter++
     
  endwhile

  close, 1
  close, 2

  keep = keep[1:*]
  nkeep = n_elements(keep)

  keep = keep[sort(keep)]
  nobj = n_elements(UNIQ(keep))

  print, ''
  print, 'All done. '+string(nkeep)+' files retained for '+string(nobj)+' unique objects'
  
end
;;

pro doit
  quantInspect, $
     'sourceLists/1.0_1.8_14.0_21.8_c99.9_pacut0_ALLSRCS.list', $
     '1.0_1.8_14.0_21.8_c99.9_pacut0_ALLSRCS_QUANTITATIVE_3031.list', $
     ZONE = 'OUT', $
     CONTAM_LEVEL = 0.3, $
     CONTAM_THRESH = 0.31
end
