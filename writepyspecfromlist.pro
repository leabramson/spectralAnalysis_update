;; Write after culling using visual inspection via visInspec.pro

pro writepysepecfromlist, inlist
  
  readcol, inlist, bluefiles, redfiles, f = 'A,A'
  nobj = n_elements(bluefiles)

  for ii = 0, nobj - 1 do begin

     field = strmid(bluefiles[ii], 0, 8)

;     outdir   = field
     boutbase = repstr(bluefiles[ii], '.fits')
     routbase = repstr(redfiles[ii], '.fits')
     
     print, ''
     print, 'WRITING 1D SPECTRA FOR : '+bluefiles[ii], ' ', redfiles[ii]
     print, ''
     
     for jj = 0, 5 do begin

        case jj of
           0: begin
              suffix = '_OPTFL.pyspec'
              region = 'optiTOT'
           end
           1: begin
              suffix = '_INNFL.pyspec'
              region = 'inner'
           end
           2: begin
              suffix = '_INTUP.pyspec'
              region = 'interUP'
           end
           3: begin
              suffix = '_INTDN.pyspec'
              region = 'interDN'
           end
           4: begin
              suffix = '_OUTUP.pyspec'
              region = 'outerUP'
           end           
           5: begin
              suffix = '_OUTDN.pyspec'
              region = 'outerDN'
           end
        endcase
        
        boutput = boutbase+suffix
        routput = routbase+suffix

        print, ' ... to : '+boutput, ' ', routput
        
        writeForPyspecfit, bluefiles[ii], boutput, $
                           CHIP = 'blue', $
                           REGION = region
        writeForPyspecfit, bluefiles[ii], repstr(boutput, '.pyspec', '.masked.pyspec'), $
                           CHIP = 'blue', $
                           REGION = region, $
                           /maskmore

        writeForPyspecfit, redfiles[ii], routput, $
                           CHIP = 'red', $
                           REGION = region
        writeForPyspecfit, redfiles[ii], repstr(routput, '.pyspec', '.masked.pyspec'), $
                           CHIP = 'red', $
                           REGION = region, $
                           /maskmore

     endfor   
     print, ''
     
  endfor

end

;; writepysepecfromlist, 'FITS_FILES_ONLY_1.0_1.8_14.0_21.8_c99.9_pacut0_ALLSRCS_INSPECTED.list'
;; writepysepecfromlist, 'FITS_FILES_ONLY_1.0_1.8_14.0_21.8_c99.9_pacut0_ALLSRCS_QUANTITATIVE.list'
;;  writepysepecfromlist, 'FITS_FILES_ONLY_1.0_1.8_14.0_21.8_c99.9_pacut0_MACS0717_QUANTITATIVE.list'
