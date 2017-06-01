function eazylookup, bandpass

  readcol, '~/LocalData/Abramson/eazyFilterDef.dat', ef, rf, $
           f = 'I,I', comment = '#', /quick, /silent

  eazyfilter = ef[where(rf eq bandpass)]
  
  RETURN, eazyfilter
end

;;
;;
;;

pro printphotforpyspec, photfile, $
                        REGION = region, $
                        OUTPUT = output, $
                        PRINT_LIMS = print_lims

  if NOT keyword_set(print_lims) then $
     print_lims = 0 $
  else $
     print_lims = 1

  data = mrdfits(photfile, 1)
  tags = tag_names(data)

  region = strupcase(region)
  case region of
     'OPTFL'  : begin
        sed = data.SED_OPTI
        err = data.SED_OPTI_ERR
     end
     'INNFL'   : begin
        sed = data.SED_INNER
        err = data.SED_INNER_ERR
     end
     'INTUP': begin
        sed = data.SED_INTER_UP
        err = data.SED_INTER_UP_ERR
     end
     'INTDN': begin
        sed = data.SED_INTER_DN
        err = data.SED_INTER_DN_ERR
     end
     'OUTUP': begin
        sed = data.SED_OUTER_UP
        err = data.SED_OUTER_UP_ERR
     end
     'OUTDN': begin
        sed = data.SED_OUTER_DN
        err = data.SED_OUTER_DN_ERR
     end
  endcase

  filts = data.FILTER
  nfilts = n_elements(filts)

  get_lun, lun
  openw, lun, output, width = 512
  for ii = 0, nfilts - 1 do begin

     ef = eazylookup(filts[ii])

     if NOT print_lims then begin  ;;  can only do no limits for now
        if sed[ii] gt 0 then $
           printf, lun, f = '(%"%i %f %f")', $
                   ef, sed[ii], err[ii]
     endif else begin
        if sed[ii] gt 0 then $
           printf, lun, f = '(%"%i %f %f")', $
                   ef, sed[ii], err[ii] $
        else $
           stop
     endelse

  endfor 
  close, lun
  free_lun, lun
  
  spawn, 'cat '+output  
  
end

pro test
  printphotforpyspec, 'MACS1149/00753_1_B_resolvedSED.fits', $
                      region = 'INNFL', output = 'test.sed'
  printphotforpyspec, 'MACS0744/00660_2_B_resolvedSED.fits', $
                      region = 'INNFL', output = 'test.sed'
end
