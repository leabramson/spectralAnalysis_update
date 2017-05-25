pro writeonedirforpyspecfit, field

  organizeonefield, field

  readcol, field+'/fileBreakdownFor1Dextraction.list', $
           fname, ID, PA, grism, comment = '#', $
           f = 'A,A,A,A'

  sID = ID[sort(ID)]
  uID = sID[UNIQ(sID)]
  ngals = n_elements(uID)

  for ii = 0, ngals - 1 do begin

     print, ''
     print, 'Writing extraction for ID : '+uId[ii]+' in field : '+field
     print, ''

     outdir = field+'/'+uID[ii]+'_1Dextractions'
     check = file_info(outdir)
     if NOT check.exists then $
        spawn, 'mkdir '+outdir

     hits = where(ID eq uID[ii], nhits)

     tfiles  = fname[hits]
     tPAs    = PA[hits]
     tgrisms = grism[hits]

     for jj = 0, nhits - 1 do begin

        if tgrisms[jj] eq 'B' then $
           CHIP = 'blue' $
        else $
           CHIP = 'red'

        outbase = repstr(tfiles[jj], '.fits')

        optiSuff    = '_OPTFL.pyspec'
        innerSuff   = '_INNFL.pyspec'
        interUpSuff = '_INTUP.pyspec'
        interDnSuff = '_INTDN.pyspec'
        outerUpSuff = '_OUTUP.pyspec'
        outerDnSuff = '_OUTDN.pyspec'

        input = field+'/'+tfiles[jj]
        
        ;; Do optimal aperture -- serves as basis for PYSPECFIT modeling
        output = outdir+'/'+outbase+optisuff
        writeForPyspecfit, input, output, $
                           CHIP = CHIP, $
                           REGION = 'optiTOT';, $
;                           MASKMORE = maskmore

        ;; Do other extractions
        output = outdir+'/'+outbase+innerSuff
        writeForPyspecfit, input, output, $
                           CHIP = CHIP, $
                           REGION = 'inner', $
                           MASKMORE = 0

        output = outdir+'/'+outbase+interUpSuff
        writeForPyspecfit, input, output, $
                           CHIP = CHIP, $
                           REGION = 'interUP', $
                           MASKMORE = 0

        output = outdir+'/'+outbase+interDnSuff
        writeForPyspecfit, input, output, $
                           CHIP = CHIP, $
                           REGION = 'interDN', $
                           MASKMORE = 0

        output = outdir+'/'+outbase+outerUpSuff
        writeForPyspecfit, input, output, $
                           CHIP = CHIP, $
                           REGION = 'outerUP', $
                           MASKMORE = 0

        output = outdir+'/'+outbase+outerDnSuff
        writeForPyspecfit, input, output, $
                           CHIP = CHIP, $
                           REGION = 'outerDN', $
                           MASKMORE = 0

        
     endfor
     
     
  endfor
  
end

;;
;;
;;

pro writeall

  writeonedirforpyspecfit, 'RXJC2248'
  writeonedirforpyspecfit, 'RXJC1347'
;  writeonedirforpyspecfit, 'MACS0416'
  writeonedirforpyspecfit, 'MACS0744'
  writeonedirforpyspecfit, 'MACS0717'
  writeonedirforpyspecfit, 'MACS1149'
  writeonedirforpyspecfit, 'MACS1423'
  writeonedirforpyspecfit, 'MACS2129'
  writeonedirforpyspecfit, 'ABEL2744'
;  writeonedirforpyspecfit, 'ABEL0370'
  
end
