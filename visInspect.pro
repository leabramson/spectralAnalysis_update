pro visInspect, inlist, outlist

  ;; Visually inspect a bunch of extracted spectra to see which have
  ;; enough good data to be worth PYSPECFITting

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
  print, 'Displaying for inspection ... '
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

        print, ' >>> Now inspecting : '+tfiles[ii]+' (B) and '+tfiles[ii+1]+' (R)'
        print, ''
        
        bdata = mrdfits(tfiles[ii], 1, /silent)
        bcontam = mrdfits(bdata.FILENAME, 'contam', /silent)
        
        rdata = mrdfits(tfiles[ii+1], 1, /silent)
        rcontam = mrdfits(rdata.FILENAME, 'contam', /silent)
        
        plot, bdata.LAMBDA, bdata.TRACE, /xsty, yran = [0,60], $
              ytitle = 'pix', xtitle = greek('lambda')+' ['+textoIDL('\AA')+']', $
              pos = [0.15,0.15,0.55,0.55]
        cgimage, bdata.SPEC2D, /over, stretch = 2
        cgtext, strmid(bdata.filename, 34, /rev), col = 'ff5500'x, charsize = 1, $
                /norm, 0.20, 0.20
        oplot, bdata.LAMBDA, bdata.INNERUP, col = 255, thick = 1
        oplot, bdata.LAMBDA, bdata.INNERDN, col = 255, thick = 1
        oplot, bdata.LAMBDA, bdata.INTERUP, col = '00aa00'x, thick = 1
        oplot, bdata.LAMBDA, bdata.INTERDN, col = '00aa00'x, thick = 1
        oplot, bdata.LAMBDA, bdata.OUTERUP, col = 'ff5500'x, thick = 1
        oplot, bdata.LAMBDA, bdata.OUTERDN, col = 'ff5500'x, thick = 1
        
        plot, bdata.LAMBDA, bdata.TRACE, /xsty, yran = [0,60], $
              xtickname = replicate(' ', 60), ytickname = replicate(' ', 60), $
              pos = [0.15,0.55,0.55,0.95], /noer
        cgimage, bcontam, /over, stretch = 2
        oplot, bdata.LAMBDA, bdata.INNERUP, col = 255, thick = 1
        oplot, bdata.LAMBDA, bdata.INNERDN, col = 255, thick = 1
        oplot, bdata.LAMBDA, bdata.INTERUP, col = '00aa00'x, thick = 1
        oplot, bdata.LAMBDA, bdata.INTERDN, col = '00aa00'x, thick = 1
        oplot, bdata.LAMBDA, bdata.OUTERUP, col = 'ff5500'x, thick = 1
        oplot, bdata.LAMBDA, bdata.OUTERDN, col = 'ff5500'x, thick = 1
        
        plot, rdata.LAMBDA, rdata.TRACE, /xsty, yran = [0,60], $
              ytickname = replicate(' ', 60), $
              xtitle = greek('lambda')+' ['+textoIDL('\AA')+']', $
              pos = [0.55,0.15,0.95,0.55], /noer
        cgimage, rdata.SPEC2D, /over, stretch = 2
        cgtext, strmid(rdata.filename, 34, /rev), col = long('0055ff'x), charsize = 1, $
                /norm, 0.60, 0.20
        oplot, rdata.LAMBDA, rdata.INNERUP, col = 255, thick = 1
        oplot, rdata.LAMBDA, rdata.INNERDN, col = 255, thick = 1
        oplot, rdata.LAMBDA, rdata.INTERUP, col = '00aa00'x, thick = 1
        oplot, rdata.LAMBDA, rdata.INTERDN, col = '00aa00'x, thick = 1
        oplot, rdata.LAMBDA, rdata.OUTERUP, col = 'ff5500'x, thick = 1
        oplot, rdata.LAMBDA, rdata.OUTERDN, col = 'ff5500'x, thick = 1
        
        plot, rdata.LAMBDA, rdata.TRACE, /xsty, yran = [0,60], $
              xtickname = replicate(' ', 60), ytickname = replicate(' ', 60), $
              pos = [0.55,0.55,0.95,0.95], /noer
        cgimage, rcontam, /over, stretch = 2
        oplot, rdata.LAMBDA, rdata.INNERUP, col = 255, thick = 1
        oplot, rdata.LAMBDA, rdata.INNERDN, col = 255, thick = 1
        oplot, rdata.LAMBDA, rdata.INTERUP, col = '00aa00'x, thick = 1
        oplot, rdata.LAMBDA, rdata.INTERDN, col = '00aa00'x, thick = 1
        oplot, rdata.LAMBDA, rdata.OUTERUP, col = 'ff5500'x, thick = 1
        oplot, rdata.LAMBDA, rdata.OUTERDN, col = 'ff5500'x, thick = 1

        print, tfiles[ii], ' ', field[counter]+'_'+string(idNo[counter], f = '(I05)')
        print, ''
        decision = ''
        read, decision, prompt = 'Fit this object? [Y/n] '
        print, ''

        if decision eq 'y' or decision eq '' or decision eq 'Y' then begin

           printf, 2, tfiles[ii], ' ', tfiles[ii+1]
           keep = [keep, dirname+'_'+string(objno, f = '(I05)')]
                     
           printf, 1, field[counter]+' ', idNo[counter], $
                   ra[counter], dec[counter], $
                   z[counter], zq[counter], $
                   pz[counter], mag[counter], strmid(repstr(tfiles[ii], '.fits', ''), 2, 1, /rev)
                      
        endif

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
;;
;;

pro doit

  visInspect, $
     'sourceLists/1.0_1.8_14.0_21.8_c99.9_pacut0_ALLSRCS.list', $
     '1.0_1.8_14.0_21.8_c99.9_pacut0_ALLSRCS_INSPECTED.list'

;  visInspect, $
;     'sourceLists/1.0_1.8_14.0_21.8_c99.9_pacut0_ALLSRCS.list', $
;     'foo.list'
end
