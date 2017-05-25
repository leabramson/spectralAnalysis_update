pro easyspecvis, bluefile, redfile, $
                 XSIZE = xsize, YSIZE = ysize

  if NOT keyword_set(XSIZE) then xsize = 600
  if NOT keyword_set(YSIZE) then ysize = 600
  
  bdata = mrdfits(bluefile, 1, /silent)
  rdata = mrdfits(redfile , 1, /silent)
  
  window, 3, xsize = xsize, ysize = ysize, retain = 2

  plot, bdata.LAMBDA, bdata.TRACE, $
        /xsty, yran = [0,60], $
        ytitle = 'pix', xtitle = greek('lambda')+' ['+textoIDL('\AA')+']', $
        pos = [0.15,0.15,0.55,0.55], charsize = 1.1
  cgimage, bdata.SPEC2D, /over, stretch = 2
  cgtext, strmid(bdata.filename, 34, /rev), col = 'ff5500'x, charsize = 1.1, $
          /norm, 0.175, 0.20, charthick = 1.5
  oplot, bdata.LAMBDA, bdata.INNERUP, col = 255, thick = 1
  oplot, bdata.LAMBDA, bdata.INNERDN, col = 255, thick = 1
  oplot, bdata.LAMBDA, bdata.INTERUP, col = '00aa00'x, thick = 1
  oplot, bdata.LAMBDA, bdata.INTERDN, col = '00aa00'x, thick = 1
  oplot, bdata.LAMBDA, bdata.OUTERUP, col = 'ff5500'x, thick = 1
  oplot, bdata.LAMBDA, bdata.OUTERDN, col = 'ff5500'x, thick = 1
  
  plot, bdata.LAMBDA, bdata.TRACE, $
        /xsty, yran = [0,60], $
        xtickname = replicate(' ', 60), ytickname = replicate(' ', 60), $
        pos = [0.15,0.55,0.55,0.95], /noer
  bcontam = mrdfits(bdata.FILENAME, 'contam')
  cgimage, bcontam, /over, stretch = 2
  oplot, bdata.LAMBDA, bdata.INNERUP, col = 255, thick = 1
  oplot, bdata.LAMBDA, bdata.INNERDN, col = 255, thick = 1
  oplot, bdata.LAMBDA, bdata.INTERUP, col = '00aa00'x, thick = 1
  oplot, bdata.LAMBDA, bdata.INTERDN, col = '00aa00'x, thick = 1
  oplot, bdata.LAMBDA, bdata.OUTERUP, col = 'ff5500'x, thick = 1
  oplot, bdata.LAMBDA, bdata.OUTERDN, col = 'ff5500'x, thick = 1
  
  plot, rdata.LAMBDA, rdata.TRACE, $
        /xsty, yran = [0,60], $
        ytickname = replicate(' ', 60), $
        xtitle = greek('lambda')+' ['+textoIDL('\AA')+']', $
        pos = [0.55,0.15,0.95,0.55], /noer, charsize = 1.1
  cgimage, rdata.SPEC2D, /over, stretch = 2
  cgtext, strmid(rdata.filename, 34, /rev), col = long('0055ff'x), charsize = 1.1, $
          /norm, 0.575, 0.20, charthick = 1.5
  oplot, rdata.LAMBDA, rdata.INNERUP, col = 255, thick = 1
  oplot, rdata.LAMBDA, rdata.INNERDN, col = 255, thick = 1
  oplot, rdata.LAMBDA, rdata.INTERUP, col = '00aa00'x, thick = 1
  oplot, rdata.LAMBDA, rdata.INTERDN, col = '00aa00'x, thick = 1
  oplot, rdata.LAMBDA, rdata.OUTERUP, col = 'ff5500'x, thick = 1
  oplot, rdata.LAMBDA, rdata.OUTERDN, col = 'ff5500'x, thick = 1
  
  plot, rdata.LAMBDA, rdata.TRACE, $
        /xsty, yran = [0,60], $
        xtickname = replicate(' ', 60), ytickname = replicate(' ', 60), $
        pos = [0.55,0.55,0.95,0.95], /noer
  rcontam = mrdfits(rdata.FILENAME, 'contam')
  cgimage, rcontam, /over, stretch = 2
  oplot, rdata.LAMBDA, rdata.INNERUP, col = 255, thick = 1
  oplot, rdata.LAMBDA, rdata.INNERDN, col = 255, thick = 1
  oplot, rdata.LAMBDA, rdata.INTERUP, col = '00aa00'x, thick = 1
  oplot, rdata.LAMBDA, rdata.INTERDN, col = '00aa00'x, thick = 1
  oplot, rdata.LAMBDA, rdata.OUTERUP, col = 'ff5500'x, thick = 1
  oplot, rdata.LAMBDA, rdata.OUTERDN, col = 'ff5500'x, thick = 1

  window, 0
  
end
