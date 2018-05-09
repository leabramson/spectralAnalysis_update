pro plotmasscomptime

;00451_2_lgn_structEvo.fits  00660_2_lgn_structEvo.fits	00900_1_lgn_structEvo.fits  01916_2_lgn_structEvo.fits
  data = mrdfits('00451_2_lgn_structEvo.fits',1)

  time = data[0].TIME
  odex = value_locate(time, data[0].TOBS) - 1

  oomass = data[-1].OMASS[odex,0]
  eoomass = data[-1].OMASS[odex,1]
  otmass = data[-1].TMASS[odex,0]
  eotmass = data[-1].TMASS[odex,1]

  oosfr = data[-1].OSFR[odex,0]
  eoosfr = data[-1].OSFR[odex,1]
  otsfr = data[-1].TSFR[odex,0]
  eotsfr = data[-1].TSFR[odex,1]
  
  xxx = [time, reverse(time)]
  yyy1 = [data[-1].OMASS[*,0] - data[-1].OMASS[*,1], reverse(data[-1].OMASS[*,0] + data[-1].OMASS[*,1])]
  yyy2 = [data[-1].TMASS[*,0] - data[-1].TMASS[*,1], reverse(data[-1].TMASS[*,0] + data[-1].TMASS[*,1])]

  yyy3 = [data[-1].OSFR[*,0] - data[-1].OSFR[*,1], reverse(data[-1].OSFR[*,0] + data[-1].OSFR[*,1])]
  yyy4 = [data[-1].TSFR[*,0] - data[-1].TSFR[*,1], reverse(data[-1].TSFR[*,0] + data[-1].TSFR[*,1])]

  plotsym, 0, /fill
  
  plot, time, data[-1].OMASS[*,0], yran = [9,12], $
        xtitle = '!18t!X [Gyr]', ytitle = 'log !18M!X!D*!N or log !18SFR!X + 8', $
        xran = [1,data[0].TOBS + 2], /nodat
  polyfill, (xxx > !X.CRANGE[0]) < !X.CRANGE[1], yyy1 > !Y.CRANGE[0]
  polyfill, (xxx > !X.CRANGE[0]) < !X.CRANGE[1], yyy2 > !Y.CRANGE[0], $
            /line_fill, col = 255, orien = 45, thick = 1, spacing = 0.05
  polyfill, (xxx > !X.CRANGE[0]) < !X.CRANGE[1], (yyy3 + 8) > !Y.CRANGE[0], col = '777777'x
  polyfill, (xxx > !X.CRANGE[0]) < !X.CRANGE[1], (yyy4 + 8) > !Y.CRANGE[0], $
            /line_fill, col = 'ff5500'x, orien = 45, thick = 1, spacing = 0.05
;  oploterror, [data[0].TOBS], [oomass], [eoomass], psym = 8, /nohat, errcol = 0
;  oploterror, [data[0].TOBS], [otmass], [eotmass], psym = 8, /nohat, errcol = 200

  stop

end
