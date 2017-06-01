;; 10^-18 erg/s/cm^2
function gethasfr, haf, z, $
                   AV = av, $
                   RV = rv, $
                   HBF = hbf, $
                   DUSTCOR = dustcor

  if NOT keyword_set(RV) then rv = 3.1
  
  if keyword_set(DUSTCOR) then begin

     ;; USE THE CARDELLI+89 LAW TO DEREDDEN SHIT
         
     hal = 1./0.6563 ;; ha line position in 1/micron
     hbl = 1./0.4863

     ya = hal - 1.82
     yb = hbl - 1.82

     a = [1, 0.17699, -0.50447, -0.02427, 0.72085, $
          0.01979, -0.77530, 0.32999]
     b = [0, 1.41338, 2.28305, 1.07233, -5.38434, $
          -0.62251, 5.30260, -2.09002]

     ka = rv * poly(ya, a) + poly(ya, b)
     kb = rv * poly(yb, a) + poly(yb, b)

     ebv = 2.5 * alog10(2.86 * hbf / haf) / (ka - kb)

     haf *= 10.^(ebv * ka / 2.5)
     
  endif

  sfr = haf * 1d-18 * dluminosity(z, /cm)^2 * 4 * !pi / 1.26d41
  
  RETURN, sfr
end
