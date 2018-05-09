function getoneID, inRA, inDEC, mastercat, $
                   matchlen = matchlen

  if NOT keyword_Set(matchlen) then matchlen = 0.5

  data = mrdfits(mastercat, 1)
  ra = data.X_WORLD
  dec = data.Y_WORLD

  spherematch, ra, dec, inRA, inDEC, $
               double(matchlen) / 3600, $
               minds, tinds, len

  print, ''
  print, 'GLASS ID : ', data[minds].NUMBER
  print, ''
  
  RETURN, data[minds].NUMBER
end
