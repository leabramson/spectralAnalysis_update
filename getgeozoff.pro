function getGeoZoff, fitsname, region

  struct = mrdfits(fitsname, 1, /silent)
  
  ctr = struct.LSF_INNER_PARS[1]
  dl  = struct.LAMBDA[1] - struct.LAMBDA[0]
  ml  = median(struct.LAMBDA)
  
  case REGION of
     'INTUP': c2 = struct.LSF_INTER_UP_PARS[1]
     'INTDN': c2 = struct.LSF_INTER_DN_PARS[1]
     'OUTUP': c2 = struct.LSF_OUTER_UP_PARS[1]
     'OUTDN': c2 = struct.LSF_OUTER_DN_PARS[1]
  endcase

  loff = (c2 - ctr) * dl
  zoff = dl / ml
  
  RETURN, zoff
end

