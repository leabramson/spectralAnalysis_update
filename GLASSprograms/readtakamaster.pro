pro readtakamaster, infile, output

;   1 NUMBER                 Running object number                                     
;   2 MAG_AUTO               Kron-like elliptical aperture magnitude                    [mag]
;   3 MAGERR_AUTO            RMS error for AUTO magnitude                               [mag]
;   4 KRON_RADIUS            Kron apertures in units of A or B                         
;   5 BACKGROUND             Background at centroid position                            [count]
;   6 X_WORLD                Barycenter position along world x axis                     [deg]
;   7 Y_WORLD                Barycenter position along world y axis                     [deg]
;   8 X_IMAGE                Object position along x                                    [pixel]
;   9 Y_IMAGE                Object position along y                                    [pixel]
;  10 THETA_IMAGE            Position angle (CCW/x)                                     [deg]
;  11 ELONGATION             A_IMAGE/B_IMAGE                                           
;  12 A_IMAGE                Profile RMS along major axis                               [pixel]
;  13 CLASS_STAR             S/G classifier output                                     
;  14 FLUX_RADIUS            Fraction-of-light radii                                    [pixel]
;  16 FLAGS                  Extraction flags                                          
;  17 ISOAREA_IMAGE          Isophotal area above Analysis threshold                    [pixel**2]
;  18 FLUX_AUTO              Flux within a Kron-like elliptical aperture                [count]
;  19 FLUXERR_AUTO           RMS error for AUTO flux                                    [count]
;  20 FLUX_APER              Flux vector within fixed circular aperture(s)              [count]
;  23 FLUXERR_APER           RMS error vector for aperture flux(es)                     [count]
;  26 FLUX_ISO               Isophotal flux                                             [count]
;  27 FLUXERR_ISO            RMS error for isophotal flux                               [count]

  readcol, infile, $
           NUMBER       , $
           MAG_AUTO     , $     
           MAGERR_AUTO  , $
           KRON_RADIUS  , $
           BACKGROUND   , $
           X_WORLD      , $
           Y_WORLD      , $
           X_IMAGE      , $
           Y_IMAGE      , $
           THETA_IMAGE  , $
           ELONGATION   , $
           A_IMAGE      , $
           CLASS_STAR   , $
           FLUX_RADIUS  , $
           FLAGS        , $
           ISOAREA_IMAGE, $
           FLUX_AUTO    , $
           FLUXERR_AUTO , $
           FLUX_APER1   , $
           FLUX_APER2   , $
           FLUX_APER3   , $
           FLUXERR_APER1, $                      
           FLUXERR_APER2, $
           FLUXERR_APER3, $
           FLUX_ISO     , $
           FLUXERR_ISO  , $
           f = $
           'L,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'D,'+$
           'D,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,'+$
           'F,',$
           /silent
  nsrcs = n_elements(number)
  
  savedata = {ID_TAKA      : 0L, $
              MAG_AUTO     : 0., $
              MAGERR_AUTO  : 0., $
              KRON_RADIUS  : 0., $
              BACKGROUND   : 0., $
              RA           : 0.d, $
              DEC          : 0.d, $
              THETA_IMAGE  : 0., $
              ELONGATION   : 0., $
              A_IMAGE      : 0., $
              CLASS_STAR   : 0., $
              FLUX_RADIUS  : 0., $
              FLAGS        : 0., $
              ISOAREA_IMAGE: 0., $
              FLUX_AUTO    : 0., $
              FLUXERR_AUTO : 0., $
              FLUX_APER1   : 0., $
              FLUX_APER2   : 0., $
              FLUX_APER3   : 0., $
              FLUXERR_APER1: 0., $
              FLUXERR_APER2: 0., $
              FLUXERR_APER3: 0., $
              FLUX_ISO     : 0., $
              FLUXERR_ISO  : 0.}

  savedata = replicate(savedata, nsrcs)

  for ii = 0, nsrcs - 1 do begin
     savedata[ii].ID_TAKA       = NUMBER[ii]        
     savedata[ii].MAG_AUTO      = MAG_AUTO[ii]      
     savedata[ii].MAGERR_AUTO   = MAGERR_AUTO[ii]   
     savedata[ii].KRON_RADIUS   = KRON_RADIUS[ii]   
     savedata[ii].BACKGROUND    = BACKGROUND[ii]    
     savedata[ii].RA            = X_WORLD[ii]       
     savedata[ii].DEC           = Y_WORLD[ii]       
     savedata[ii].THETA_IMAGE   = THETA_IMAGE[ii]   
     savedata[ii].ELONGATION    = ELONGATION[ii]    
     savedata[ii].A_IMAGE       = A_IMAGE[ii]       
     savedata[ii].CLASS_STAR    = CLASS_STAR[ii]    
     savedata[ii].FLUX_RADIUS   = FLUX_RADIUS[ii]   
     savedata[ii].FLAGS         = FLAGS[ii]         
     savedata[ii].ISOAREA_IMAGE = ISOAREA_IMAGE[ii] 
     savedata[ii].FLUX_AUTO     = FLUX_AUTO[ii]     
     savedata[ii].FLUXERR_AUTO  = FLUXERR_AUTO[ii]  
     savedata[ii].FLUX_APER1    = FLUX_APER1[ii]    
     savedata[ii].FLUX_APER2    = FLUX_APER2[ii]    
     savedata[ii].FLUX_APER3    = FLUX_APER3[ii]    
     savedata[ii].FLUXERR_APER1 = FLUXERR_APER1[ii]
     savedata[ii].FLUXERR_APER2 = FLUXERR_APER2[ii]
     savedata[ii].FLUXERR_APER3 = FLUXERR_APER3[ii]
     savedata[ii].FLUX_ISO      = FLUX_ISO[ii]
     savedata[ii].FLUXERR_ISO   = FLUXERR_ISO[ii]
  endfor

  mwrfits, savedata, output, /create

end
