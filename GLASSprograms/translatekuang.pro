pro translatekuang, fitsIn, fitsOut

  d1 = mrdfits(fitsIn, 1)
  tags = tag_names(d1)

  ;; Compensate for known Kuang catalog key issues.
  
  if total(tags eq 'Z_CHI_A') eq 1 then $
     chitag = 'Z_CHI_A' $
  else $
     chitag = 'CHI_A'
  chipos = (where(tags eq chitag))[0]

  if total(tags eq 'Z_L68') eq 1 then begin
     e99ltag = 'Z_L99'
     e99utag = 'Z_U99'
     e95ltag = 'Z_L95'
     e95utag = 'Z_U95'
     e68ltag = 'Z_L68'
     e68utag = 'Z_U68'
  endif else begin
     e99ltag = 'L99'
     e99utag = 'U99'
     e95ltag = 'L95'
     e95utag = 'U95'
     e68ltag = 'L68'
     e68utag = 'U68'
  endelse
  e68lpos = (where(tags eq e68ltag))[0]
  e68upos = (where(tags eq e68utag))[0]
  e95lpos = (where(tags eq e95ltag))[0]
  e95upos = (where(tags eq e95utag))[0]
  e99lpos = (where(tags eq e99ltag))[0]
  e99upos = (where(tags eq e99utag))[0]
  
  
  if total(tags eq 'Z_NFILT') eq 1 then $
     nfilttag = 'Z_NFILT' $
  else $
     nfilttag = 'NFILT'
  nfiltpos = (where(tags eq nfilttag))[0]

  if total(tags eq 'F435W_FLUX_RADIUS') eq 1 then $
     radinc = 1 $
  else $
     radinc = 0

  if radinc then begin

     tags_to_keep = ['Z_A',$
                     'Z_M1', $
                     chitag, $
                     e68ltag, $
                     e68utag, $
                     e95ltag, $
                     e95utag, $
                     e99ltag, $
                     e99utag, $
                     nfilttag, $
                     'ALPHA_J2000', $
                     'DELTA_J2000', $
                     'NUMBER', $
                     'F435W_MAG_AUTO', $
                     'F606W_MAG_AUTO', $
                     'F814W_MAG_AUTO', $
                     'F105W_MAG_AUTO', $
                     'F125W_MAG_AUTO', $
                     'F140W_MAG_AUTO', $
                     'F160W_MAG_AUTO', $
                     'F435W_FLUX_RADIUS', $
                     'F606W_FLUX_RADIUS', $
                     'F814W_FLUX_RADIUS', $
                     'F105W_FLUX_RADIUS', $
                     'F125W_FLUX_RADIUS', $
                     'F140W_FLUX_RADIUS', $
                     'F160W_FLUX_RADIUS']

     d2 = {KUANG_PZ_A        : 0., $
           KUANG_PZ_M1       : 0., $
           KUANG_PZ_CHI_A    : 0., $
           KUANG_PZ_L68      : 0., $
           KUANG_PZ_U68      : 0., $
           KUANG_PZ_L95      : 0., $
           KUANG_PZ_U95      : 0., $
           KUANG_PZ_L99      : 0., $
           KUANG_PZ_U99      : 0., $
           KUANG_PZ_NFILT    : 0,  $
           KF435W_MAG_AUTO   : 0., $
           KF606W_MAG_AUTO   : 0., $
           KF814W_MAG_AUTO   : 0., $
           KF105W_MAG_AUTO   : 0., $
           KF125W_MAG_AUTO   : 0., $
           KF140W_MAG_AUTO   : 0., $
           KF160W_MAG_AUTO   : 0., $
           KF435W_FLUX_RADIUS: 0., $
           KF606W_FLUX_RADIUS: 0., $
           KF814W_FLUX_RADIUS: 0., $
           KF105W_FLUX_RADIUS: 0., $
           KF125W_FLUX_RADIUS: 0., $
           KF140W_FLUX_RADIUS: 0., $
           KF160W_FLUX_RADIUS: 0., $
           RA                : 0.d, $
           DEC               : 0.d, $
           KUANG_ID          : 0L}

     d2.KF435W_FLUX_RADIUS = d1.F435W_FLUX_RADIUS
     d2.KF606W_FLUX_RADIUS = d1.F606W_FLUX_RADIUS
     d2.KF814W_FLUX_RADIUS = d1.F814W_FLUX_RADIUS
     d2.KF105W_FLUX_RADIUS = d1.F105W_FLUX_RADIUS
     d2.KF125W_FLUX_RADIUS = d1.F125W_FLUX_RADIUS
     d2.KF140W_FLUX_RADIUS = d1.F140W_FLUX_RADIUS
     d2.KF160W_FLUX_RADIUS = d1.F160W_FLUX_RADIUS

  endif else begin

     tags_to_keep = ['Z_A',$
                     'Z_M1', $
                     chitag, $
                     e68ltag, $
                     e68utag, $
                     e95ltag, $
                     e95utag, $
                     e99ltag, $
                     e99utag, $
                     nfilttag, $
                     'ALPHA_J2000', $
                     'DELTA_J2000', $
                     'NUMBER', $
                     'F435W_MAG_AUTO', $
                     'F606W_MAG_AUTO', $
                     'F814W_MAG_AUTO', $
                     'F105W_MAG_AUTO', $
                     'F125W_MAG_AUTO', $
                     'F140W_MAG_AUTO', $
                     'F160W_MAG_AUTO']

     d2 = {KUANG_PZ_A        : 0., $
           KUANG_PZ_M1       : 0., $
           KUANG_PZ_CHI_A    : 0., $
           KUANG_PZ_L68      : 0., $
           KUANG_PZ_U68      : 0., $
           KUANG_PZ_L95      : 0., $
           KUANG_PZ_U95      : 0., $
           KUANG_PZ_L99      : 0., $
           KUANG_PZ_U99      : 0., $
           KUANG_PZ_NFILT    : 0,  $
           KF435W_MAG_AUTO   : 0., $
           KF606W_MAG_AUTO   : 0., $
           KF814W_MAG_AUTO   : 0., $
           KF105W_MAG_AUTO   : 0., $
           KF125W_MAG_AUTO   : 0., $
           KF140W_MAG_AUTO   : 0., $
           KF160W_MAG_AUTO   : 0., $
           RA                : 0.d, $
           DEC               : 0.d, $
           KUANG_ID          : 0L}

  endelse

  d2 = replicate(d2, n_elements(d1))

  d2.KUANG_PZ_A         = d1.Z_A     
  d2.KUANG_PZ_M1        = d1.Z_M1    
  d2.KUANG_PZ_CHI_A     = d1.(chipos) 
  d2.KUANG_PZ_L68       = d1.(e68lpos)    
  d2.KUANG_PZ_U68       = d1.(e68upos)    
  d2.KUANG_PZ_L95       = d1.(e95lpos)       
  d2.KUANG_PZ_U95       = d1.(e95upos)    
  d2.KUANG_PZ_L99       = d1.(e99lpos)    
  d2.KUANG_PZ_U99       = d1.(e99upos)    
  d2.KUANG_PZ_NFILT     = d1.(nfiltpos)    
  d2.KF435W_MAG_AUTO    = d1.F435W_MAG_AUTO   
  d2.KF606W_MAG_AUTO    = d1.F606W_MAG_AUTO   
  d2.KF814W_MAG_AUTO    = d1.F814W_MAG_AUTO   
  d2.KF105W_MAG_AUTO    = d1.F105W_MAG_AUTO   
  d2.KF125W_MAG_AUTO    = d1.F125W_MAG_AUTO   
  d2.KF140W_MAG_AUTO    = d1.F140W_MAG_AUTO   
  d2.KF160W_MAG_AUTO    = d1.F160W_MAG_AUTO   
  d2.RA                 = d1.ALPHA_J2000
  d2.DEC                = d1.DELTA_J2000
  d2.KUANG_ID           = d1.NUMBER
  
  mwrfits, d2, fitsOut, /create  
  
end
;translateKuang, 'KUANG_PZS/CTRS/hst_macs0717_hffclash_psfmatch_60mas_photoz.fits', 'M0717_Kuang_photoz.fits'
