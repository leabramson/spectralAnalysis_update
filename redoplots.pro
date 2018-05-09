;; Program to redo most of the plots in the paper, e.g., after
;; refitting spectra

pro redoplots

  plotuvj, 'study5_lgn_Bestfits.list', 'study5_exp_Bestfits.list'
  plotsizemass, 'study5_lgn_Bestfits.list', 'study5_exp_Bestfits.list', $
                /circularize
  plotsfms, 'study5_lgn_Bestfits.list', 'study5_exp_Bestfits.list', $
            nsig = 2
  plotprofiles, 'study5_lgn_Bestfits.list', 'study5_exp_Bestfits.list'
  plotsfhs
  btslope, 'study5_lgn_Structs.list', 'study5_lgn_Bestfits.list', dt = 1

  plot1evoprof_folded, '00900_1_lgn_bestfit.fits','00900_1_lgn_recons.fits', '00900_1_evostuff.eps'
  plot1evoprof_folded, '00660_2_lgn_bestfit.fits', '00660_2_lgn_recons.fits', '00660_2_evostuff.eps'
  plot1evoprof_folded, '00451_2_lgn_bestfit.fits', '00451_2_lgn_recons.fits', '00451_2_evostuff.eps'
  plot1evoprof_folded, '01916_2_lgn_bestfit.fits', '01916_2_lgn_recons.fits', '01916_2_evostuff.eps'

  plotlgnpars, 'study5_lgn_Chains.list'
  plotexppars, 'study5_exp_Chains.list'
  
end
