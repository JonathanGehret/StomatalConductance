# function [flux, tleaf_dif] = TleafFunc (physcon, atmos, leaf, flux, tleaf_val)

TleafFunc = function(physcon, atmos, leaf, flux, tleaf_val){

# Calculate leaf temperature and fluxes for an input leaf temperature
# (tleaf_val) and compare the new temperature to the prior
# temperature. This function returns a value tleaf_dif = 0 when leaf
# temperature does not change between iterations.
  
  
  
# call necessary functions:
source("LeafBoundaryLayer.R")
source("LeafPhotosynthesis.R")
source("LeafTemperature.R")

if (tleaf_val < 0){
  stop ('TleafFunc error')
}

# --- Current value for leaf temperature

flux$tleaf = tleaf_val;

# --- Leaf boundary layer conductances

flux = LeafBoundaryLayer (physcon, atmos, leaf, flux);
# Output
#   flux$gbh          ! Leaf boundary layer conductance, heat (mol/m2 leaf/s)
#   flux$gbv          ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
#   flux$gbc          ! Leaf boundary layer conductance, CO2 (mol CO2/m2 leaf/s)
# -------------------------------------------------------------------------

#test values

flux$gbh = 0.638
# gbw typically between (~1–4 mol m–2 s–1) 
flux$gbv = 0.702
flux$gbc = 0.517

# --- Leaf photosynthesis and stomatal conductance

flux = LeafPhotosynthesis (physcon, atmos, leaf, flux);

# --- Leaf temperature and energy fluxes

flux = LeafTemperature (physcon, atmos, leaf, flux);

# ------------------------------------------------------
# Input/ouput
#   flux$tleaf       ! Leaf temperature (K)
#
# Output
#   flux$rnet        ! Leaf net radiation (W/m2 leaf)
#   flux$lwrad       ! Longwave radiation emitted from leaf (W/m2 leaf)
#   flux$shflx       ! Leaf sensible heat flux (W/m2 leaf)
#   flux$lhflx       ! Leaf latent heat flux (W/m2 leaf)
#   flux$etflx       ! Leaf transpiration flux (mol H2O/m2 leaf/s)
# ------------------------------------------------------

# therefore, for leaf.temperature for now:
# flux$tleaf = 
# flux$rnet =      
# flux$lwrad =    
# flux$shflx =    
# flux$lhflx =    
# flux$etflx =   

# --- Compare with prior value for leaf temperature

tleaf_dif = flux$tleaf - tleaf_val;

flux_tleaf_dif = list(flux,tleaf_dif)

return(flux_tleaf_dif)
}
