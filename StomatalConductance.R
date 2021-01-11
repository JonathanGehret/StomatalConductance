# Supplemental program 12.2

# write as function?


params = list()
physcon = list()
inputatmos = list()
flux = flux()
leaf = list()
ground = list()


# -------------------------------------------------------------------------

#   atmos$patm        ! Atmospheric pressure (Pa)
#   atmos$rhomol      ! Molar density mol/m3)
#   atmos$wind        ! Wind speed (m/s)
#   atmos$tair        ! Air temperature (K)
#   atmos$o2air         ! Atmospheric O2 (mmol/mol)
#   atmos$co2air        ! Atmospheric CO2 (umol/mol)
#   atmos$eair          ! Vapor pressure of air (Pa)
#   atmos$irsky         ! Atmospheric longwave radiation (W/m2)
#   atmos$swr           ! SW radiation (W/m2)
#   atmos$qair          ! specific humidity (kg/kg)

Hainich5Days = read.csv("5_days.csv",header = T, dec = ",", sep = ";")
# inputatmos$eair = Hainich5Days$ea * 1000
atmos$tair = Hainich5Days$TA_F + 273.15
atmos$co2air = Hainich5Days$CO2
atmos$relhum = Hainich5Days$RH
atmos$wind = Hainich5Days$WS_F
atmos$patm = Hainich5Days$PA_F * 1000
atmos$irsky = Hainich5Days$LW_IN_F
atmos$swr = Hainich5Days$SW_IN_F
atmos$eair = 0.61094*exp((17.625*atmos$tair)/(atmos$tair+243.04))
atmos$qair = physcon$mmh2o / physcon$mmdry * atmos$eair / (atmos$patm - (1 - physcon$mmh2o/physcon$mmdry) * atmos$eair);
atmos$o2air = 0.209 * 1000;           # Atmospheric  and O2 (mmol/mol)
atmos$rhomol = atmos$patm / (physcon$rgas * atmos$tair); # Molar density (mol/m3)
atmos$rhoair = atmos$rhomol * physcon$mmdry * (1 - (1 - physcon$mmh2o/physcon$mmdry) * atmos$eair / atmos$patm); # Air density (kg/m3)
atmos$mmair = atmos$rhoair / atmos$rhomol;     # Molecular mass of air (kg/mol)
atmos$cpair = physcon$cpd * (1 + (physcon$cpw/physcon$cpd - 1) * atmos$qair) * atmos$mmair; # Specific heat of air at constant pressure (J/mol/K)
atmos$swsky[params$vis] = 0.5 * atmos$swr;   # short wave sky
atmos$swsky[params$nir] = 0.5 * atmos$swr;   # short wave sky

#   flux$apar           ! Leaf absorbed PAR (umol photon/m2 leaf/s)

flux$apar = Hainich5Days$PPFD_IN - Hainich5Days$PPFD_OUT

##### --- Physical constants ####

physcon$grav = 9.80665;               # Gravitational acceleration (m/s2)
physcon$tfrz = 273.15;                # Freezing point of water (K)
physcon$sigma = 5.67e-08;             # Stefan-Boltzmann constant (W/m2/K4)
physcon$mmdry = 28.97 / 1000;         # Molecular mass of dry air (kg/mol)
physcon$mmh2o = 18.02 / 1000;         # Molecular mass of water (kg/mol)
physcon$cpd = 1005;                   # Specific heat of dry air at constant pressure (J/kg/K)
physcon$cpw = 1846;                   # Specific heat of water vapor at constant pressure (J/kg/K)
physcon$rgas = 8.31446;               # Universal gas constant (J/K/mol)
physcon$visc0 = 13.3e-06;             # Kinematic viscosity at 0C and 1013.25 hPa (m2/s)
physcon$Dh0 = 18.9e-06;               # Molecular diffusivity (heat) at 0C and 1013.25 hPa (m2/s)
physcon$Dv0 = 21.8e-06;               # Molecular diffusivity (H2O) at 0C and 1013.25 hPa (m2/s)
physcon$Dc0 = 13.8e-06;               # Molecular diffusivity (CO2) at 0C and 1013.25 hPa (m2/s)

# Test values boundary layer conductance

flux$gbh = 0.638
# gbw typically between (~1–4 mol m–2 s–1) 
flux$gbv = 0.702
flux$gbc = 0.517


# Test value leaf temperature
flux$tleaf = atmos$tair;

#####

StomatalConductance = function(flux,leaf,params,physcon,atmos){

# -------------------------------------------------------------------------
# Calculate leaf gas exchange coupled with the leaf energy budget for C3
# and C4 plants. Leaf temperature is calculated from the energy balance.


# call these:
source("satvap.R")
source("LeafPhysiologyParams.R")
source("LeafFluxes.R")

# --- Waveband indices for visible and near-infrared

params$vis = 1; params$nir = 2

# --- Set leaf physiology variables

# Stomatal conductance: 0 = Medlyn model. 1 = Ball-Berry model. 2 = WUE optimization
# Photosynthetic pathway: 1 = C3. 0 = C4
# Photosynthesis co-limitation: 0 = no. 1 = yes
# leaf$gstyp = 1; leaf$c3psn = 1; leaf$colim = 1;

# Leaf physiological parameters

leaf = LeafPhysiologyParams(params,physcon,leaf);

testlist = c()

loop_i = list()
for (i in 1:240){
  
  #params = list()
  #physcon = list()
  loop_i$co2air_i = atmos$co2air[i]
  loop_i$tair_i = atmos$tair[i]
  testfunction(tair_i)
  loop_i$tleaf = flux$tleaf[i]
  #flux_i = flux[i]
  #leaf = list()
  #ground = list()
  # leafphysiologyparams if vcmaxse, jmaxse and rdse dependent on atmos$tair
  flux = LeafPhotosynthesis(physcon, atmos, leaf, flux)
  testlist[i] = atmos$co2air_i
  

  return(print(testlist))
}


# --- Flux calculations for 20 leaves with dleaf = 1 - 20 cm

for (p in 1:20) {
  
  leaf$dleaf = p / 100;
  
  # --- Initial leaf temperature
  
  flux$tleaf = atmos$tair;
  
  # --- Leaf temperature, energy fluxes, photosynthesis, and stomatal conductance
  
  flux = LeafFluxes (physcon, atmos, leaf, flux);
  
  # --- Save data for output
  # either create vectores x1...x10, add a value each itereation to each one and combine in the end
  # using x1 = c()....x10 = c()
  # and x1[p] = leaf$dleaf... x10[p] = flux$gs;
  # or add all of them to a data.frame from the beginning, like here
  
  LeafFluxes_output = data.frame()
  
  x1 = leaf$dleaf * 100;             # m -> cm
  x2 = flux$apar;
  x3 = flux$tleaf - physcon$tfrz;    # K -> oC
  x4 = flux$qa;
  x5 = flux$lhflx;
  x6 = flux$etflx * 1000;            # mol H2O/m2/s -> mmol H2O/m2/s
  x7 = flux$an;
  x8 = flux$an / flux$etflx * 0.001; # mmol CO2 / mol H2O
  x9 = flux$gbh;
  x10 = flux$gs;
  
  LeafFluxes_output[p] = c(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)
}

# --- Plot data

plot(LeafFluxes_output[x1],LeafFluxes_output[x3], xlab = 'Leaf dimension (cm)', ylab = 'Leaf temperature (°C)')

# --- Write data to output file

# A = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10];

write.table(LeafFluxes_output,"data.txt")
# filename = 'data.txt';
# fileID = fopen(filename,'w');
#fprintf(fileID,'#10.3f #10.3f #10.3f #10.3f #10.3f #10.3f #10.3f #10.3f #10.5f #10.5f\n', A);
# fclose(fileID);


# return xxx

 }

}