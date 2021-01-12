# Supplemental program 12.2

# write as function?


#params = list()
#physcon = list()
#atmos = list()
#flux = list()
#leaf = list()
#ground = list()


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
atmos$time = Hainich5Days$Date.Time
atmos$tair_i = Hainich5Days$TA_F + 273.15
atmos$co2air_i = Hainich5Days$CO2
atmos$relhum = Hainich5Days$RH
atmos$wind = Hainich5Days$WS_F
atmos$patm = Hainich5Days$PA_F * 1000
atmos$irsky = Hainich5Days$LW_IN_F
atmos$swr = Hainich5Days$SW_IN_F

atmos$eair = 0.61094*exp((17.625*atmos$tair_i)/(atmos$tair_i+243.04))
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
flux$tleaf_i = atmos$tair_i;

#####

StomatalConductanceLeafPhotosynthesis = function(flux,leaf,params,physcon,atmos){

# -------------------------------------------------------------------------
# Calculate leaf gas exchange coupled with the leaf energy budget for C3
# and C4 plants. Leaf temperature is calculated from the energy balance.


# call these:
source("satvap.R")
source("LeafPhysiologyParams.R")
source("LeafFluxes.R")
source("LeafPhotosynthesis.R")

# --- Waveband indices for visible and near-infrared

params$vis = 1; params$nir = 2

# --- Set leaf physiology variables

# Stomatal conductance: 0 = Medlyn model. 1 = Ball-Berry model. 2 = WUE optimization
# Photosynthetic pathway: 1 = C3. 0 = C4
# Photosynthesis co-limitation: 0 = no. 1 = yes
# leaf$gstyp = 1; leaf$c3psn = 1; leaf$colim = 1;

# Leaf physiological parameters

leaf = LeafPhysiologyParams(params,physcon,leaf);

output_an = c()
output_gs = c()

for (i in 1:240){
  atmos$co2air = atmos$co2air_i[i]
  atmos$tair = atmos$tair_i[i]
  flux$tleaf = flux$tleaf_i[i]
  
  # leafphysiologyparams if vcmaxse, jmaxse and rdse dependent on atmos$tair
  # leaf$rdse = 
  leaf$vcmaxse = 668.39 - 1.07 * atmos$tair
  leaf$jmaxse = 659.7 - 0.75 * atmos$tair
  
  
  LeafPhotosynthesis = function(physcon, atmos, leaf, flux){
    
    #source("sp_12_02.R")  
    source("hybrid_root_ci.R")
    source("satvap.R")
    source("CiFunc.R")
    library("signal")
    library("pracma")
    #source("LeafBoundaryLayer.R") 
    
    # for testing purposes:
    # 1.: get physcon, atmos,leaf from sp_12_02.R (leaf via LeafPhysiologyParams.R)
    #params = list()
    #physcon = list()
    #atmos = list()
    #flux = list()
    #leaf = list()
    
    #helplist = sp_12_02(flux,leaf,params,physcon,atmos)
    
    # 2. set initial leaf temperature:
    #flux$tleaf = atmos$tair;
    # 3. get flux$gbv, flux$gbc, flux$apar from LeafBoundaryLayer.R
    # flux = LeafBoundaryLayer(physcon, atmos, leaf, flux)
    # testing with these values return gs = 0.01, which is plausible
    
    # Calculate leaf photosynthesis using one of two methods:
    # (1) Calculate leaf photosynthesis and stomatal conductance
    # by solving for the value of Ci that satisfies the
    # metabolic, stomatal constraint, and diffusion equations.
    # This is used with the Ball-Berry style stomatal models.
    # (2) Calculate leaf photosynthesis for a specified stomatal
    # conductance. Then calculate Ci from the diffusion equation.
    # This is used with the WUE stomatal optimization.
    
    # ------------------------------------------------------
    # Input
    #   physcon$tfrz        ! Freezing point of water (K)
    #   physcon$rgas        ! Universal gas constant (J/K/mol)
    #   atmos$co2air        ! Atmospheric CO2 (umol/mol)
    #   atmos$eair          ! Vapor pressure of air (Pa)
    #   leaf$gstyp          ! Stomatal conductance: 0 = Medlyn. 1 = Ball-Berry. 2 = WUE optimization
    #   leaf$c3psn          ! Photosynthetic pathway: 1 = C3. 0 = C4 plant
    #   leaf$vcmax25        ! Maximum carboxylation rate at 25C (umol/m2/s)
    #   leaf$jmax25         ! Maximum electron transport rate at 25C (umol/m2/s)
    #   leaf$rd25           ! Leaf respiration rate at 25C (umol CO2/m2/s)
    #   leaf$kc25           ! Michaelis-Menten constant for CO2 at 25C (umol/mol)
    #   leaf$ko25           ! Michaelis-Menten constant for O2 at 25C (mmol/mol)
    #   leaf$cp25           ! CO2 compensation point at 25C (umol/mol)
    #   leaf$kcha           ! Activation energy for Kc (J/mol)
    #   leaf$koha           ! Activation energy for Ko (J/mol)
    #   leaf$cpha           ! Activation energy for Cp (J/mol)
    #   leaf$vcmaxha        ! Activation energy for Vcmax (J/mol)
    #   leaf$jmaxha         ! Activation energy for Jmax (J/mol)
    #   leaf$rdha           ! Activation energy for Rd (J/mol)
    #   leaf$vcmaxhd        ! Deactivation energy for Vcmax (J/mol)
    #   leaf$jmaxhd         ! Deactivation energy for Jmax (J/mol)
    #   leaf$rdhd           ! Deactivation energy for Rd (J/mol)
    #   leaf$vcmaxse        ! Entropy term for Vcmax (J/mol/K)
    #   leaf$jmaxse         ! Entropy term for Jmax (J/mol/K)
    #   leaf$rdse           ! Entropy term for Rd (J/mol/K)
    #   leaf$vcmaxc         ! Vcmax scaling factor for high temperature inhibition (25 C = 1.0)
    #   leaf$jmaxc          ! Jmax scaling factor for high temperature inhibition (25 C = 1.0)
    #   leaf$rdc            ! Rd scaling factor for high temperature inhibition (25 C = 1.0)
    #   leaf$phi_psii       ! Quantum yield of PS II
    #   leaf$theta_j        ! Empirical curvature parameter for electron transport rate
    #   leaf$kp25_c4        ! C4: Initial slope of CO2 response curve at 25C (mol/m2/s)
    #   leaf$g0             ! Ball-Berry minimum leaf conductance (mol H2O/m2/s)
    #   leaf$g1             ! Ball-Berry slope of conductance-photosynthesis relationship
    #   flux$gbv            ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
    #   flux$gbc            ! Leaf boundary layer conductance, CO2 (mol CO2/m2 leaf/s)
    #   flux$apar           ! Leaf absorbed PAR (umol photon/m2 leaf/s)
    #   flux$tleaf          ! Leaf temperature (K)
    #
    # Input or output (depending on method)
    #   flux$gs             ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    #
    # Output
    #   flux$vcmax          ! Maximum carboxylation rate (umol/m2/s)
    #   flux$jmax           ! Maximum electron transport rate (umol/m2/s)
    #   flux$cp             ! CO2 compensation point (umol/mol)
    #   flux$kc             ! Michaelis-Menten constant for CO2 (umol/mol)
    #   flux$ko             ! Michaelis-Menten constant for O2 (mmol/mol)
    #   flux$je             ! Electron transport rate (umol/m2/s)
    #   flux$kp_c4          ! C4: Initial slope of CO2 response curve (mol/m2/s)
    #   flux$ac             ! Leaf Rubisco-limited gross photosynthesis (umol CO2/m2 leaf/s)
    #   flux$aj             ! Leaf RuBP regeneration-limited gross photosynthesis (umol CO2/m2 leaf/s)
    #   flux$ap             ! Leaf product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m2 leaf/s)
    #   flux$ag             ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
    #   flux$an             ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    #   flux$rd             ! Leaf respiration rate (umol CO2/m2 leaf/s)
    #   flux$cs             ! Leaf surface CO2 (umol/mol)
    #   flux$ci             ! Leaf intercellular CO2 (umol/mol)
    #   flux$hs             ! Leaf fractional humidity at surface (dimensionless)
    #   flux$vpd            ! Leaf vapor pressure deficit at surface (Pa)
    # ------------------------------------------------------
    
    # --- Adjust photosynthetic parameters for temperature
    
    # if (leaf$c3psn == 1) 1 = C3
    
    # C3 temperature response
    
    ft = function(tl, ha) {exp(ha/(physcon$rgas*(physcon$tfrz+25)) * (1-(physcon$tfrz+25)/tl));}
    fth = function(tl, hd, se, fc) {fc / (1 + exp((-hd+se*tl)/(physcon$rgas*tl)));}
    
    flux$kc = leaf$kc25 * ft(flux$tleaf, leaf$kcha);
    flux$ko = leaf$ko25 * ft(flux$tleaf, leaf$koha);
    flux$cp = leaf$cp25 * ft(flux$tleaf, leaf$cpha);
    
    t1 = ft(flux$tleaf, leaf$vcmaxha);
    t2 = fth(flux$tleaf, leaf$vcmaxhd, leaf$vcmaxse, leaf$vcmaxc);
    flux$vcmax = leaf$vcmax25 * t1 * t2;
    
    t1 = ft(flux$tleaf, leaf$jmaxha);
    t2 = fth(flux$tleaf, leaf$jmaxhd, leaf$jmaxse, leaf$jmaxc);
    flux$jmax = leaf$jmax25 * t1 * t2;
    
    t1 = ft(flux$tleaf, leaf$rdha);
    t2 = fth(flux$tleaf, leaf$rdhd, leaf$rdse, leaf$rdc);
    flux$rd = leaf$rd25 * t1 * t2;
    
    # --- Electron transport rate for C3 plants
    
    # Solve the polynomial: aquad*Je^2 + bquad*Je + cquad = 0
    # for Je. Correct solution is the smallest of the two roots.
    
    # if (leaf$c3psn == 1) 1 = C3
    
    qabs = 0.5 * leaf$phi_psii * flux$apar;
    aquad = leaf$theta_j;
    bquad = -(qabs + flux$jmax);
    cquad = qabs * flux$jmax;
    pcoeff = c(aquad,bquad,cquad);
    proots = roots(pcoeff); 
    # polyroot(R) = roots(m), but reversed order
    # pcoeff = c(cquad,bquad,aquad);
    # proots = polyroot(pcoeff)
    #is.complex(proots[1])
    #is.complex(proots[2])
    proots[1] = Re(proots[1]);
    proots[2] = Re(proots[2]);
    flux$je = min(proots[1], proots[2]);
    
    # --- Ci calculation
    
    # if (leaf$gstyp <= 1) 1 = Ball-Berry
    
    # Initial estimates for Ci
    
    #if (leaf$c3psn == 1) 1 = C3
    ci0 = 0.7 * atmos$co2air;
    ci1 = ci0 * 0.99;
    
    # Solve for Ci: Use CiFunc to iterate photosynthesis calculations
    # until the change in Ci is < tol. Ci has units umol/mol
    
    tol = 0.1;                 # Accuracy tolerance for Ci (umol/mol)
    # func_name = 'CiFunc';      # The function name
    
    flux_dummy = hybrid_root_ci (physcon, atmos, leaf, flux, ci0, ci1, tol); 
    flux = flux_dummy[[1]] # careful about lists!
    #flux = flux[[1]]
    flux$ci = flux_dummy[[2]];
    
    # --- Relative humidity and vapor pressure at leaf surface
    
    esat = satvap ((flux$tleaf-physcon$tfrz));
    flux$hs = (flux$gbv * atmos$eair + flux$gs * esat) / ((flux$gbv + flux$gs) * esat);
    flux$vpd = max(esat - flux$hs*esat, 0.1);
    
    # --- Make sure iterative solution is correct
    
    if (flux$gs < 0) {
      stop ('LeafPhotosynthesis: negative stomatal conductance')
    }
    
    # Compare with Ball-Berry model. The solution blows up with low eair. In input
    # data, eair should be > 0.05*esat to ensure that hs does not go to zero.
    
    #if (leaf.gstyp == 1){
    gs_err = leaf$g1 * max(flux$an, 0) * flux$hs / flux$cs + leaf$g0;
    if (abs(flux$gs-gs_err)*1e06 > 1e-04) {
      fprintf('gs = #15.4f\n', flux$gs)
      fprintf('gs_err = #15.4f\n', gs_err)
      stop ('LeafPhotosynthesis: failed Ball-Berry error check')
    }
    
    # Compare with Medlyn model. The solutions blows up with vpd = 0. The
    # quadratic calcuation of gsw in CiFunc constrains vpd > 50 Pa, so this
    # comparison is only valid for those conditions.
    
    # if (leaf$gstyp == 0)
    #  if ((esat - atmos$eair) > 50)
    #    gs_err = 1.6 * (1 + leaf$g1 / sqrt(flux$vpd*0.001)) * max(flux$an, 0) / flux$cs + leaf$g0;
    #if (abs(flux$gs-gs_err)*1e06 > 1e-04)
    #  fprintf('gs = #15.4f\n', flux$gs)
    #fprintf('gs_err = #15.4f\n', gs_err)
    #error ('LeafPhotosynthesis: failed Medlyn error check')
    #end
    #end
    #end
    
    # Compare with diffusion equation: An = (ca - ci) * gleaf
    
    an_err = (atmos$co2air - flux$ci) / (1 / flux$gbc + 1.6 / flux$gs);
    if (flux$an > 0 & abs(flux$an-an_err) > 0.01){
      fprintf('An = #15.4f\n', flux$an)
      fprintf('An_err = #15.4f\n', an_err)
      stop ('LeafPhotosynthesis: failed diffusion error check')
    }
    
    return(flux)
    
}
  flux = LeafPhotosynthesis(physcon, atmos, leaf, flux)
  output_an[i] = flux$an
  output_gs[i] = flux$gs
  #output_time[i] = atmos$tl
}

#atmos$time == time

output_angs = data.frame(
  Time = atmos$time,
  an = output_an,
  gs = output_gs
  )

write.csv(output_angs,file = "testoutputs1")

plot(output_an ~ output_gs)
plot(output_an ~ atmos$tair_i)
plot(output_an ~ atmos$eair)
plot(output_an ~ Hainich5Days$NIGHT)
plot(output_an ~ Hainich5Days$TIMESTAMP_END, type = "l")

plot(output_gs ~ output_gs)
plot(output_gs ~ atmos$tair_i)
plot(output_gs ~ atmos$eair)
plot(output_gs ~ Hainich5Days$NIGHT)
plot(output_gs ~ Hainich5Days$TIMESTAMP_END, type = "l")



# --- Flux calculations for 20 leaves with dleaf = 1 - 20 cm

# for (p in 1:20) {
  
  #leaf$dleaf = p / 100;
  
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
#}

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