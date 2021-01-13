
d <- read.table("testoutputs1", header = T, sep=",", dec = ".")
Hainich5Days = read.csv("5_days.csv",header = T, dec = ",", sep = ";")


plot(d$an ~ c(1:240), xlab = "Timesteps", ylab = "Leaf net photosynthesis (umol CO2/m2 leaf/s)")
plot(d$gs ~ c(1:240), xlab = "Timesteps", ylab = "Leaf stomatal conductance (mol H2O/m2 leaf/s)")
plot(d$an ~ Hainich5Days$TA_F, xlab = "Temperature (°C)", ylab = "Leaf net photosynthesis (umol CO2/m2 leaf/s)")
plot(d$gs ~ Hainich5Days$TA_F, xlab = "Temperature (°C)", ylab = "Leaf stomatal conductance (mol H2O/m2 leaf/s)")
plot(d$an ~ Hainich5Days$CO2, xlab = "Ambient CO2 concentration (µmol mol-1)", ylab = "Leaf net photosynthesis (umol CO2/m2 leaf/s)" )
plot(d$gs ~ Hainich5Days$CO2, xlab = "Ambient CO2 concentration (µmol mol-1)", ylab = "Leaf stomatal conductance (mol H2O/m2 leaf/s)" )

plot(Hainich5Days$TA_F ~ c(1:240), xlab ="Timesteps", ylab = "Temperature (°C)")
