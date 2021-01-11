gesamtfuncion = function(flux,leaf){

flux
atmos

source("leafboundarylayer.R")

for i in 5days:{
  
  LeafBoundaryLayer(flux,leaf) # returns flux_ci_x
  
  flux = flux_ci_x[1,]
  SoilMoisture(flux)
  
  LeafPhotosynthesis()
}  
  

return(flux)