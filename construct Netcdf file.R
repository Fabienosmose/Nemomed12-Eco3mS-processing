# Script pour créer un nouveau fichier Netcdf modèle biogéochimique Nemomed12-Eco3m avec moyenne sur 8 ans
# Warning ! les variables planctons sont en mmolC donc il faut diviser par la surface des cellules pour avoir des mmolC/m²
# Author : Fabien Moullec
# Date: 09/09/2017
##

rm(list = ls())

# Load libraries
library(ncdf4)
library(plyr)
library(abind)


# To load the surface cell
ncin <- nc_open("C:/Users/Fabien/Documents/Scripts R/Nemomed12-Eco3m/Data/Nemo outputs by steps/1/20060101-20060115-OsmoseInputs.nc")
surface = ncvar_get(ncin, "dxdy_t")

# Récupérer les variables dans chaque Netcdf
# Set the path

j <- 24
wd <- paste0('C:/Users/Fabien/Documents/Scripts R/Nemomed12-Eco3m/Data/Nemo outputs by steps/',j)
setwd(dir = wd)

files = list.files(pattern='nc',full.names=TRUE)

# Loop over files

# Picophyto
tab.1 <- array(NA,c(567,264,7))
for(i in seq_along(files)){
  cat("i",i, "\n")
  ncin = nc_open(files[i])
  v1 = ncvar_get(ncin,'picophyto')
  var.1 <- v1/surface 
  var.1[is.na(var.1)]<- NA
  tab.1[,,i]=var.1
}

assign(paste0("picophy.",j), (aaply(tab.1, c(1,2),mean, na.rm=T)))

# Nanophyto 
tab.2 <- array(NA,c(567,264,7))
for(i in seq_along(files)) {
  cat("i",i, "\n")
  ncin = nc_open(files[i])
  v2 = ncvar_get(ncin,'nanophyto')
  var.2 <- v2/surface 
  var.2[is.na(var.2)]<- NA
  tab.2[,,i]=var.2
}
assign(paste0("nanophy.",j), (aaply(tab.2, c(1,2),mean, na.rm=T)))
  
# Microphyto
tab.3 <- array(NA,c(567,264,7))
for(i in seq_along(files)) {
  ncin = nc_open(files[i])
  v3 = ncvar_get(ncin,'microphyto')
  var.3 <- v3/surface
  var.3[is.na(var.3)]<- NA
  tab.3[,,i]=var.3
}
assign(paste0("microphyto.",j),(aaply(tab.3, c(1,2),mean, na.rm=T)))

# Nanozoo
tab.4 <- array(NA,c(567,264,7))
for(i in seq_along(files)) {
  ncin = nc_open(files[i])
  v4 = ncvar_get(ncin,'nanozoo')
  var.4 <- v4/surface
  var.4[is.na(var.4)]<- NA
  tab.4[,,i]=var.4
}
assign(paste0("nanozo.",j),(aaply(tab.4, c(1,2),mean, na.rm=T)))

# Microzoo
tab.5 <- array(NA,c(567,264,7))
for(i in seq_along(files)) {
  ncin = nc_open(files[i])
  v5 = ncvar_get(ncin,'microzoo')
  var.5 <- v5/surface
  var.5[is.na(var.5)]<- NA
  tab.5[,,i]=var.5
}
assign(paste0("microzo.",j),(aaply(tab.5, c(1,2),mean, na.rm=T)))

# Mesozoo
tab.6 <- array(NA,c(567,264,7))
for(i in seq_along(files)) {
  ncin = nc_open(files[i])
  v6 = ncvar_get(ncin,'mesozoo')
  var.6 <- v6/surface
  var.6[is.na(var.6)]<- NA
  tab.6[,,i]=var.6
}
assign(paste0("mesozo.",j), (aaply(tab.6, c(1,2),mean, na.rm=T)))

# Regroupement des arrays en un seul
bb <- abind(picophy.24, nanophy.24, microphyto.24, nanozo.24, microzo.24, mesozo.24, along=3) # To change

#changer les coeficient de conversion de mmolC/m² à tonne/km²
cv1=0.1790541 # soit 179.0541 mgWW/mmolC ou 993.75 mgww/mmolN
cv2=0.1790541 
cv3=0.1790541
cv4=0.15 
cv5=0.15
cv6=0.15
# soit 150 mgWW/mmolC ou 832.5 mgww/mmolN

bbb<-abind(bb[,,1]*cv1,bb[,,2]*cv2,bb[,,3]*cv3,bb[,,4]*cv4,bb[,,5]*cv5, bb[,,6]*cv6, along = 3)

assign(paste0("plancton.",j),bbb)

###############################################################################################################################
############### Construction Netcdf par pas de temps ########################################################

wd <- "C:/Users/Fabien/Documents/Scripts R/Nemomed12-Eco3m/Data/with mean"
setwd(dir = wd)

# Définition des dimensions
xpos = ncdim_def("x","",1:567,unlim=FALSE) 
ypos = ncdim_def("y","",1:264,unlim=FALSE) 
temps = ncdim_def("time","",1,unlim=FALSE) 
ltl = ncdim_def("ltl","",1:6,unlim=FALSE) 

# Définition des variables
mv<--99
dx = ncvar_def("dx","",list(xpos,ypos),mv,longname="longitude")
dy = ncvar_def("dy","",list(xpos,ypos),mv,longname="latitude")
ltl_biomass = ncvar_def("ltl_biomass","",list(xpos,ypos,ltl, temps),mv,longname="ltl_biomass")

# Création du nouveau netcdf
b = nc_create("eco3m_med.24.nc",list(dx,dy,ltl_biomass)) # To change

# Remplissage des variables
setwd("C:/Users/Fabien/Documents/Scripts R/Nemomed12-Eco3m/Data/Nemo outputs by steps/1")
ncname <- "20060101-20060115-OsmoseInputs"
ncfname <- paste(ncname, ".nc", sep = "")
nc <- nc_open(ncfname)

#
ncvar_put(b, dx, ncvar_get(nc,varid="longitude_t"), verbose=TRUE)
ncvar_put(b, dy ,ncvar_get(nc,varid="latitude_t"), verbose=TRUE)
ncvar_put(b, ltl_biomass ,plancton.24, verbose=TRUE) # To change

# Fermer et sauvergarder le Netcdf
nc_close(b)
