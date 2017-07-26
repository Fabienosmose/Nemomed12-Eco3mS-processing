# Script to interpolate Nemo-Eco3m
# Author: Fabien Moullec
# Date: 22/02/2017

rm(list=ls())

# Load libraries
library(ncdf4)
library(raster)
library(rgdal)
library(gtools)
library(matlab)
library(RColorBrewer)
library(fields)
library(gtools)

# load raster de la mediterranée
rast.med <- readGDAL("C:/Users/Fabien/Documents/Scripts R/création des inputs osmose/grille osmose/rastmed.tif")
rast.med <- raster(rast.med)

# Check a NetCDF file 
setwd("~/Scripts R/Nemomed12-Eco3m/Data/with mean")
ncname <- "24-eco3m_med.24"
ncfname <- paste(ncname, ".nc", sep = "")
# open a NetCDF file
ncin <- nc_open(ncfname) ; print(ncin)
names(ncin$var)
print(paste("The file has",ncin$nvars,"variables,",ncin$ndims,"dimensions and",ncin$natts,"NetCDF attributes"))

# Récupérer les longitudes et latitudes
la <- ncvar_get(ncin, "dy"); lo <- ncvar_get(ncin, "dx")
image.plot(la);image.plot(lo)

la.vec <- as.vector(la)
lo.vec <- as.vector(lo)

# Création d'une boucle sur tous les netcdf (24 soit un par pas de temps) et sur les 6 types de plancton
tab <- array(NA,c(165,380,6,24))

#setwd("C:/Users/Fabien/Documents/Scripts R/Nemomed12-Eco3m/Data/with mean/")
files <- list.files("C:/Users/Fabien/Documents/Scripts R/Nemomed12-Eco3m/Data/with mean")
files <- mixedsort(files)# permet de trier les fichiers Netcdf dans l'ordre de 1 à 24

for (j in seq_along(files)){ 
  
  cat("j:",j,"\n" )
  
  ncin <- nc_open(files[j])
  
  dat <- array(NA,c(165,380,6))
  
  for (i in 1:6){
    plk = ncvar_get(nc=ncin, "ltl_biomass")[,,i]
    plk <- as.vector(plk)
    
    # Création d'un spatial dataframe
    df <- cbind(lo.vec, la.vec, plk)
    df <- as.data.frame(df)
    coordinates(df) <- ~ lo.vec+la.vec
    
    # create an empty raster
    reso <- 1/12
    rast_base <- raster(xmn=min(lo.vec),xmx=max(lo.vec),ymn=min(la.vec), ymx=max(la.vec), resolution= reso)
    
    # Rasterize nemo
    rast.bgc <- rasterize(df, rast_base,"plk", fun=mean)
    
    # Put projection of nemo as osmose
    proje <- "+init=epsg:3035 +proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
    rast.bgc.proj <- projectRaster(rast.bgc, crs = proje)
    
    # resample
    nemo.plk <- resample(rast.bgc.proj,rast.med, method ='bilinear')
    
    fab <- matrix(nemo.plk@data@values,165,380, byrow = T)
    fab_t <- fab[nrow(fab):1,] # to reverse the order of the rows
    dat[,,i]=fab_t
    
  }
  
  tab[,,,j]=dat
  
}

tabl <- aperm(tab, c(2,1,3,4))

#### To build a unique Netcdf file

setwd("C:/Users/Fabien/Documents/Scripts R/Nemomed12-Eco3m/Data") # Enregistre le Netcdf dans ce répertoire

# Définition des dimensions
xpos = ncdim_def("nx","",1:380,unlim=FALSE) 
ypos = ncdim_def("ny","",1:165,unlim=FALSE) 
temps = ncdim_def("Time","",1:24,unlim=FALSE) 
ltl = ncdim_def("ltl","",1:6,unlim=FALSE) 

# Définition des variables
# Attention dans les Netcdf, les dimensions sont inversées par rapport à R
mv<--99
dx = ncvar_def(name="longitude",units="",dim=list(xpos,ypos),missval=mv,longname="longitude")
dy = ncvar_def(name="latitude",units="",dim=list(xpos,ypos),missval=mv,longname="latitude")
ltl_biomass = ncvar_def("ltl_biomass","",list(xpos,ypos,ltl, temps),mv,longname="ltl_biomass")
time = ncvar_def("time", list(temps), longname = "time", units = "time step")

# Création du nouveau ncdf
b = nc_create("eco3m_med_v11.nc",list(dx,dy,ltl_biomass,time))

# Pour alimenter les variables dx et dy, on récupère les long et lat à partir du raster créé
(coordinates(nemo.plk)[,1])-> x
(coordinates(nemo.plk)[,2])-> y
ldx <- matrix(x,380,165, byrow = F)
ldy <- matrix(y,380,165, byrow = F)

# For time
ti <- 1:24

# Remplissage des variables
ncvar_put(b, dx, ldx, verbose = TRUE)
ncvar_put(b, dy ,ldy, verbose = TRUE)
ncvar_put(b, ltl_biomass ,tabl, verbose = TRUE)
ncvar_put(b, time, ti, verbose = TRUE)

# Fermer et sauvergarder le Netcdf
nc_close(b)

#######################################################################
######################## TO DOUBLE CHECK ##############################
#######################################################################

# Check a NetCDF file 
setwd("C:/Users/Fabien/Documents/Scripts R/Nemomed12-Eco3m/Data")
ncname <- "eco3m_med_v11"
ncfname <- paste(ncname, ".nc", sep = "")
# open a NetCDF file
ncin <- nc_open(ncfname) ; print(ncin)
names(ncin$var)
print(paste("The file has",ncin$nvars,"variables,",ncin$ndims,"dimensions and",ncin$natts,"NetCDF attributes"))


latitude = ncvar_get(ncin, "latitude")
latitude <- as.vector(latitude)
longitude = ncvar_get(ncin, "longitude")
longitude <- as.vector(longitude)

plk = ncvar_get(ncin, "ltl_biomass")[,,3,10]
plk <- flipud(t(plk))

projec <- "+init=epsg:3035 +proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
r=raster(plk,  xmn=min(longitude)-5000, xmx=max(longitude)+5000, ymn=min(latitude)-5000, ymx=max(latitude)+5000,crs=projec)
# rast <-projectRaster(r, crs = projec)
# rainbow color scheme
cols = rev(colorRampPalette(brewer.pal(11, 'Spectral'))(255))
#setting margins for plotting
par(mar=c(2,2,1,1))
#plot(r, col=cols, main='picophy time step 1')
plot(r, col=tim.colors(10000), axes=T)

