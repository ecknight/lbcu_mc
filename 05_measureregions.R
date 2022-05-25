library(tidyverse)
library(sf)
library(sp)
library(adehabitatHR)
library(raster)

#1. Set number of clusters----
n <- 4

#2. Read in clustered tracking data----
dat <- read.csv("Data/LBCUKDEClusters.csv") %>% 
  dplyr::filter(nclust==n, id!=129945787)

#3. Make sp object----
dat.sp <- dat %>% 
  st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
  mutate(ID = paste0(season, "-", kdecluster)) %>% 
  dplyr::select(ID, geometry) %>% 
  as_Spatial()

#4. Calculate KDE----
kd <- kernelUD(dat.sp, grid = 1000, extent=2, h="href", same4all=FALSE)

#5. Rasterize----
list <- names(kd)
for(i in 1:length(list)){
  kd.r <- raster(kd[[i]])
  raster::writeRaster(kd.r, paste0("gis/kde_", list[i], ".tif"), overwrite=TRUE)
}