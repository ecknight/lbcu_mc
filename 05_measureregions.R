library(tidyverse)
library(sf)
library(sp)
library(adehabitatHR)
library(raster)

#TO DO: FOR REGIONS, DECIDE ABOUT BOOTSTRAPS OR NOT####

#PART A: REGIONS####

#1. Set number of clusters----
n <- 3

#2. Read in clustered tracking data and assign----
dat <- read.csv("Data/LBCUKDEClusters.csv") %>% 
  dplyr::filter(nclust==n,
                boot==1,
                !is.na(X))

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
  raster::writeRaster(kd.r, paste0("gis/region/kde_", list[i], ".tif"), overwrite=TRUE)
}

#6. Get shp of 95% isopleth----
kd.sp <- getverticeshr(kd, percent=95) %>% 
  st_as_sf() %>% 
  st_transform(crs=4326)
write_sf(kd.sp, "gis/kde.shp")

#PART B: INDIVIDUALS####

#1. Import & wrangle data----
dat.ind <- read.csv("Data/LBCU_FilteredData_Segmented.csv") %>% 
  dplyr::filter(id %in% dat$id,
                !(season %in% c("springmig", "fallmig") & stopover==0),
                winter!="wintermig") %>% 
  mutate(seasoncluster = ifelse(winter!="other", str_sub(winter, -1, -1), stopovercluster),
         seasoncluster = ifelse(is.na(seasoncluster), 0, as.numeric(seasoncluster)),
         kdeid = paste(id, season, seasoncluster, sep="-")) %>% 
  left_join(dat %>% 
              dplyr::select(id, kdecluster) %>% 
              unique()) %>% 
  group_by(kdeid) %>% 
  mutate(n=n()) %>% 
  ungroup() %>% 
  dplyr::filter(n >= 5)

#2. Make sp object----
dat.ind.sp <- dat.ind %>% 
  st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
  mutate(ID = paste(kdecluster, kdeid, sep="-")) %>% 
  dplyr::select(ID, geometry) %>% 
  as_Spatial()

#3. Calculate KDE----
kd.ind <- kernelUD(dat.ind.sp, grid = 100, extent=2, h="href", same4all=FALSE)

#5. Rasterize----
list <- names(kd.ind)
for(i in 1:length(list)){
  kd.r <- raster(kd.ind[[i]])
  raster::writeRaster(kd.r, paste0("gis/individual/kde_", list[i], ".tif"), overwrite=TRUE)
}
