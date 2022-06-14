library(tidyverse)
library(sf)
library(sp)
library(adehabitatHR)
library(raster)

#TO DO: FOR REGIONS, DECIDE ABOUT BOOTSTRAPS OR NOT####

#PART A: REGIONS####

#1. Set number of clusters----
n <- 2

#2. Read in clustered tracking data and assign----
dat <- read.csv("Data/LBCUKDEClusters.csv") %>% 
  dplyr::filter(nclust==n,
                boot==1,
                !is.na(X))

#3. Make sp object----
dat.sp <- dat %>% 
  st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
  mutate(ID = paste0(kdecluster, "-", season, "-region")) %>% 
  dplyr::select(ID, geometry) %>% 
  as_Spatial()

#4. Calculate KDE----
kd <- kernelUD(dat.sp, grid = 1000, extent=2, h="href", same4all=FALSE)

#5. Rasterize----
list <- names(kd)
for(i in 1:length(list)){
  kd.r <- raster(kd[[i]])
  kd.r.1 <- (kd.r-minValue(kd.r))/(maxValue(kd.r)-minValue(kd.r))
  raster::writeRaster(kd.r.1, paste0("gis/raster/kde_", list[i], ".tif"), overwrite=TRUE)
}

#6. Get shp of 95% isopleth----
kd.sp <- getverticeshr(kd, percent=95) %>% 
  st_as_sf() %>% 
  st_transform(crs=4326) %>% 
  dplyr::select(-area)

write_sf(kd.sp, "gis/shp/kde_region.shp")

#PART B: INDIVIDUALS####

#1. Import & wrangle data----
dat.ind <- read.csv("Data/LBCU_FilteredData_Segmented.csv") %>% 
  dplyr::filter(id %in% dat$id,
                !(season %in% c("springmig", "fallmig") & stopover==0),
                winter!="wintermig") %>% 
  mutate(seasoncluster = ifelse(winter!="other", str_sub(winter, -1, -1), stopovercluster),
         seasoncluster = ifelse(is.na(seasoncluster), 0, as.numeric(seasoncluster)),
         kdeid = paste(season, id, seasoncluster, sep="-")) %>% 
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
  dplyr::select(ID, geometry)

#3. Set up loop----
inds <- unique(dat.ind.sp$ID)

kd.ind.sp <- data.frame()
for(i in c(1:10)){
  
  dat.i <- dat.ind.sp %>% 
    dplyr::filter(ID==inds[i]) %>% 
    as_Spatial()
  
  #4. Calculate KDE----
  kd.ind <- kernelUD(dat.i, grid = 1000, extent=2, h="href", same4all=FALSE)
  
  #5. Rasterize----
  kd.r <- raster(kd.ind[[1]])
  kd.r.1 <- (kd.r-minValue(kd.r))/(maxValue(kd.r)-minValue(kd.r))
  raster::writeRaster(kd.r.1, paste0("gis/raster/kde_", inds[i], ".tif"), overwrite=TRUE)
  
  #6. Get shp of 95% isopleth----
  kd.sp.i <- try(getverticeshr(kd.ind, percent=95) %>% 
    st_as_sf() %>% 
    st_transform(crs=4326))
  
  if(class(kd.sp.i)[1]=="sf"){
    kd.ind.sp <- rbind(kd.ind.sp, kd.sp.i)
  }
  else{
    file.remove(paste0("gis/raster/kde_", inds[i], ".tif"))
  }
  
  print(paste0("Finished individual ", i, " of ", length(inds)))
  
}

write_sf(kd.ind.sp, "gis/shp/kde_individual.shp")

#7. Put the 95% isopleth shp together----
kd.all.sp <- rbind(kd.ind.sp, kd.sp)
write_sf(kd.all.sp, "gis/shp/kde.shp")