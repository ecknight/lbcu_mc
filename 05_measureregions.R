library(tidyverse)
library(sf)
library(sp)
library(adehabitatHR)
library(raster)

#TO DO: FOR REGIONS, DECIDE ABOUT BOOTSTRAPS OR NOT####

#PART A: REGIONS####

#1. Set number of clusters----
n <- 2

#2. Get dominant cluster id for each bird----
clust <- read.csv("Data/LBCUKDEClusters.csv") %>% 
  dplyr::filter(nclust==n,
                !is.na(X)) %>% 
  group_by(id) %>% 
  summarize(kdecluster = round(mean(kdecluster))) %>% 
  ungroup()

#3. Import & wrangle daily locations----
#filter out migration
dat <- read.csv("Data/LBCU_FilteredData_Segmented.csv") %>% 
  dplyr::filter(id %in% dat$id,
                !(season %in% c("springmig", "fallmig") & stopover==0),
                winter!="wintermig") %>% 
  mutate(seasoncluster = ifelse(winter!="other", str_sub(winter, -1, -1), stopovercluster),
         seasoncluster = ifelse(is.na(seasoncluster), 0, as.numeric(seasoncluster)),
         kdeid = paste(season, id, year, seasoncluster, sep="-")) %>% 
  left_join(clust) %>% 
  group_by(kdeid) %>% 
  mutate(n=n()) %>% 
  ungroup() %>% 
  dplyr::filter(n >= 5)

#2. Make sp object----
dat.sp <- dat %>% 
  st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
  mutate(ID = paste(kdecluster, kdeid, sep="-")) %>% 
  dplyr::select(ID, geometry)

#3. Set up loop----
inds <- unique(dat.sp$ID)

kd.sp <- data.frame()
for(i in c(1:length(inds))){
  
  dat.i <- dat.sp %>% 
    dplyr::filter(ID==inds[i]) %>% 
    as_Spatial()
  
  #4. Calculate KDE----
  kd <- kernelUD(dat.i, grid = 1000, extent=2, h="href", same4all=FALSE)
  
  #5. Rasterize----
  kd.r <- raster(kd[[1]])
  kd.r.1 <- (kd.r-minValue(kd.r))/(maxValue(kd.r)-minValue(kd.r))
  raster::writeRaster(kd.r.1, paste0("gis/raster/kde_", inds[i], ".tif"), overwrite=TRUE)
  
  #6. Get shp of 95% isopleth----
  kd.sp.i <- try(getverticeshr(kd, percent=50) %>% 
    st_as_sf() %>% 
    st_transform(crs=4326))
  
  if(class(kd.sp.i)[1]=="sf"){
    kd.sp <- rbind(kd.sp, kd.sp.i)
  }
  else{
    file.remove(paste0("gis/raster/kde_", inds[i], ".tif"))
  }
  
  print(paste0("Finished individual ", i, " of ", length(inds)))
  
}

kd.out <- kd.sp %>% 
  mutate(area = round(area))

write_sf(kd.out, "gis/shp/kde_individual.shp")