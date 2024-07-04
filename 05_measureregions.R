library(tidyverse)
library(sf)
library(sp)
library(adehabitatHR)
library(raster)
library(lme4)
library(MuMIn)

#A. INDIVIDUAL HOME RANGES####

#1. Set number of clusters----
n <- c("3", "expert", "flyway")

#2. Get dominant cluster id for each bird----
#mean works because it's always just between 2 clusters
clust <- read.csv("Data/LBCUKDEClusters.csv") %>% 
  dplyr::filter(nclust %in% n,
                !is.na(X)) |> 
  dplyr::select(id, year, nclust, group) |> 
  unique()

#3. Import & wrangle daily locations----
#filter out migration (not stopovers)
#take out some egregious misclassifications to fix overestimates of home range
dat <- read.csv("Data/LBCU_FilteredData_Segmented.csv") %>% 
  dplyr::filter(id %in% clust$id,
                !(season %in% c("springmig", "fallmig") & stopover==0),
                winter!="wintermig") %>% 
  mutate(seasoncluster = ifelse(winter!="other", str_sub(winter, -1, -1), stopovercluster),
         seasoncluster = ifelse(is.na(seasoncluster), 0, as.numeric(seasoncluster)),
         kdeid = paste(season, id, year, seasoncluster, sep="-")) %>% 
  inner_join(clust, multiple="all") %>% 
  group_by(nclust, kdeid) %>% 
  mutate(n=n()) %>% 
  ungroup() %>% 
  dplyr::filter(n >= 5) %>% 
  dplyr::filter(!(id==1418951627 & year==2020 & season=="breed" & X > -12030000 & Y < 6050000),
                !(id==46768108 & year==2015 & season=="breed" & Y < 5800000))

#4. Make sp object----
dat.sp <- dat %>% 
  st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
  mutate(ID = paste(nclust, group, kdeid, sep="-")) %>% 
  dplyr::select(ID, geometry)

#5. Set up loop----
inds <- unique(dat.sp$ID)

kd.sp <- data.frame()
for(i in c(1:length(inds))){
  
  dat.i <- dat.sp %>% 
    dplyr::filter(ID==inds[i]) %>% 
    as_Spatial()
  
  #6. Calculate KDE----
  kd <- kernelUD(dat.i, grid = 1000, extent=2, h="href", same4all=FALSE)
  
  #7. Rasterize----
  kd.r <- raster(kd[[1]])
  kd.r.1 <- (kd.r-minValue(kd.r))/(maxValue(kd.r)-minValue(kd.r))
  raster::writeRaster(kd.r.1, paste0("gis/raster_ind/kde_", inds[i], ".tif"), overwrite=TRUE)
  
  #8. Get shp of 50% isopleth----
  kd.sp.i <- try(getverticeshr(kd, percent=50) %>% 
    st_as_sf() %>% 
    st_transform(crs=4326))
  
  if(class(kd.sp.i)[1]=="sf"){
    kd.sp <- rbind(kd.sp, kd.sp.i)
  }
  else{
    file.remove(paste0("gis/raster_ind/kde_", inds[i], ".tif"))
  }
  
  print(paste0("Finished individual ", i, " of ", length(inds)))
  
}

#9. save----
kd.out <- kd.sp %>% 
  mutate(area = round(area)/100, 
         rad = round(sqrt(area/pi), 1)) %>% 
  separate(id, into=c("nclust", "group", "season", "bird", "year", "cluster"), remove=FALSE)

write_sf(kd.out, "gis/shp/kde_individual.shp")
kd.out <- read_sf("gis/shp/kde_individual.shp")

#10. Seasonal averages----
kd.sum <- kd.out %>% 
  dplyr::select(-geometry) %>% 
  data.frame() %>% 
#  group_by(group, season) %>% 
#  group_by(season) %>% 
  summarize(area.mn = mean(area), 
            area.sd = sd(area)) %>% 
#  ungroup() %>% 
  mutate(rad = sqrt(area.mn/pi))
kd.sum

#11. Peak at effects of cluster & season on homerange area----
l1 <- lme4::lmer(area ~ group*season + (1|bird), kd.out, na.action = "na.fail")
MuMIn::dredge(l1)
summary(l1)

ggplot(filter(kd.out, area < 10000)) +
  geom_boxplot(aes(x=season, y=log(area), fill=factor(group))) +
  facet_wrap(~nclust)

#B. REGIONS####

#1. Set up loop----
regs <- dat %>% 
  dplyr::select(nclust, group, season) %>% 
  unique() %>% 
  mutate(ID = paste(nclust, group, season))

isos <- seq(5, 95, 5)

kd.sp <- data.frame()
for(i in c(1:nrow(regs))){
  
  #2. Get centroid for each bird----
  dat.i <- dat %>% 
    inner_join(regs[i,]) %>% 
    mutate(ID = paste(nclust, group, season)) %>% 
    group_by(ID, id) %>% 
    summarize(X=mean(X),
              Y=mean(Y)) %>% 
    ungroup() %>% 
    st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
    dplyr::select(-id) %>% 
    as_Spatial()
  
  #3. Calculate KDE----
  kd <- kernelUD(dat.i, grid = 1000, extent=2, h="href", same4all=FALSE)
  
  #4. Rasterize----
  kd.r <- raster(kd[[1]])
  kd.r.1 <- (kd.r-minValue(kd.r))/(maxValue(kd.r)-minValue(kd.r))
  raster::writeRaster(kd.r.1, paste0("gis/raster_reg/kde_", regs$ID[i], ".tif"), overwrite=TRUE)
  
  #5. Get shp of isopleths----
  for(j in 1:length(isos)){
    
    kd.sp.i <- try(getverticeshr(kd, percent=isos[j]) %>% 
                     st_as_sf() %>% 
                     st_transform(crs=4326))
    
    if(class(kd.sp.i)[1]=="sf"){
      kd.sp <- rbind(kd.sp, kd.sp.i %>% 
                       mutate(iso = isos[j]))
    }
    else{
      file.remove(paste0("gis/raster/kde_", regs[i], ".tif"))
    }
    
  }

  print(paste0("Finished region ", i, " of ", nrow(regs)))
  
}

#6. save----
kd.out <- kd.sp %>% 
  mutate(area = round(area)/100, 
         rad = round(sqrt(area/pi), 1)) %>% 
  separate(id, into=c("nclust", "group", "season"), remove=TRUE)

write_sf(kd.out, "gis/shp/kde_region.shp")
