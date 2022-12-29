library(tidyverse)
library(sf)
library(vegan)
library(ggmap)

options(scipen=9999)

#1. Read in data and filter out individuals that don't have locations for both seasons, individuals on atlantic coast
dat <- read.csv("Data/LBCU_FilteredData_Segmented.csv") %>% 
  dplyr::filter(!id %in% c(46768277, 33088, 129945787, 46770723, 46769927, 86872))

#1. Breeding ground means----
breed.mn <- dat %>% 
  dplyr::filter(season=="breed") %>% 
  group_by(id, year) %>% 
  summarize(lat = mean(lat),
            lon = mean(lon),
            days=n()) %>% 
  ungroup() %>% 
  mutate(season="breed",
         cluster=1)

#2. Wintering ground means----
winter.mn <- dat %>% 
  dplyr::filter(season=="winter", winter!="wintermig") %>% 
  mutate(year = ifelse(doy < 130, year-1, year),
         cluster = as.numeric(str_sub(winter, -1, -1))) %>% 
  group_by(id, year, cluster) %>% 
  summarize(lat = mean(lat),
            lon = mean(lon),
            days=n()) %>% 
  ungroup() %>% 
  mutate(season="winter")

#3. Stopover means----
stop.mn <- dat %>% 
  dplyr::filter(stopover==1) %>% 
  group_by(id, season, year, stopovercluster) %>% 
  summarize(lat = mean(lat),
            lon = mean(lon),
            days=n()) %>% 
  rename(cluster = stopovercluster) %>% 
  group_by(id, season, year) %>% 
  mutate(cluster = row_number()) %>% 
  ungroup()

#4. Put together----
mn <- rbind(breed.mn, winter.mn, stop.mn)
table(mn$season, mn$cluster)

#5. Get distance to coast-----
coast <- read_sf("gis/gshhg-shp-2.3.7/GSHHS_shp/l/GSHHS_l_L1.shp") %>% 
  st_make_valid() %>% 
  st_cast("LINESTRING")

mn.sf <- mn %>% 
  st_as_sf(coords=c("lon", "lat"), crs=4326) 

mn.near <- mn.sf %>% 
  st_nearest_feature(coast)

mn.coast <- data.frame(distance = as.numeric(st_distance(mn.sf, coast[mn.near,], by_element = TRUE))) %>% 
  cbind(mn)

#check
ggplot(mn.coast) +
  geom_point(aes(x=lon, y=lat, colour=distance))

#6. Transform to utm----
mn.utm <- st_as_sf(mn.coast, coords=c("lon", "lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(mn.coast)

write.csv(mn.utm, "Data/LBCUMCLocations.csv", row.names = FALSE)

#6. Visualize----
ggplot(mn.utm %>% dplyr::filter(season!="fallmig")) +
  geom_path(aes(x=X, y=Y, group=id)) +
  geom_point(aes(x=X, y=Y, colour=season))

ids <- unique(mn.utm$id)

for(i in 1:length(ids)){
  mn.i <- mn.utm %>% 
    dplyr::filter(id==ids[i]) %>% 
    arrange(year, season, cluster)
  
  ggplot(mn.i) +
    geom_path(aes(x=X, y=Y)) +
    geom_point(aes(x=X, y=Y, colour = season), size=3, alpha = 0.7) +
    scale_colour_viridis_d() +
    facet_wrap(~year)
  
   ggsave(filename=paste0("Figs/MC/", ids[i], ".jpeg"))
  
}
