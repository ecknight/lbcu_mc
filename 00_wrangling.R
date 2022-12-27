library(tidyverse)
library(sf)
library(vegan)
library(ggmap)
library(rgee)

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

#5. Get elevation----
ee_Initialize(gcs=TRUE)
ee_check()

mn.sf <- st_as_sf(mn, coords=c("lon", "lat"), crs=4326)
mn.ee <- sf_as_ee(mn.sf)

dem <- ee$Image("USGS/GMTED2010")

mn.dem <- ee_extract(
  x = dem,
  y = mn.ee,
  sf = FALSE
)

#6. Get ecoregion----
#Shapefile from https://www.epa.gov/eco-research/ecoregions

eco <- read_sf("gis/na_cec_eco_l1/NA_CEC_Eco_Level1.shp")

mn.eco <- mn.sf %>% 
  st_transform(raster::crs(eco)) %>% 
  st_intersection(eco) %>% 
  data.frame() %>% 
  dplyr::select(-NA_L1KEY, -NA_L1NAME, -geometry, -Shape_Leng, -Shape_Area) %>% 
  rename(ecoreg = NA_L1CODE)

#7. Get distance to coast-----
coast <- read_sf("gis/gshhg-shp-2.3.7/GSHHS_shp/l/GSHHS_l_L1.shp") %>% 
  st_make_valid() %>% 
  st_cast("LINESTRING")

mn.near <- st_nearest_feature(mn.sf, coast)
mn.coast <- data.frame(distance = as.numeric(st_distance(mn.sf, coast[mn.near,], by_element = TRUE))) %>% 
  cbind(mn)

ggplot(mn.coast) +
  geom_point(aes(x=lon, y=lat, colour=distance))

#8. Transform to utm----
mn.utm <- st_as_sf(mn, coords=c("lon", "lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(mn) %>% 
  left_join(mn.dem) %>% 
  left_join(mn.coast) %>% 
  rename(elevation = be75)

write.csv(mn.utm, "Data/LBCUMCLocations.csv", row.names = FALSE)

#9. Visualize----
ggplot(mn.utm) +
  geom_path(aes(x=X, y=Y, group=id)) +
  geom_point(aes(x=X, y=Y, colour=season))

ids <- unique(mn.utm$id)

for(i in 1:length(ids)){
  mn.i <- mn.utm %>% 
    dplyr::filter(id==ids[i]) %>% 
    arrange(year, season, cluster)
  
  ggplot(mn.i) +
#    geom_path(aes(x=X, y=Y)) +
    geom_point(aes(x=X, y=Y, colour = season), size=3, alpha = 0.7) +
    scale_colour_viridis_d() +
    facet_wrap(~year)
  
#   ggsave(filename=paste0("Figs/MC/", ids[i], ".jpeg"))
  
}
