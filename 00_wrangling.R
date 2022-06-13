library(tidyverse)
library(sf)
library(vegan)
library(data.table)
library(ggmap)

options(scipen=9999)

#1. Read in data and filter out individuals that don't have locations for both seasons
dat <- read.csv("Data/LBCU_FilteredData_Segmented.csv") %>% 
  dplyr::filter(!id %in% c(46768277, 33088, 129945787))

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

#3. Put together----
mn <- rbind(breed.mn, winter.mn, stop.mn)

mn.utm <- st_as_sf(mn, coords=c("lon", "lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(mn)

table(mn.utm$season, mn.utm$cluster)

write.csv(mn.utm, "Data/LBCUMCLocations.csv", row.names = FALSE)

#4. Visualize----
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
    geom_point(aes(x=X, y=Y, colour = days), size=3, alpha = 0.7) +
    scale_colour_viridis_c() +
    facet_wrap(~year)
  
  ggsave(filename=paste0("Figs/Days/", ids[i], ".jpeg"))
  
}
