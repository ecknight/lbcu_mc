library(tidyverse)
library(sf)
library(vegan)
library(data.table)
library(ggmap)

options(scipen=9999)

#NOTE: This ignores birds with multiple wintering grounds and variation between years for now. We're just interested in broad patterns for exploration

#1. Read in data and filter out individuals that don't have locations for both seasons
dat <- read.csv("Data/LBCU_FilteredData_Segmented.csv") %>% 
  dplyr::filter(!id %in% c(46768277))

#1. Breeding ground means----
breed.mn <- dat %>% 
  dplyr::filter(season=="breed") %>% 
  group_by(id) %>% 
  summarize(lat = mean(lat),
            lon = mean(lon)) %>% 
  ungroup() %>% 
  mutate(season="breed",
         cluster=1) 

#2. Wintering ground means----
winter.mn <- dat %>% 
  dplyr::filter(season=="winter") %>% 
  group_by(id) %>% 
  summarize(lat = mean(lat),
            lon = mean(lon)) %>% 
  ungroup() %>% 
  mutate(season="winter",
         cluster=1)

#3. Stopover means----
stop.mn <- dat %>% 
  dplyr::filter(stopover==1) %>% 
  group_by(id, season, stopovercluster) %>% 
  summarize(lat = mean(lat),
            lon = mean(lon)) %>% 
  ungroup() %>% 
  rename(cluster = stopovercluster)

#3. Put together----
mn <- rbind(breed.mn, winter.mn, stop.mn)

mn.utm <- st_as_sf(mn, coords=c("lon", "lat"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(mn)

write.csv(mn.utm, "Data/LBCUMCLocations.csv", row.names = FALSE)

#4. Visualize----
ggplot(mn.utm) +
  geom_path(aes(x=X, y=Y, group=id)) +
  geom_point(aes(x=X, y=Y, colour=season))