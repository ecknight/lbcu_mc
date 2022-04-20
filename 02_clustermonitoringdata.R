library(tidyverse)
library(ppclust)

options(scipen=9999)

#1. Import tracking data with training clusters----
dat.raw <- read.csv("Data/LBCUKDEClusters.csv") 

#2. Import bbs monitoring data----
load("~/Library/Application Support/bbsBayes/bbs_raw_data.RData")
bbs <- bbs_data[["bird"]] %>% 
  dplyr::filter(AOU==2640)
routes <- dat %>% 
  dplyr::select(RouteDataID) %>%
  unique() %>% 
  left_join(bbs_data[["route"]]) %>% 
  group_by(countrynum, statenum, Route, Latitude, Longitude) %>%
  summarize(n = n()) %>% 
  dplyr::filter(!is.na(Latitude)) %>% 
  ungroup()

bbs.sf <- routes %>% 
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326) %>% 
  st_transform(crs=3857)


#3. Set up loop through # of clusters---
clusters <- c(2:10)

for(j in 1:length(clusters)){
  
  #4. Wrangle----
  dat.j <- dat.raw %>% 
    dplyr::filter(nclust==clusters[j],
                  season=="breed") 
  mem.j <- dat.j %>% 
    mutate(value = 1) %>% 
    pivot_wider(id_cols=id, names_from=kdecluster, values_from=value, values_fill=0) %>% 
    dplyr::select(-id) %>% 
    data.frame()
  
  #5. Fuzzy c-means cluster with cluster id as initial membership #----
  fcm.j <- fcm(dat.j[,c("X", "Y")], centers=clusters[j], memberships=mem.j)
  
  #6. Apply model to BBS routes----
  
  
  #7. Save out membership of each BBS route----
  
}

