library(tidyverse)
library(ppclust)
library(sf)
library(data.table)

options(scipen=9999)

#1. Import tracking data with training clusters----
track.raw <- read.csv("Data/LBCUKDEClusters.csv") %>% 
  mutate(type = "track") 

#2. Import bbs monitoring data----
load("~/Library/Application Support/bbsBayes/bbs_raw_data.RData")

bbs <- bbs_data[["bird"]] %>% 
  dplyr::filter(AOU==2640)

routes <- bbs %>% 
  dplyr::select(RouteDataID) %>%
  unique() %>% 
  left_join(bbs_data[["route"]]) %>% 
  group_by(countrynum, statenum, Route, Latitude, Longitude) %>%
  summarize(n = n()) %>% 
  dplyr::filter(!is.na(Latitude)) %>% 
  ungroup()

bbs.sf <- routes %>% 
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(routes) %>% 
  mutate(id=paste(countrynum, statenum, Route, sep="-"),
         type="bbs")

#3. Set up loop through # of clusters---
clusters <- c(2:10)

fcm.out <- list()
for(j in 1:length(clusters)){
  
  #4. Wrangle tracking data----
  track.j <- track.raw %>% 
    dplyr::filter(nclust==clusters[j],
                  season=="breed")
  
  #5. Put together with BBS data----
  dat.j <- track.j %>% 
    dplyr::select(id, X, Y, type) %>% 
    rbind(bbs.sf %>% 
            dplyr::select(id, X, Y, type))
  
  #6. Create membership matrix with known membership for tracking data and unknown for bbs----
  mem.j <- track.j %>% 
    mutate(value = 1) %>% 
    pivot_wider(id_cols=id, names_from=kdecluster, values_from=value, values_fill=0) %>% 
    dplyr::select(-id) %>% 
    data.frame() %>% 
    rbind(expand.grid(id=bbs.sf$id, kdecluster=rep(seq(1, clusters[j], 1))) %>% 
            mutate(value=1/clusters[j]) %>% 
            pivot_wider(id_cols=id, names_from=kdecluster, values_from=value) %>% 
            dplyr::select(-id) %>% 
            data.frame())

  #7. Fuzzy c-means cluster with cluster id as initial membership #----
  fcm.j <- fcm(dat.j[,c("X", "Y")], centers=clusters[j], memberships=mem.j)

  #8. Keep fcm membership of each BBS route and tracking data point----
  fcm.out[[j]] <- dat.j %>% 
    mutate(fcmcluster=fcm.j$cluster,
           nclust=clusters[j])
  
  print(paste0("Finished cluster ", clusters[j], " of ", length(clusters)))
  
}

#9. Collapse the results----
fcm.all <- rbindlist(fcm.out)

#10. Check for agreement with kde for tracking data----
fcm.track <- fcm.all %>% 
  dplyr::filter(type=="track") %>% 
  left_join(track.raw %>% 
              dplyr::filter(season=="breed") %>% 
              mutate(id=as.character(id)))

cor.test(fcm.track$fcmcluster, fcm.track$kdecluster)

ggplot(fcm.track) +
  geom_point(aes(x=kdecluster, y=fcmcluster)) +
  facet_wrap(~nclust)

#11. Visualize----
ggplot(fcm.all) +
  geom_point(aes(x=X, y=Y, colour=factor(fcmcluster), pch=type)) +
  geom_point(data=fcm.track, aes(x=X, y=Y, colour=factor(kdecluster)), size=5, alpha=0.5) +
  facet_wrap(~nclust)

#NOPE THIS DIDN"T WORK AT ALL. REVISIT CLASSIFICATION

#12. Save out----
fcm <- rbindlist(fcm.out)