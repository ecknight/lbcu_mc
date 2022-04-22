library(tidyverse)
library(class)
library(sf)
library(data.table)
library(adehabitatHR)

options(scipen=9999)

#TO DO: FIGURE OUT KNN HYPERPARAMETER####

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
clusters <- c(2:6,8:9)

knn.out <- list()
knn.area <- data.frame()
for(j in 1:length(clusters)){
  
  #4. Wrangle tracking data----
  track.j <- track.raw %>% 
    dplyr::filter(nclust==clusters[j],
                  season=="breed")
  
  #5. Put together with BBS data----
  bbs.j <- bbs.sf %>% 
            dplyr::select(id, X, Y, type)

  #7. Fuzzy c-means cluster with cluster id as initial membership #----
  knn.j <- knn(train=track.j[,c("X", "Y")], test=bbs.j[,c("X", "Y")], cl=track.j$kdecluster, k=1, prob=TRUE)

  #8. Keep knn membership of each BBS route and tracking data point----
  knn.out[[j]] <- data.frame(knncluster=knn.j,
                        knnprob=attr(knn.j, "prob"),
                        nclust=clusters[j]) %>% 
    cbind(bbs.j)
  
  #9. Calculate MCP area for BBSbayes-----
  sp.j <- SpatialPointsDataFrame(coords=cbind(bbs.j$X, bbs.j$Y), 
                                   data=data.frame(ID=knn.j),
                                   proj4string = CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"))
  
  mp.j <- mcp(sp.j[,1], percent=100)
  
  knn.area <- data.frame(mp.j) %>% 
    rename(knncluster=id) %>% 
    mutate(nclust=clusters[j]) %>% 
    rbind(knn.area)
  
  print(paste0("Finished cluster ", clusters[j], " of ", length(clusters)))
  
}

#10. Collapse the results----
knn.all <- rbindlist(knn.out)

#11. Visualize----
ggplot(knn.all) +
  geom_point(aes(x=X, y=Y, colour=factor(knncluster))) +
  geom_point(data=filter(track.raw, season=="breed"), aes(x=X, y=Y, colour=factor(kdecluster)), pch=21, fill="white") +
  facet_wrap(~nclust)

#12. Save out----
write.csv(knn.all, "LBCUBBSClusters.csv", row.names = FALSE)
write.csv(knn.area, "area_weight.csv", row.names = FALSE)
