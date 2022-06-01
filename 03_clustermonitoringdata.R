library(tidyverse)
library(class)
library(sf)
library(data.table)
library(adehabitatHR)

options(scipen=9999)

#TO DO: FIGURE OUT KNN HYPERPARAMETER####

#1. Import tracking data with training clusters----
track.raw <- read.csv("Data/LBCUKDEClusters.csv")

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

#3. Set # of clusters---
clusters <- 3

#4. Set up bootstrap loop----
boot <- max(track.raw$boot)

knn.out <- list()
for(i in 1:boot){
  
  #5. Wrangle tracking data----
  track.i <- track.raw %>% 
    dplyr::filter(nclust==clusters,
                  boot==i,
                  season=="breed")
  
  #6. Put together with BBS data----
  bbs.i <- bbs.sf %>% 
    dplyr::select(id, X, Y, type)
  
  #7. KNN with cluster id as initial membership #----
  knn.i <- knn(train=track.i[,c("X", "Y")], test=bbs.i[,c("X", "Y")], cl=track.i$kdecluster, k=1, prob=TRUE)
  
  #8. Keep knn membership of each BBS route and tracking data point----
  knn.out[[i]] <- data.frame(knncluster=knn.i,
                             knnprob=attr(knn.i, "prob"),
                             boot=i) %>% 
    cbind(bbs.i)
  
  print(paste0("Finished bootstrap ", i, " of ", boot))
  
}

#9. Collapse results----
knn.all <- rbindlist(knn.out)

#10. Look at variation in cluster membership----
knn.sum <- knn.all %>% 
  group_by(id, knncluster) %>% 
  summarize(n=n()) %>% 
  ungroup()
summary(knn.sum$n)

knn.99 <- knn.sum %>% 
  dplyr::filter(n < 100)
table(knn.99$id)
#Only 46 points that have 1 instance of cluster variation

#11. Pick mean dominant cluster ID for each route----
knn.final <- knn.all %>% 
  group_by(id, X, Y, type) %>% 
  summarize(knncluster = round(mean(as.numeric(knncluster)))) %>% 
  ungroup()

#14. Visualize----
track.final <- track.raw %>% 
  dplyr::filter(season=="breed", nclust==3, boot==1)

ggplot(knn.final) +
  geom_point(aes(x=X, y=Y, colour=factor(knncluster))) +
  geom_point(data=track.final, aes(x=X, y=Y, colour=factor(kdecluster)), pch=21, fill="white")

#15. Calculate MCP area for BBSbayes-----
sp <- SpatialPointsDataFrame(coords=cbind(knn.final$X, knn.final$Y), 
                               data=data.frame(ID=knn.final$knncluster),
                               proj4string = CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"))

mp <- mcp(sp[,1], percent=100)

knn.area <- data.frame(mp) %>% 
  rename(knncluster=id)

#15. Save out----
write.csv(knn.all, "Data/LBCUBBSClustersAll.csv", row.names = FALSE)
write.csv(knn.area, "Data/area_weight.csv", row.names = FALSE)
write.csv(knn.final, "Data/LBCUBBSClusters.csv", row.names = FALSE)
