library(tidyverse)
library(class)
library(sf)
library(data.table)
library(adehabitatHR)
library(rpart)
library(party)
library(klaR)
library(caret)

options(scipen=9999)

#A. PREAMBLE####

#1. Import tracking data with training clusters----
track.raw <- read.csv("Data/LBCUKDEClusters.csv")

#2. Import bbs monitoring data----
load("~/Library/Application Support/bbsBayes/bbs_raw_data.RData")

bbs <- bbs_data[["bird"]] %>% 
  dplyr::filter(AOU==2640,
                Year >= 1970)

routes <- bbs %>% 
  dplyr::select(RouteDataID) %>%
  unique() %>% 
  left_join(bbs_data[["route"]]) %>% 
  group_by(countrynum, statenum, Route, Latitude, Longitude) %>%
  summarize(n = n()) %>% 
  dplyr::filter(!is.na(Latitude)) %>% 
  ungroup()

bbs.utm <- routes %>% 
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(routes) %>% 
  mutate(id=paste(countrynum, statenum, Route, sep="-"),
         type="bbs")

#3. Determine if tracking data near enough to bbs route----

#calculate distance between tracking points and bbs routes
track.sf <- track.raw %>% 
  dplyr::filter(season=="breed") %>% 
  dplyr::select(lat, lon) %>% 
  unique() %>% 
  st_as_sf(coords=c("lon", "lat"), crs=4326) 

bbs.sf <- routes %>% 
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326)

bbs.near <- bbs.sf %>% 
  st_nearest_feature(track.sf)

bbs.dist <- data.frame(distance = as.numeric(st_distance(bbs.sf, track.sf[bbs.near,], by_element = TRUE))) %>% 
  cbind(bbs.utm) %>% 
  mutate(distance = distance/1000)

#visualize
ggplot() + 
  geom_point(data=bbs.dist, aes(x=Longitude, y=Latitude, colour=distance)) +
  geom_point(data=dplyr::filter(track.raw, season=="breed"), aes(x=lon, y=lat), size=3, colour="grey") + 
  scale_colour_viridis_c()

hist(bbs.dist$distance)

#remove points < 95% (500 km)
bbs.use <- bbs.dist %>% 
  dplyr::filter(distance < quantile(distance, 0.95))

#visualize again
ggplot() + 
  geom_point(data=bbs.use, aes(x=Longitude, y=Latitude, colour=distance)) +
  geom_point(data=dplyr::filter(track.raw, season=="breed"), aes(x=lon, y=lat), size=3, colour="grey") + 
  scale_colour_viridis_c()

#save for trend estimation
write.csv(bbs.use, "Data/BBSRoutesToUse.csv", row.names=FALSE)

#4. Set # of clusters---
clusts <- c("3", "expert", "flyway")

#B. PREDICT TO BBS DATA####

#1. Set up cluster loop----
knn.out <- list()
for(i in 1:length(clusts)){
  
  set.seed(i)
  
  #2. Wrangle tracking data----
  track.i <- track.raw %>% 
    dplyr::filter(nclust==clusts[i],
                  season=="breed")
  
  #3. Put together with BBS data----
  bbs.i <- bbs.use %>% 
    dplyr::select(id, X, Y, type)
  
  #4. KNN----
  knn.i <- knn(train=track.i[,c("X", "Y")], test=bbs.i[,c("X", "Y")], cl=track.i$group, k=10, prob=TRUE)
  
  knn.out[[i]] <- data.frame(knncluster=knn.i,
                             knnprob=attr(knn.i, "prob"),
                             nclust=clusts[i]) %>% 
    cbind(bbs.i)
  
  print(paste0("Finished cluster ", i, " of ", length(clusts)))
  
}

#5. Collapse results----
knn.all <- rbindlist(knn.out)

#6. Visualize----
track.final <- track.raw %>% 
  dplyr::filter(season=="breed",
                nclust %in% clusts)

ggplot(knn.all) +
  geom_point(aes(x=X, y=Y, colour=factor(knncluster))) +
  geom_point(data=track.final, aes(x=X, y=Y, colour=factor(group)), pch=21, fill="white", size=3) +
  facet_wrap(~nclust)

#9. Calculate MCP area for BBSbayes-----
sp <- SpatialPointsDataFrame(coords=cbind(knn.all$X, knn.all$Y), 
                               data=data.frame(ID=paste0(knn.all$nclust, "_", knn.all$knncluster)),
                               proj4string = CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"))

mp <- mcp(sp[,1], percent=100)

knn.area <- data.frame(mp) %>% 
  separate(id, into=c("knncluster", "nclust"))

#15. Save out----
write.csv(knn.area, "Data/area_weight.csv", row.names = FALSE)
write.csv(knn.all, "Data/LBCUBBSClusters.csv", row.names = FALSE)
