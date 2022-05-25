library(tidyverse)
library(sf)
library(sp)
library(adehabitatHR)

#1. Set number of clusters----
n <- 4

#2. Read in clustered tracking data----
dat <- read.csv("Data/LBCUKDEClusters.csv") %>% 
  dplyr::filter(nclust==n, id!=129945787)

#3. Make sp object----
dat.sp <- dat %>% 
  st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
  mutate(ID = paste0(season, "-", kdecluster)) %>% 
  dplyr::select(ID, geometry) %>% 
  as_Spatial()

#4. Calculate KDE----
kd <- kernelUD(dat.sp, grid = 1000, extent=2, h="href", same4all=FALSE)

#5. Calculate area----
kd.95 <- getverticeshr(kd, percent=95, unin="m", unout="km2") %>% 
  st_as_sf() %>% 
  mutate(rad = sqrt(area/3.14)) %>% 
  separate(id, into=c("season", "kdecluster"), remove=FALSE) %>% 
  mutate(kdecluster =as.numeric(kdecluster))
kd.50 <- getverticeshr(kd, percent=50, unin="m", unout="km2") %>% 
  st_as_sf() %>% 
  mutate(rad = sqrt(area/3.14)) %>% 
  separate(id, into=c("season", "kdecluster"), remove=FALSE) %>% 
  mutate(kdecluster =as.numeric(kdecluster))

mcp.95 <- mcp(dat.sp, percent=95, unin="m", unout="km2") %>% 
  st_as_sf() %>% 
  separate(id, into=c("season", "kdecluster"), remove=FALSE) %>% 
  mutate(kdecluster =as.numeric(kdecluster))
mcp.75 <- mcp(dat.sp, percent=75, unin="m", unout="km2") %>% 
  st_as_sf() %>% 
  separate(id, into=c("season", "kdecluster"), remove=FALSE) %>% 
  mutate(kdecluster =as.numeric(kdecluster))
mcp.50 <- mcp(dat.sp, percent=50, unin="m", unout="km2") %>% 
  st_as_sf() %>% 
  separate(id, into=c("season", "kdecluster"), remove=FALSE) %>% 
  mutate(kdecluster =as.numeric(kdecluster))

#6. Get centroids----
kd.cen <- getverticeshr(kd, percent=0.1, unin="m", unout="km2") %>% 
  st_as_sf() %>% 
  st_centroid() %>% 
  st_coordinates() %>% 
  cbind(kd.area)

#7. Buffer----
kd.buff <- kd.cen %>% 
  data.frame() %>% 
  st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
  st_buffer(dist=kd.cen$rad*1000)
  
#8. Plot----
bbs <- read.csv("Data/LBCUBBSClusters.csv") %>% 
  dplyr::filter(nclust==n) %>% 
  mutate(season="breed")

plot.kde <- ggplot() +
  geom_sf(data=kd.95, aes(fill=factor(kdecluster)), alpha=0.4) +
  geom_sf(data=kd.50, aes(fill=factor(kdecluster)), alpha=0.7) +
  geom_point(data=dat, aes(x=X, y=Y, colour=factor(kdecluster))) +
#  geom_point(data=bbs, aes(x=X, y=Y, colour=factor(knncluster))) + 
  geom_point(data=kd.cen, aes(x=X, y=Y), colour="black", size=2) +
  facet_wrap(~season) +
  theme(legend.position="none")

plot.mcp <- ggplot() +
  geom_sf(data=mcp.95, aes(fill=factor(kdecluster)), alpha=0.4) +
#  geom_sf(data=mcp.75, aes(fill=factor(kdecluster)), alpha=0.7) +
    geom_point(data=dat, aes(x=X, y=Y, colour=factor(kdecluster))) +
  #  geom_point(data=bbs, aes(x=X, y=Y, colour=factor(knncluster))) + 
  facet_wrap(~season) +
  theme(legend.position="none")

plot.buff <- ggplot() +
  geom_point(data=dat, aes(x=X, y=Y, colour=factor(kdecluster))) +
#  geom_point(data=bbs, aes(x=X, y=Y, colour=factor(knncluster))) + 
  geom_sf(data=kd.buff, aes(fill=factor(kdecluster)), alpha=0.3) +
#  geom_point(data=kd.cen, aes(x=X, y=Y), colour="black", size=2) +
  facet_wrap(~season) +
  theme(legend.position="none")

grid.arrange(plot.kde, plot.mcp, plot.buff, ncol=3, nrow=1)

#using centroids with a gaussian filter may cause high spatial overlap (but does that matter?)
  
#9. Save out the two options for GEE----
write_sf(kd.area, "Data/LBCURegions_KDE.shp")
write_sf(kd.buff, "Data/LBCURegions_Buffer.shp")
