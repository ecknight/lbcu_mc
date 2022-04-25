library(tidyverse)
library(meanShiftR)
library(maptree)
library(gridExtra)
library(sf)
library(sp)
library(MigConnectivity)
library(raster)
library(vegan)
library(data.table)

options(scipen=9999)


#1. Do eastern populations make more stopovers?----
dat <- read.csv("Data/LBCUMCLocations.csv")

stop.n <- dat %>% 
  dplyr::filter(season %in% c("fallmig", "springmig")) %>% 
  group_by(id, season) %>% 
  summarize(n=n()) %>% 
  ungroup() %>% 
  left_join(dat)

ggplot(stop.n) +
  geom_jitter(aes(x=X, y=n)) +
  geom_smooth(aes(x=X, y=n)) +
  facet_wrap(~season)

#TO DO: TRY WITH STOPOVER####

#1. Import data----
dat <- read.csv("Data/LBCUMCLocations.csv") 

#2. Wrangle data----

#ID birds with breeding and wintering ground locations
dat.n <- dat  %>% 
  dplyr::filter(season %in% c("breed", "winter")) %>% 
  group_by(id, season) %>% 
  summarize(n=n()) %>% 
  group_by(id) %>% 
  summarize(n=n()) %>% 
  dplyr::filter(n==2) %>% 
  ungroup()

#Filter
bw <- dat %>% 
  dplyr::filter(season %in% c("breed", "winter"),
                id %in% dat.n$id) %>% 
  group_by(id, season) %>% 
  summarize(X=mean(X),
            Y=mean(Y)) %>% 
  ungroup() %>% 
  pivot_wider(id_cols=id, names_from=season, values_from=X:Y)

#3. K means clustering----

#3a. Try with both breeding & wintering locations
bw.clust <- bw %>% 
  dplyr::select("X_breed", "Y_breed", "X_winter", "Y_winter")

set.seed(1234)
clust2 <- kmeans(bw.clust, 2)
bw$bwcluster2 <- clust2$cluster
clust3 <- kmeans(bw.clust, 3)
bw$bwcluster3 <- clust3$cluster
clust4 <- kmeans(bw.clust, 4)
bw$bwcluster4 <- clust4$cluster
clust5 <- kmeans(bw.clust, 5)
bw$bwcluster5 <- clust5$cluster
clust6 <- kmeans(bw.clust, 6)
bw$bwcluster6 <- clust6$cluster
clust7 <- kmeans(bw.clust, 7)
bw$bwcluster7 <- clust7$cluster
clust8 <- kmeans(bw.clust, 8)
bw$bwcluster8 <- clust8$cluster
clust9 <- kmeans(bw.clust, 9)
bw$bwcluster9 <- clust9$cluster

breed <- ggplot(bw) +
  geom_point(aes(x=X_breed, y=Y_breed, colour=factor(bwcluster3)), show.legend = FALSE)
winter <- ggplot(bw) +
  geom_point(aes(x=X_winter, y=Y_winter, colour=factor(bwcluster3)), show.legend = FALSE)
grid.arrange(breed, winter, nrow=2)
ggplot(bw) +
  geom_point(aes(x=X_breed, y=Y_breed, colour=factor(bwcluster3)))
ggplot(bw) +
  geom_point(aes(x=X_breed, y=Y_breed, colour=factor(bwcluster4)))
ggplot(bw) +
  geom_point(aes(x=X_breed, y=Y_breed, colour=factor(bwcluster5)))
ggplot(bw) +
  geom_point(aes(x=X_breed, y=Y_breed, colour=factor(bwcluster6)))
ggplot(bw) +
  geom_point(aes(x=X_breed, y=Y_breed, colour=factor(bwcluster7)))
ggplot(bw) +
  geom_point(aes(x=X_breed, y=Y_breed, colour=factor(bwcluster8)))

ggplot(bw) +
  geom_point(aes(x=X_winter, y=Y_winter, colour=factor(bwcluster2)))
ggplot(bw) +
  geom_point(aes(x=X_winter, y=Y_winter, colour=factor(bwcluster3)))
ggplot(bw) +
  geom_point(aes(x=X_winter, y=Y_winter, colour=factor(bwcluster4)))

#3b. Try with just wintering locations
w.clust <- bw %>% 
  dplyr::select("X_winter", "Y_winter")

set.seed(1234)
clust2 <- kmeans(w.clust, 2)
bw$wcluster2 <- clust2$cluster
clust3 <- kmeans(w.clust, 3)
bw$wcluster3 <- clust3$cluster
clust4 <- kmeans(w.clust, 4)
bw$wcluster4 <- clust4$cluster
clust5 <- kmeans(w.clust, 5)
bw$wcluster5 <- clust5$cluster

ggplot(bw) +
  geom_point(aes(x=X_breed, y=Y_breed, colour=factor(wcluster2)))
ggplot(bw) +
  geom_point(aes(x=X_breed, y=Y_breed, colour=factor(wcluster3)))
ggplot(bw) +
  geom_point(aes(x=X_breed, y=Y_breed, colour=factor(wcluster4)))

ggplot(bw) +
  geom_point(aes(x=X_winter, y=Y_winter, colour=factor(wcluster2)))
ggplot(bw) +
  geom_point(aes(x=X_winter, y=Y_winter, colour=factor(wcluster3)))
ggplot(bw) +
  geom_point(aes(x=X_winter, y=Y_winter, colour=factor(wcluster4)))

#4. Meanshift----

mat1 <- matrix(bw$X_winter)
mat2 <- matrix(bw$Y_winter)
mat <- cbind(mat1, mat2)

shift <- meanShift(mat,
                   algorithm="KDTREE",
                   bandwidth=c(1000000,1000000))

bw$shift <- shift[[1]]

ggplot(bw) +
  geom_point(aes(x=X_winter, y=Y_winter, colour=factor(shift)))

#5. Hierarchical clustering----
bw.clust <- bw %>% 
  dplyr::select("X_breed", "Y_breed", "X_winter", "Y_winter")

bw.dist <- dist(bw.clust)

bw.hclust <- hclust(bw.dist, method="complete")
bw.kgs <- kgs(bw.hclust, bw.dist)

bw$hclust <- cutree(bw.hclust, 9)

breed <- ggplot(bw) +
  geom_point(aes(x=X_breed, y=Y_breed, colour=factor(hclust)), show.legend = FALSE)
winter <- ggplot(bw) +
  geom_point(aes(x=X_winter, y=Y_winter, colour=factor(hclust)), show.legend = FALSE)

grid.arrange(breed, winter, nrow=2)

#6. Estimate connectivity between pairs----
methods <- c("bwcluster2", "bwcluster3", "bwcluster4", "bwcluster5", "bwcluster6", "bwcluster7", "bwcluster8", "bwcluster9")

i.list <- list()
for(i in 1:length(methods)){
  
  clusters <- bw %>% 
    dplyr::select(id, X_breed, X_winter, Y_breed, Y_winter, methods[i]) %>% 
    rename(method = methods[i])
  nclust <- max(clusters$method)
  pairs <- expand.grid(1:nclust, 1:nclust) %>% 
    dplyr::filter(Var1!=Var2)
  
  j.list <- list()
  for(j in 1:nrow(pairs)){
    
    pairs.j <- c(pairs$Var1[j], pairs$Var2[j])
    
    bw.j <- clusters %>% 
      dplyr::filter(method %in% pairs.j)
    
    breed.mat.j <- bw.j %>% 
      dplyr::select(X_breed, Y_breed) %>% 
      vegdist("euclidean")
    
    winter.mat.j <- bw.j %>% 
      dplyr::select(X_winter, Y_winter) %>% 
      vegdist("euclidean")
    
    mantel.j <- mantel(breed.mat.j, winter.mat.j)
    
    j.list[[j]] <- data.frame(r = mantel.j[["statistic"]],
                              p = mantel.j[["signif"]],
                              n=nrow(bw.j),
                              cluster1 = pairs$Var1[j],
                              cluster2 = pairs$Var2[j],
                              method=methods[i])
    
  }
  
  i.list[[i]] <- rbindlist(j.list)
  
}

m.list <- rbindlist(i.list, fill=TRUE)

m.filter <- m.list %>% 
  dplyr::select(r, p, n, method) %>% 
  unique()

m.summary <- m.filter %>% 
  dplyr::filter(!is.na(r)) %>% 
  group_by(method) %>% 
  summarize(mean = mean(r),
            sd = sd(r),
            max = max(r),
            min = min(r)) %>% 
  ungroup()
m.summary

#7. Estimate connectivity within clusters----

methods <- c("bwcluster2", "bwcluster3", "bwcluster4", "bwcluster5", "bwcluster6", "bwcluster7", "bwcluster8")

i.list <- list()
for(i in 1:length(methods)){
  
  clusters <- bw %>% 
    dplyr::select(id, X_breed, X_winter, Y_breed, Y_winter, methods[i]) %>% 
    rename(method = methods[i])
  nclust <- max(clusters$method)
  
  j.list <- list()
  for(j in 1:nclust){
    
    bw.j <- clusters %>% 
      dplyr::filter(method==j)
    
    breed.mat.j <- bw.j %>% 
      dplyr::select(X_breed, Y_breed) %>% 
      vegdist("euclidean")
    
    winter.mat.j <- bw.j %>% 
      dplyr::select(X_winter, Y_winter) %>% 
      vegdist("euclidean")
    
    mantel.j <- mantel(breed.mat.j, winter.mat.j)
    
    j.list[[j]] <- data.frame(r = mantel.j[["statistic"]],
                              p = mantel.j[["signif"]],
                              n=nrow(bw.j),
                              cluster = j,
                              method=methods[i])
    
  }
  
  i.list[[i]] <- rbindlist(j.list)
  
}

m.list <- rbindlist(i.list)

m.summary <- m.list %>% 
  dplyr::filter(!is.na(r)) %>% 
  group_by(method) %>% 
  summarize(mean = mean(r),
            sd = sd(r),
            max = max(r),
            min = min(r)) %>% 
  ungroup()
m.summary

#7. Plot----
breed <- ggplot(bw) +
  geom_point(aes(x=X_breed, y=Y_breed, colour=factor(hclust)))
winter <- ggplot(bw) +
  geom_point(aes(x=X_winter, y=Y_winter, colour=factor(hclust)))

hist <- ggplot(m.filter) +
  geom_histogram(aes(x=r)) +
  facet_wrap(~method)

matrix <- ggplot(m.list) +
  geom_raster(aes(x=cluster1, y=cluster2, fill=r)) +
  scale_fill_viridis_c()

clustering <- grid.arrange(breed, hist, winter, matrix, nrow=2, ncol=2)
ggsave(clustering, filename = "Figs/Clustering.jpg", width = 12, height=8)


#8. Look at BBS distribution----
load("~/Library/Application Support/bbsBayes/bbs_raw_data.RData")
dat <- bbs_data[["bird"]] %>% 
  dplyr::filter(AOU==2640)
routes <- dat %>% 
  dplyr::select(RouteDataID) %>%
  unique() %>% 
  left_join(bbs_data[["route"]]) %>% 
  group_by(countrynum, statenum, Route, Latitude, Longitude) %>%
  summarize(n = n()) %>% 
  dplyr::filter(!is.na(Latitude)) %>% 
  ungroup()

dat.sf <- routes %>% 
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326) %>% 
  st_transform(crs=3857)

plot(dat.sf)

dat.sp <- as_Spatial(dat.sf)

ggplot(bw) +
  geom_point(data=data.frame(st_coordinates(dat.sf)), aes(x=X, y=Y), colour="red") +
  geom_point(aes(x=X_breed, y=Y_breed))


#9. Spatial overlap----
library(adehabitatHR)
library(sp)

trends <- read.csv("bbsBayesModels/Trends_gamye.csv")

bw.sp <- SpatialPointsDataFrame(coords=cbind(bw$X_breed, bw$Y_breed), 
                                data=data.frame(ID=bw$bwcluster9),
                                proj4string = CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"))

mp <- mcp(bw.sp[,1], percent=100)

plot(mp)
#plot(dat.sp, col=dat.sp$n, pch=2, add=TRUE)
plot(bw.sp, col=bw.sp$ID, add=TRUE)

#10. DTW----
library(dtwclust)

#Migrations only
all.df <- read.csv("Data/LBCU_FilteredData_Segmented.csv") %>% 
  dplyr::filter(!id %in% c(46768277, 33088),
                segment!="stationary") %>% 
  mutate(season=case_when(season=="breed" & segment=="arrive" ~ "springmig",
                          season=="breed" & segment=="depart" ~ "fallmig",
                          season=="winter" & segment=="arrive" ~ "fallmig", 
                          season=="winter" & segment=="depart" ~ "springmig",
                          !is.na(season) ~ season),
         legid = paste(id, year, season, sep="-"))

all.df.s <- all.df %>% 
  dplyr::filter(season=="springmig") %>% 
  arrange(legid) %>% 
  dplyr::select(legid, X, Y)

all.df.f <- all.df %>% 
  dplyr::filter(season=="fallmig") %>% 
  arrange(legid) %>% 
  dplyr::select(legid, X, Y)

all.s <- group_split(all.df.s, legid, .keep=FALSE)
names(all.s) <- unique(all.df.s$legid)

all.f <- group_split(all.df.f, legid, .keep=FALSE)
names(all.f) <- unique(all.df.f$legid)

#By day of year
all.df <- read.csv("Data/LBCU_FilteredData_Segmented.csv") %>% 
  dplyr::filter(!id %in% c(46768277, 33088)) %>% 
  mutate(legid = paste(id, year, sep="-"))

all.df.s <- all.df %>% 
  dplyr::filter(doy > 35 & doy < 220) %>% 
  arrange(legid) %>% 
  dplyr::select(legid, X, Y)

all.s <- group_split(all.df.s, legid, .keep=FALSE)
names(all.s) <- unique(all.df.s$legid)

clust.s <- tsclust(all.s, k = 3L:9L,
                   distance = "dtw_basic", centroid = "pam",
                   seed = 94L)
names(clust.s) <- paste0("k_", 3L:9L)
sapply(clust.s, cvi, type = "internal")

clust.s.4 <- tsclust(all.s, k = 4L,
                     distance = "dtw_basic", centroid = "pam",
                     seed = 94L)
plot(clust.s.4)

all.clust.s <- data.frame(legid = unique(all.df.s$legid),
                          clust = clust.s.4@cluster) %>% 
  full_join(all.df %>% 
              dplyr::filter(doy > 35 & doy < 220))

ggplot(all.clust.s) +
  geom_line(aes(x=X, y=Y, group=legid, colour=factor(clust)))

ggplot(all.clust.s %>% dplyr::filter(season=="breed")) +
  geom_point(aes(x=X, y=Y, colour=factor(clust)))

ggplot(all.clust.s %>% dplyr::filter(season=="winter")) +
  geom_point(aes(x=X, y=Y, colour=factor(clust)))


#11. K-means with stopovers----
#How to deal with multiple years of data??? Randomly select one year and then bootstrap?
yr <- dat %>% 
  dplyr::filter(id %in% dat.n$id) %>% 
  dplyr::select(id, year) %>% 
  group_by(id) %>% 
  sample_n(1) %>% 
  ungroup()

#Take out cluster #4, there's only 1 bird season with
bmw.1 <- inner_join(dat, yr) %>% 
  dplyr::filter(cluster < 4) %>% 
  mutate(clusterid = paste0(season, "_", cluster))  %>% 
  dplyr::select(id, clusterid, X, Y)

#Calculate means for NAs
nas <- inner_join(dat, yr) %>% 
  dplyr::filter(cluster < 4) %>% 
  mutate(clusterid = paste0(season, "_", cluster)) %>% 
  group_by(clusterid) %>% 
  summarize(X=mean(X),
            Y=mean(Y)) %>% 
  ungroup()

bmw.0 <- expand.grid(clusterid=nas$clusterid, id=yr$id) %>% 
  anti_join(bmw.1) %>% 
  left_join(nas)

bmw <- rbind(bmw.1, bmw.0) %>% 
  dplyr::select(id, clusterid, X, Y) %>% 
  pivot_wider(id_cols=id, names_from=clusterid, values_from=X:Y) 

#3a. Try with both breeding & wintering locations
bmw.clust <- bmw %>% 
  dplyr::select(-id)

set.seed(1234)
clust3 <- kmeans(bmw.clust, 3)
bmw$cluster3 <- clust3$cluster
clust4 <- kmeans(bmw.clust, 4)
bmw$cluster4 <- clust4$cluster
clust8 <- kmeans(bmw.clust, 8)
bmw$cluster8 <- clust8$cluster

ggplot(bmw) +
  geom_point(aes(x=X_breed_1, y=Y_breed_1, colour=factor(cluster4)))
ggplot(bmw) +
  geom_point(aes(x=X_winter_1, y=Y_winter_1, colour=factor(cluster4)), show.legend = FALSE)







#12. Trying to get MC to work----
library(MigConnectivity)
