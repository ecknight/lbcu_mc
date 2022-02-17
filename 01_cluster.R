library(tidyverse)
library(meanShiftR)
library(maptree)
library(gridExtra)

options(scipen=9999)

#TO DO: TRY WITH STOPOVER####

#1. Import data----
dat <- read.csv("Data/LBCUMCLocations.csv") 

#2. Wrangle data----

#ID birds with breeding and wintering ground lcoations
dat.n <- dat  %>% 
  dplyr::filter(season %in% c("breed", "winter")) %>% 
  group_by(id) %>% 
  summarize(n=n()) %>% 
  dplyr::filter(n==2) %>% 
  ungroup()

#Filter
bw <- dat %>% 
  dplyr::filter(season %in% c("breed", "winter"),
                id %in% dat.n$id) %>% 
  dplyr::select(id, season, X, Y) %>% 
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

breed <- ggplot(bw) +
  geom_point(aes(x=X_breed, y=Y_breed, colour=factor(bwcluster8)), show.legend = FALSE)
winter <- ggplot(bw) +
  geom_point(aes(x=X_winter, y=Y_winter, colour=factor(bwcluster8)), show.legend = FALSE)
grid.arrange(breed, winter, nrow=2)

ggplot(bw) +
  geom_point(aes(x=X_breed, y=Y_breed, colour=factor(bwcluster3)))
ggplot(bw) +
  geom_point(aes(x=X_breed, y=Y_breed, colour=factor(bwcluster4)))

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

bw$hclust <- cutree(bw.hclust, 8)

breed <- ggplot(bw) +
  geom_point(aes(x=X_breed, y=Y_breed, colour=factor(hclust)), show.legend = FALSE)
winter <- ggplot(bw) +
  geom_point(aes(x=X_winter, y=Y_winter, colour=factor(hclust)), show.legend = FALSE)

grid.arrange(breed, winter, nrow=2)