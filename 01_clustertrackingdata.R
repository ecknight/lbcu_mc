library(tidyverse)
library(ClustImpute)
library(data.table)

options(scipen=9999)

#TO DO: SCALE FROM 0 TO 1

#1. Import data----
dat.raw <- read.csv("Data/LBCUMCLocations.csv")

#2. Wrangle data to longest duration winter & stopover location for each individual----
dat.days <- dat.raw %>% 
  group_by(id, year, season) %>% 
  arrange(-days) %>% 
  mutate(cluster.days = row_number()) %>% 
  ungroup()

dat.days %>% 
  group_by(season, cluster.days) %>% 
  summarize(mean = mean(days), 
            sd = sd(days),
            max = max(days),
            min = min(days))

ggplot(dat.days) +
  geom_histogram(aes(x=days)) +
  facet_grid(cluster.days ~ season)

# dat <- dat.days %>% 
# #  dplyr::filter(season %in% c("breed", "winter")) %>% 
#   dplyr::filter(cluster.days==1)

dat <- dat.days

table(dat$season)

#3. Use only birds with known breeding & wintering location---
#ID birds with breeding and wintering ground locations
dat.n <- dat  %>% 
  dplyr::filter(season %in% c("breed", "winter")) %>% 
  group_by(id, season) %>% 
  summarize(n=n()) %>% 
  group_by(id) %>% 
  summarize(n=n()) %>% 
  dplyr::filter(n==2) %>% 
  ungroup()

#4. Scale variables from 0 to 1 for each season----
dat.scale <- dat %>% 
  dplyr::filter(id %in% dat.n$id) %>% 
  group_by(season) %>% 
  mutate(Xs = (X-min(X))/(max(X)-min(X)),
         Ys = (Y-min(Y))/(max(Y)-min(Y)),
         elevs = (elevation-min(elevation))/(max(elevation)-min(elevation)),
         dists = (distance-min(distance))/(max(distance)-min(distance))) %>% 
  ungroup()

#4. Set # of bootstraps & set up loop----
boot <- 10

dat.kde <- list()
set.seed(1234)
for(i in 1:boot){
  
  #5. Pick one point for each season for each individual & make it wide----
  dat.i <- dat.scale %>% 
    dplyr::filter(season %in% c("breed", "winter")) %>% 
    group_by(id, season) %>% 
    sample_n(1) %>% 
#    dplyr::select(id, season, Xs, Ys, dists) %>% 
    dplyr::select(id, season, Ys, dists) %>% 
    pivot_wider(id_cols=id, names_from=season, values_from=Ys:dists) %>%
    data.frame()
  
  #8. KDE with incomplete data clustering----
  clusters <- c(2:6)
  
  kde.cluster <- list()
  for(j in 1:length(clusters)){
    
    dat.j <- dat.i %>% 
      dplyr::select(-id)
    kde.j <- ClustImpute(dat.j, clusters[j], nr_iter=100)
    kde.cluster[[j]] <- data.frame(clusters = kde.j$clusters,
                                   nclust = clusters[j],
                                   id = dat.i$id)
    
  }
  
  dat.kde[[i]] <- rbindlist(kde.cluster) %>% 
    pivot_wider(id_cols=id, names_from=nclust, values_from=clusters, names_prefix="kde_") %>% 
    left_join(dat.i) %>% 
    pivot_longer(names_to="nclust", values_to="kdecluster", cols=kde_2:kde_6, names_prefix="kde_") %>% 
    mutate(nclust = as.numeric(nclust),
           boot = i)
  
print(paste0("Finished bootstrap ", i, " of ", boot))
  
}

#8. Wrangle output----
dat.out <- rbindlist(dat.kde) %>% 
  # pivot_longer(cols=c(Xs_breed:dists_winter),
  #              names_to="season", values_to="value") %>% 
  pivot_longer(cols=c(Ys_breed:dists_winter),
               names_to="season", values_to="value") %>% 
  separate(season, into=c("coord", "season")) %>% 
  pivot_wider(names_from=coord, values_from=value)

#9. Look at variation----
dat.sum <- dat.out %>% 
  group_by(id, season, nclust, kdecluster) %>% 
  summarize(n=n()) %>% 
  ungroup()

ggplot(dat.sum) +
  geom_histogram(aes(x=n)) +
  facet_grid(season~nclust)

#10. Get main cluster id----
dat.mean <- dat.out %>% 
#  dplyr::filter(season=="breed") %>% 
  group_by(id, season, nclust) %>% 
  # summarize(X = mean(Xs, na.rm=FALSE), 
  #           Y = mean(Ys, na.rm=FALSE),
  #           elevation = mean(elevs, na.rm=FALSE),
  #           distance = mean(dists, na.rm=FALSE),
  #           kdecluster = round(mean(kdecluster))) %>% 
  summarize(Ys = mean(Ys, na.rm=FALSE),
            dists = mean(dists, na.rm=FALSE),
            kdecluster = round(mean(kdecluster))) %>% 
  ungroup() %>% 
  left_join(dat.scale)

#11. Look at # of individuals per cluster----
dat.n <- dat.mean %>% 
  dplyr::filter(season=="breed") %>% 
  group_by(nclust, kdecluster) %>% 
  summarize(n=n()) %>% 
  dplyr::filter(n < 5)

#12. Remove nclust with < 5 individuals in a cluster----
dat.final <- dat.out %>% 
  dplyr::filter(!nclust %in% dat.n$nclust)

#13. Plot----
ggplot(dat.mean %>% dplyr::filter(!nclust %in% dat.n$nclust)) +
  geom_point(aes(x=X, y=Y, colour=factor(kdecluster), size=dists)) +
  facet_grid(season ~ nclust, scales="free_y")

ggsave(filename="Figs/KDE_stopovers.jpeg", width=18, height = 10)

#14. Save----
write.csv(dat.final, "Data/LBCUKDEClusters.csv", row.names=FALSE)
#dat.out <- read.csv("Data/LBCUKDEClusters.csv")


