library(tidyverse)
library(ClustImpute)
library(data.table)

options(scipen=9999)

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

dat <- dat.days %>% 
#  dplyr::filter(season %in% c("breed", "winter")) %>% 
  dplyr::filter(cluster.days==1)

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

#4. Set # of bootstraps & set up loop----
boot <- 100

dat.kde <- list()
set.seed(1)
for(i in 1:boot){
  
  #5. Pick one point for each season for each individual and make it wide----
  
  dat.i <- dat %>% 
    dplyr::filter(id %in% dat.n$id) %>% 
    group_by(id, season) %>% 
    sample_n(1) %>% 
    ungroup() %>% 
    dplyr::select(id, season, X, Y) %>% 
    pivot_wider(id_cols=id, names_from=season, values_from=X:Y) %>% 
    data.frame()
  
  #6. KDE with incomplete data clustering----
  clusters <- c(2:9)
  
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
    pivot_longer(names_to="nclust", values_to="kdecluster", cols=kde_2:kde_9, names_prefix="kde_") %>% 
    mutate(nclust = as.numeric(nclust),
           boot = i)
  
  #7. Visualize----
  # ggplot(dat.kde) +
  #   geom_point(aes(x=X, y=Y, colour=factor(kdecluster))) +
  #   facet_grid(season ~ nclust)
  # 
  # ggsave(filename="Figs/KDE_stopovers.jpeg", width=18, height = 10)
  
print(paste0("Finished bootstrap ", i, " of ", boot))
  
}

#8. Save----
dat.out <- rbindlist(dat.kde) %>% 
  pivot_longer(cols=c(X_breed:Y_winter),
               names_to="season", values_to="value") %>% 
  separate(season, into=c("coord", "season")) %>% 
  pivot_wider(names_from=coord, values_from=value)
  
write.csv(dat.out, "Data/LBCUKDEClusters.csv", row.names=FALSE)

