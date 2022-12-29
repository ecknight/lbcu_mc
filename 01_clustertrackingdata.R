library(tidyverse)
library(ClustImpute)
library(data.table)

options(scipen=9999)

#1. Import data----
dat <- read.csv("Data/LBCUMCLocations.csv") %>% 
  rename(seasoncluster=cluster)

#2. Use only birds with known breeding & wintering location---
#ID birds with breeding and wintering ground locations
dat.n <- dat  %>% 
  dplyr::filter(season %in% c("breed", "winter")) %>% 
  group_by(id, season) %>% 
  summarize(n=n()) %>% 
  group_by(id) %>% 
  summarize(n=n()) %>% 
  dplyr::filter(n==2) %>% 
  ungroup()

#3. Set # of bootstraps & set up loop----
boot <- 100

dat.kde <- list()
set.seed(1234)
for(i in 1:boot){
  
  #4. Pick one point for each season for each individual & make it wide----
  dat.i <- dat %>% 
    dplyr::filter(id %in% dat.n$id) %>% 
    group_by(id, season) %>% 
    sample_n(1) %>% 
    ungroup()
  
  #Select columns for clustering
  dat.j <- dat.i %>% 
    dplyr::select(id, season, X, Y) %>% 
    pivot_wider(id_cols=id, names_from=season, values_from=X:Y) %>%
    data.frame()
  
  dat.k <- dat.j %>% 
    dplyr::select(-id)
  
  #5. KDE with incomplete data clustering----
  clusters <- c(2:6)
  
  kde.cluster <- list()
  for(j in 1:length(clusters)){
    kde.j <- ClustImpute(dat.k, clusters[j], nr_iter=100)
    kde.cluster[[j]] <- data.frame(group = kde.j$clusters,
                                   nclust = clusters[j],
                                   id = dat.j$id)
  }
  
  #6. Add manual clusters----
  dat.man <- dat.i  %>% 
    dplyr::filter(season=="winter") %>% 
    mutate(group = case_when(X > -10960000 & distance < 100000 ~ 4,
                             lon < -105 & lon > -108 & distance < 10000 ~ 2,
                             lon < -108 & lon > -118 ~ 2,
                             lon < -118 ~ 1,
                             !is.na(lon) ~ 3)) %>% 
    dplyr::select(id, group) %>% 
    unique() %>% 
    full_join(dat.i) %>% 
    data.frame() %>% 
    mutate(nclust="manual")
  
  #7. Put together----
  dat.kde[[i]] <- rbindlist(kde.cluster) %>% 
    pivot_wider(id_cols=id, names_from=nclust, values_from=group, names_prefix="kde_") %>% 
    left_join(dat.i) %>% 
    pivot_longer(names_to="nclust", values_to="group", cols=kde_2:kde_6, names_prefix="kde_") %>% 
    rbind(dat.man) %>% 
    mutate(boot = i)
  
print(paste0("Finished bootstrap ", i, " of ", boot))
  
}

#8. Wrangle output----
dat.out <- rbindlist(dat.kde)

#9. Look at variation----
dat.sum <- dat.out %>% 
  group_by(id, season, nclust, group) %>% 
  summarize(n=n()) %>% 
  ungroup()

ggplot(dat.sum) +
  geom_histogram(aes(x=n)) +
  facet_grid(season~nclust)

#10. Remove nclust with < 5 individuals in a cluster----
# dat.n <- dat.out %>% 
#   dplyr::filter(season=="winter") %>% 
#   group_by(nclust, group) %>% 
#   summarize(n=n()) %>% 
#   dplyr::filter(n/100 < 5)
# 
# dat.final <- dat.out %>% 
#   dplyr::filter(!nclust %in% dat.n$nclust)

#11. Plot----
ggplot(dat.out) +
  geom_point(aes(x=X, y=Y, colour=factor(group))) +
  facet_grid(season ~ nclust, scales="free_y") + 
  scale_colour_viridis_d()

ggsave(filename="Figs/KDE_stopovers.jpeg", width=18, height = 10)

#12. Save----
write.csv(dat.out, "Data/LBCUKDEClusters.csv", row.names=FALSE)
