library(tidyverse)
library(ClustImpute)
library(data.table)

options(scipen=9999)

#1. Import data----
dat <- read.csv("Data/LBCU_FilteredData_Segmented.csv") %>% 
  dplyr::filter(!id %in% c(46768277, 33088, 129945787, 46770723, 46769927))

#2. Use only birds with known breeding & wintering location---
#ID birds with breeding and wintering ground locations
birds <- dat  %>% 
  dplyr::filter(season %in% c("breed", "winter")) %>% 
  group_by(id, season) %>% 
  summarize(n=n()) %>% 
  group_by(id) %>% 
  summarize(n=n()) %>% 
  dplyr::filter(n==2) %>% 
  ungroup()

dat.birds <- dat %>% 
  dplyr::filter(id %in% birds$id)

#3. Set # of bootstraps & set up loop----
boot <- 100

dat.kde <- list()
set.seed(1)
for(i in 1:boot){
  
  #4. Pick day of year to start at----
  day.i <- round(runif(1, 1, 365))
  
  #5. Pick one year of data for each bird----
  dat.i <- dat.birds %>% 
    arrange(id, date) %>% 
    group_by(id) %>% 
    mutate(day = row_number()) %>% 
    dplyr::filter(doy==day.i) %>% 
    sample_n(1) %>% 
    ungroup() %>% 
    dplyr::select(id, day) %>% 
    rename(start = day) %>% 
    mutate(end = start + 364) %>% 
    right_join(dat.birds) %>% 
    arrange(id, date) %>% 
    group_by(id) %>% 
    mutate(day = row_number()) %>% 
    dplyr::filter(day >= start & day <= end) %>% 
    mutate(day = row_number()) %>% 
    ungroup() %>% 
    dplyr::select(id, day, X, Y) %>% 
    pivot_wider(id_cols=id, names_from=day, values_from=X:Y) %>% 
    data.frame()
  
  #6. KDE with incomplete data clustering----
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
  
  #7. Save loop output----
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
  pivot_longer(cols=c(X_1:Y_365),
               names_to="day", values_to="value") %>% 
  separate(day, into=c("coord", "day")) %>% 
  pivot_wider(names_from=coord, values_from=value)

write.csv(dat.out, "Data/LBCUKDEClusters_year.csv", row.names=FALSE)

#9. Look at group membership----
locs <- read.csv("Data/LBCUMCLocations.csv")

dat.clust <- dat.out %>% 
  dplyr::select(id, nclust, kdecluster, boot) %>% 
  unique() %>% 
  left_join(locs)

ggplot(dat.clust) +
  geom_point(aes(x=X, y=Y, colour=factor(kdecluster))) +
  facet_grid(season~nclust) +
  scale_colour_viridis_d()

#10. Look at variation----
dat.sum <- dat.out %>% 
  group_by(id, nclust, kdecluster) %>% 
  summarize(n=n()) %>% 
  ungroup()

ggplot(dat.sum) +
  geom_histogram(aes(x=n)) +
  facet_grid(~nclust)