library(tidyverse)
library(ClustImpute)
library(data.table)

options(scipen=9999)

#1. Import data----
dat.raw <- read.csv("Data/LBCUMCLocations.csv") 

#NOTE: Bootstrapping could occur in step 2 or step 3. Should probably do step 2

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

#4. Pick one year of data for each individual and make it wide----
set.seed(1234)
dat.i <- dat.days %>% 
  dplyr::filter(id %in% dat.n$id) %>% 
  dplyr::select(id, year) %>% 
  unique() %>% 
  group_by(id) %>% 
  sample_n(1) %>% 
  ungroup() %>% 
  left_join(dat)  %>% 
  dplyr::select(id, year, season, X, Y) %>% 
  pivot_wider(id_cols=id:year, names_from=season, values_from=X:Y) %>% 
  data.frame()

#5. KDE with incomplete data clustering----
clusters <- c(2:10)

kde.cluster <- list()
for(j in 1:length(clusters)){
  
  dat.j <- dat.i %>% 
    dplyr::select(-id, -year)
  kde.j <- ClustImpute(dat.j, clusters[j], nr_iter=100)
  kde.cluster[[j]] <- data.frame(clusters = kde.j$clusters,
                                 nclust = clusters[j],
                                 id = dat.i$id,
                                 year = dat.i$year)
  
}

dat.kde <- rbindlist(kde.cluster) %>% 
  pivot_wider(id_cols=id:year, names_from=nclust, values_from=clusters, names_prefix="kde_") %>% 
  left_join(dat) %>% 
  pivot_longer(names_to="nclust", values_to="kdecluster", cols=kde_2:kde_10, names_prefix="kde_") %>% 
  mutate(nclust = as.numeric(nclust))

#6. Visualize----
ggplot(dat.kde) +
  geom_point(aes(x=X, y=Y, colour=factor(kdecluster))) +
  facet_grid(season ~ nclust)

ggsave(filename="Figs/KDE_stopovers.jpeg", width=18, height = 10)

#7. Save clusters----
write.csv(dat.kde, "Data/LBCUKDEClusters.csv", row.names=FALSE)

