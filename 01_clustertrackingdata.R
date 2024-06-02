library(tidyverse)
library(ClustImpute)
library(data.table)
library(missMDA)
library(FactoMineR)
library(sf)

options(scipen=9999)

#TO DO: CLUSTER SPATIALLY AND TAKE PROMINANT POINT ACROSS DAYS#####

#1. Import data----
raw <- read.csv("Data/LBCUMCLocations.csv")

dat <- raw %>% 
  dplyr::filter(!id %in% c(46768277, 33088, 129945787, 46770723, 46769927, 86872)) %>% 
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

dat.use <- dat |> 
  dplyr::filter(id %in% dat.n$id)

#3. Get single seasonal location for each individual----
#This is currently just the single location with most days, could cluster spatially and sum days across years
dat.main <- dat.use |> 
  group_by(id, season) |> 
  dplyr::filter(days==max(days)) |> 
  sample_n(1) |> 
  ungroup()

#4. Make wide----
dat.wide <- dat.main |> 
  select(id, season, X, Y) |> 
  pivot_wider(id_cols=id, names_from=season, values_from=X:Y) %>%
  data.frame()

dat.clust <- dat.wide |> 
  dplyr::select(-id)

#5. KDE with incomplete data clustering----

clusters <- c(2:5)

kde.cluster <- list()
for(j in 1:length(clusters)){
  kde.j <- ClustImpute(dat.clust, clusters[j], nr_iter=100)
  kde.cluster[[j]] <- data.frame(group = kde.j$clusters,
                                 nclust = clusters[j],
                                 id = dat.wide$id)
}

dat.out <- rbindlist(kde.cluster) %>% 
  pivot_wider(id_cols=id, names_from=nclust, values_from=group, names_prefix="kde_") %>% 
  left_join(dat.main) |> 
  pivot_longer(names_to="nclust", values_to="group", cols=kde_2:kde_5, names_prefix="kde_")

#6. Check nclust have at least 5 individuals----
table(dat.out$nclust, dat.out$group)

#7. Add expert clusters----
dat.expert <- dat.main  %>% 
  dplyr::filter(season=="winter") %>% 
  mutate(group = case_when(X > -10960000 & distance < 100000 ~ 4,
                           lon < -105 & lon > -108 & distance < 10000 ~ 2,
                           lon < -108 & lon > -118 ~ 2,
                           lon < -118 ~ 1,
                           !is.na(lon) ~ 3)) %>% 
  mutate(nclust="expert") |> 
  dplyr::select(id, group, nclust) |> 
  left_join(dat.main)

#8. Add flyway clusters----

#Get shapefile
sf_use_s2(FALSE)
flyway <- read_sf("Data/Atlas Regions/Final_globalregions.shp") |> 
  dplyr::filter(region %in% c("Pacific", "Midcontinent", "Atlantic")) |> 
  st_make_valid()

#Intersect
dat.fly <- dat.expert |> 
  dplyr::filter(season!="winter") |> 
  st_as_sf(coords=c("lon", "lat"), crs=4326, remove=FALSE) |> 
  st_intersection(flyway) |> 
  st_drop_geometry()

#Look at individuals with more than one flyway
dat.multi <- dat.fly |> 
  dplyr::select(id, region) |> 
  unique() |> 
  group_by(id) |> 
  summarize(n=n()) |> 
  ungroup() |> 
  dplyr::filter(n > 1) |> 
  left_join(dat.fly)

ggplot() +
  geom_sf(data=flyway, aes(group=region), fill=NA) +
  geom_line(data=dat.multi, aes(x=lon, y=lat, group=id)) + 
  geom_point(data=dat.multi, aes(x=lon, y=lat, colour=season), size=2) +
  xlim(c(-120, -100)) +
  ylim(c(30, 50))

#Pick the region of migration 

#9. Put together----
dat.final <- rbind(dat.out, dat.expert, dat.fly)

#10. Plot----
ggplot(dat.final) +
  geom_point(aes(x=X, y=Y, colour=factor(group))) +
  facet_grid(season ~ nclust, scales="free_y") + 
  scale_colour_viridis_d()

ggsave(filename="Figs/KDE_stopovers.jpeg", width=18, height = 10)

#11. Save----
write.csv(dat.final, "Data/LBCUKDEClusters.csv", row.names=FALSE)
