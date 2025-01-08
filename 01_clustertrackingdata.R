library(tidyverse)
library(ClustImpute)
library(data.table)
library(downloader)
library(sf)

#1. Import data----
dat <- read.csv("Data/LBCUMCLocations.csv")

#2. Use only birds with known breeding & wintering location---
#ID birds with breeding and wintering ground locations
dat.n <- dat  |> 
  dplyr::filter(season %in% c("breed", "winter")) |> 
  group_by(id, season) |> 
  summarize(n=n()) |> 
  group_by(id) |> 
  summarize(n=n()) |> 
  dplyr::filter(n==2) |> 
  ungroup()

dat.use <- dat |> 
  dplyr::filter(id %in% dat.n$id)

#3. Get single seasonal location for each individual----
#This is the single location with most days
set.seed(1234)
dat.main <- dat.use |>
  group_by(id, season) |>
  dplyr::filter(days==max(days)) |>
  sample_n(1) |>
  ungroup()

#4. Make wide----
dat.wide <- dat.main |> 
  select(id, season, X, Y) |> 
  pivot_wider(id_cols=id, names_from=season, values_from=X:Y) |>
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

dat.out <- rbindlist(kde.cluster) |> 
  pivot_wider(id_cols=id, names_from=nclust, values_from=group, names_prefix="kde_") |> 
  left_join(dat.main, multiple="all") |> 
  pivot_longer(names_to="nclust", values_to="group", cols=kde_2:kde_5, names_prefix="kde_")

#6. Check nclust have at least 5 individuals----
table(dat.out$nclust, dat.out$group)

#7. Add expert clusters----

#Download coastal shapefile - only do this once
temp <- tempfile()
download("https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/gshhg-shp-2.3.7.zip", temp)
unzip(zipfile=temp, exdir="Data")
unlink(temp)

#Get distance to coast
coast <- read_sf("Data/GSHHS_shp/l/GSHHS_l_L1.shp") |> 
  st_make_valid() |> 
  st_cast("LINESTRING")

dat.sf <- dat.main |> 
  unique() |> 
  st_as_sf(coords=c("lon", "lat"), crs=4326, remove=FALSE) 

dat.near <- dat.sf |> 
  st_nearest_feature(coast)

dat.coast <- data.frame(distance = as.numeric(st_distance(dat.sf, coast[dat.near,], by_element = TRUE))) |> 
  cbind(dat.main)

#check
ggplot(dat.coast) +
  geom_point(aes(x=lon, y=lat, colour=distance))

#assign cluster
dat.expert <- dat.coast  |> 
  dplyr::filter(season=="winter") |> 
  mutate(group = case_when(X > -10960000 & distance < 100000 ~ 4,
                           lon < -105 & lon > -108 & distance < 10000 ~ 2,
                           lon < -108 & lon > -118 ~ 2,
                           lon < -118 ~ 1,
                           !is.na(lon) ~ 3)) |> 
  mutate(nclust="expert") |> 
  dplyr::select(id, group, nclust) |> 
  left_join(dat.main, multiple = "all")

#8. Add flyway clusters----

#Get shapefile - available as part of code and supplementary data package on Zenodo 10.5281/zenodo.14607317
sf_use_s2(FALSE)
flyway <- read_sf("Data/Atlas Regions/Final_globalregions.shp") |> 
  dplyr::filter(region %in% c("Pacific", "Midcontinent", "Atlantic")) |> 
  st_make_valid()

#Intersect using breeding location
dat.fly <- dat.main |> 
  dplyr::filter(season=="breed") |> 
  st_as_sf(coords=c("lon", "lat"), crs=4326, remove=FALSE) |> 
  st_intersection(flyway) |> 
  st_drop_geometry() |> 
  dplyr::select(id, region) |> 
  unique() |> 
  left_join(dat.main, multiple="all") |> 
  mutate(group = ifelse(region=="Pacific", 1, 2),
         nclust = "flyway") |> 
  dplyr::select(-region)

#9. Put together----
dat.final <- rbind(dat.out, dat.expert, dat.fly)

#10. Visualize----
ggplot(dat.final) +
  geom_point(aes(x=X, y=Y, colour=factor(group))) +
  facet_grid(season ~ nclust, scales="free_y") + 
  scale_colour_viridis_d()

#11. Save----
write.csv(dat.final, "Data/LBCUKDEClusters.csv", row.names=FALSE)
