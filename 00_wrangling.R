library(tidyverse)
library(sf)

#1. Read in data----
#remove individuals on the Atlantic coast or that didn't complete a full season
dat <- read.csv("Data/LBCU_FilteredData_Segmented_Manuscript.csv")

#2. Breeding ground means----
breed.mn <- dat |> 
  dplyr::filter(season=="breed") |> 
  group_by(study, id, year) |> 
  summarize(lat = mean(lat),
            lon = mean(lon),
            days=n()) |> 
  ungroup() |> 
  mutate(season="breed",
         cluster=1)

#3. Wintering ground means----
winter.mn <- dat |> 
  dplyr::filter(season=="winter", winter!="wintermig") |> 
  mutate(year = ifelse(doy < 130, year-1, year),
         cluster = as.numeric(str_sub(winter, -1, -1))) |> 
  group_by(study, id, year, cluster) |> 
  summarize(lat = mean(lat),
            lon = mean(lon),
            days=n()) |> 
  ungroup() |> 
  mutate(season="winter")

#4. Stopover means----
stop.mn <- dat |> 
  dplyr::filter(stopover==1) |> 
  group_by(study, id, season, year, stopovercluster) |> 
  summarize(lat = mean(lat),
            lon = mean(lon),
            days=n()) |> 
  rename(cluster = stopovercluster) |> 
  group_by(id, season, year) |> 
  mutate(cluster = row_number()) |> 
  ungroup()

#5. Put together----
mn <- rbind(breed.mn, winter.mn, stop.mn)
table(mn$season, mn$cluster)

#6. Transform to utm----
mn.utm <- st_as_sf(mn, coords=c("lon", "lat"), crs=4326) |> 
  st_transform(crs=3857) |> 
  st_coordinates() |> 
  cbind(mn) |> 
  rename(seasoncluster=cluster)

write.csv(mn.utm, "Data/LBCUMCLocations.csv", row.names = FALSE)
