library(tidyverse)
library(vegan)
library(lme4)
library(MuMIn)
library(sf)
library(adehabitatLT)
library(data.table)

options(scipen=99999)

setwd("G:/My Drive/SMBC")

#1. Get id for each bird----
clust <- read.csv("Data/LBCUKDEClusters.csv") %>%
  dplyr::filter(!is.na(X),
                nclust == "3") |> 
  dplyr::select(id, group) |> 
  unique()

#2. Read in daily estimate data----
raw <- read.csv("Data/LBCU_FilteredData_Segmented.csv") %>% 
              dplyr::filter(!id %in% c(46768277, 33088, 129945787, 46770723, 46769927, 86872)) %>% 
  inner_join(clust, multiple="all") %>% 
  arrange(id, date)

ggplot(raw |> sample_n(1000)) + 
  geom_point(aes(x=X, y=Y, colour=factor(group))) +
  facet_grid(~season)

#3. Wrangle migration duration----
dat.dur <- raw %>% 
  dplyr::filter(segment %in% c("depart", "arrive")) %>% 
  mutate(season = case_when(season=="breed" ~ "springmig",
                            season=="winter" ~ "fallmig",
                            !is.na(season) ~ season)) %>% 
  pivot_wider(id_cols=c(id, season, year), names_from=segment, values_from=doy) %>% 
  group_by(season) %>% 
  mutate(duration = arrive - depart, 
         arrive2 = arrive - min(arrive, na.rm=TRUE) + 1,
         depart2 = depart - min(depart, na.rm=TRUE) + 1) %>% 
  dplyr::filter(!is.na(duration)) %>% 
  ungroup() %>% 
  dplyr::filter(duration < 70)
#Remove two outliers that are due to gaps in transmission

ggplot(dat.dur) +
  geom_histogram(aes(x=duration))

ggplot(dat.dur) +
  geom_histogram(aes(x=depart2))

ggplot(dat.dur) +
  geom_histogram(aes(x=arrive2))

#4. Wrangle migration distance----

#Calculate distance
traj <- as.ltraj(xy=raw[,c("X", "Y")],
                 id=raw$id,
                 date=raw$date,
                 typeII=FALSE,
                 proj4string = CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"))
names(traj) <- unique(raw$id)

dat.traj <- rbindlist(traj, idcol=TRUE) %>% 
  rename(id='.id', X=x, Y=y) %>% 
  mutate(id=as.numeric(id)) %>% 
  arrange(id, date) %>% 
  dplyr::select(dist) %>% 
  cbind(raw %>% 
          arrange(id, date)) %>% 
  mutate(dist = dist/1000)

#migrations that weren't completed (from data cleaning script) plus last individual who didn't have data for the first leg of migration
dist.remove <- data.frame(legid=c("1418934379-2021-1spring",
                                  "99900-2021-2fall",
                                  "1146533212-2021-1spring",
                                  "77639184-2016-2fall",
                                  "479364105-2019-1spring",
                                  "46768108-2016-2fall",
                                  "188150741-2021-1spring",
                                  "172070515-2019-2fall",
                                  "1418896449-2021-1spring",
                                  "1418878943-2021-1spring",
                                  "172070319-2016-2fall",
                                  "281981414-2018-1spring",
                                  "290347351-2019-2fall",
                                  "877974163-2020-1spring",
                                  "46768277-2014-2fall",
                                  "145698291-2021-1spring",
                                  "77638376-2015-2fall",
                                  "99902-2021-2fall",
                                  "890834800-2020-1spring",
                                  "1953212690-2011-1spring",
                                  "290352179-2020-1spring")) %>% 
  separate(legid, into=c("id", "year", "season")) %>% 
  mutate(id=as.numeric(id),
         year=as.numeric(year),
         season=ifelse(season=="1spring", "springmig", "fallmig"))

dat.dist <- dat.traj %>% 
  dplyr::filter(season %in% c("fallmig", "springmig") |
                  lead(segment)=="depart" |
                  segment=="arrive") %>% 
  mutate(season = case_when(segment=="stationary" & season=="winter" ~ "springmig",
                            segment=="stationary" & season=="breed" ~ "fallmig",
                            segment=="arrive" & season=="winter" ~ "fallmig",
                            segment=="arrive" & season=="breed" ~ "springmig",
                            !is.na(season) ~ season)) %>% 
  anti_join(dist.remove) %>% 
  group_by(id, year, season) %>% 
  summarize(dist = sum(dist)) %>% 
  ungroup()

hist(dat.dist$dist)
  
#5. Wrangle migration rate----
dat.rate <- dat.dist %>% 
  inner_join(dat.dur) %>% 
  mutate(rate = dist/duration)

ggplot(dat.rate) +
  geom_histogram(aes(x=rate))

#6. Wrangle number of stopovers/home ranges----
dat.mig <- raw %>% 
  dplyr::filter(season %in% c("springmig", "fallmig")) %>% 
  group_by(id, year, season) %>% 
  summarize(n = max(stopovercluster, na.rm=TRUE)) %>% 
  ungroup() %>% 
  mutate(n = ifelse(is.infinite(n), 0, n))

hist(dat.mig$n)

#7. Wrangle stopover duration---
dat.stop <- raw %>% 
  dplyr::filter(season %in% c("springmig", "fallmig"),
                stopover==1) %>% 
  group_by(id, year, season) %>% 
  summarize(n = n()) %>% 
  ungroup()

hist(dat.stop$n)

#8. Wrangle home range size----
dat.hr <- read_sf("gis/shp/kde_individual.shp") |> 
  st_make_valid() %>% 
  st_centroid() %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  dplyr::select(X, Y) %>% 
  cbind(read_sf("gis/shp/kde_individual.shp")) %>% 
  mutate(group = as.numeric(group),
         bird = as.numeric(bird),
         year = as.numeric(year)) %>% 
  left_join(clust %>% 
              rename(bird=id)) %>% 
  data.frame() %>% 
  dplyr::select(-geometry) %>% 
  dplyr::filter(area < 20000) %>% 
  mutate(area.s = scale(area))

hist(dat.hr$area)

#9. Wrangle number of wintering home ranges----
dat.wint <- dat.hr %>% 
  dplyr::filter(season=="winter") %>% 
  group_by(season, bird, year) %>% 
  summarize(n=n(),
            X = mean(X),
            Y = mean(Y)) %>% 
  ungroup()

hist(dat.wint$n)

#9. Put all the data together for plotting----
dat.all <- dat.mig %>% 
  mutate(var = "Stopovers") %>% 
  rename(val=n) %>% 
  rbind(dat.hr %>% 
          dplyr::select(season, bird, year, area) %>% 
          rename(id=bird,
                 val=area) %>% 
          mutate(var = "HRarea")) %>% 
  rbind(dat.wint %>% 
          dplyr::select(bird, year, n) %>% 
          rename(id=bird,
                 val=n) %>% 
          mutate(var = "WinterHRs",
                 season="winter")) |> 
  rbind(dat.rate %>% 
          dplyr::select(-dist, -duration) |> 
          pivot_longer(depart:rate, values_to = "val", names_to="var")) %>% 
  rbind(dat.stop %>% 
          dplyr::select(id, year, season, n) %>% 
          rename(val=n) %>% 
          mutate(var = "stopoverduration"))

#10. Add nclust back in-----
dat.out <- read.csv("Data/LBCUKDEClusters.csv") |> 
  dplyr::filter(nclust %in% c("3", "expert", "flyway")) |> 
  dplyr::select(id, year, season, nclust, group) |> 
  inner_join(dat.all, multiple="all") |> 
  unique()

write.csv(dat.out, "Data/MovementBehaviours.csv", row.names = FALSE)
