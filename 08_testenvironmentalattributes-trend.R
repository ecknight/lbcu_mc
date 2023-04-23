library(tidyverse)
library(brms)

options(scipen=99999)

#1. Get bbs cluster assignment data----
clust <- read.csv("Data/LBCUBBSClusters.csv") %>% 
  dplyr::filter(nclust %in% c("3", "manual")) %>% 
  rename(group=knncluster, id.route = id) %>% 
  mutate(region = case_when(nclust=="3" & group==1 ~ "central",
                            nclust=="3" & group==2 ~ "west",
                            nclust=="3" & group==3 ~ "east",
                            nclust=="manual" & group==1 ~ "central",
                            nclust=="manual" & group==2 ~ "west",
                            nclust=="manual" & group==3 ~ "east - inland",
                            nclust=="manual" & group==4 ~ "east - coastal")) %>%
  dplyr::select(id.route, nclust, region)

#2. Link cluster to bbs monitoring data----
load("~/Library/Application Support/bbsBayes/bbs_raw_data.RData")

bbs.bird <- bbs_data[["bird"]] %>% 
  dplyr::filter(AOU==2640,
                Year >= 2007)

bbs <- bbs_data$route %>% 
  left_join(bbs.bird) %>% 
  mutate(count = ifelse(is.na(SpeciesTotal), 0, SpeciesTotal),
         id.route=paste(countrynum, statenum, Route, sep="-")) %>% 
  inner_join(clust) %>% 
  rename(route = Route, lat = Latitude, lon = Longitude, year = Year) %>% 
  dplyr::select(id.route, year, count, nclust, region)

#3. Link attribute data to bbs data----

raw <- read.csv("Data/LBCU_environvars.csv") %>% 
  mutate(conv = covcrop + covbuilt,
         region = case_when(kdecluster==1 ~ "central", 
                            kdecluster==2 ~ "west",
                            kdecluster==3 ~ "east")) %>% 
  dplyr::select(-kdecluster, -season, -year, -cluster) %>% 
  separate(id, into=c("kdecluster", "season", "birdid", "year", "cluster"), remove=FALSE)

raw <- read.csv("Data/MovementBehaviours.csv") %>% 
  dplyr::filter(!var %in% c("arrive2", "depart2")) %>% 
  group_by(nclust, region, id, year, season, var) %>% 
  mutate(row = row_number()) %>% 
  ungroup() %>% 
  mutate(year = year+1)

springt <- raw %>% 
  dplyr::filter(season=="springmig") %>% 
  mutate(year = year-1,
         season = paste0(season, "t"))

dat <- rbind(raw, springt) %>%  
  rename(id.bird = id) %>% 
  inner_join(bbs) %>% 
  group_by(season, var) %>% 
  mutate(val = as.numeric(scale(val))) %>% 
  ungroup() %>% 
  pivot_wider(names_from=c("season", "var"), values_from="val", values_fill=NA) %>% 
  mutate(years = (year-2007)/13)
