library(tidyverse)
library(lubridate)
library(brms)
library(performance)
library(usdm)

options(scipen=99999)

#Variables
#Spring vars - t, t+1
#Breed vars - t+1
#Fall vars - t+1
#Winter vars - t+1

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
  mutate(years = (year-2007)/13)

#4. Check for covariance----
covs <- dat %>% 
  pivot_wider(names_from=c("season", "var"), values_from="val", values_fill=NA, names_sort=TRUE) %>% 
  dplyr::select(breed_HRarea:winter_WinterHRs) %>% 
  data.frame()%>% 
  dplyr::select(-fallmig_dist,
                -springmig_dist,
                -springmigt_dist)

covs.cor <- data.frame(cor(covs, use="na.or.complete"))
covs.cor

vif(covs)

#5. MLE tests for univariate significance----
loops <- dat %>% 
  dplyr::select(nclust, region, season, var) %>% 
  unique()

priors <- c(prior(normal(0,10), class = "Intercept"),
            prior(normal(0,10), class = "b", coef ="years"),
            prior(normal(0,10), class = "b", coef ="val"))

mod.df <- data.frame()
for(i in 2:nrow(loops)){
  
  dat.i <- dat %>% 
    dplyr::filter(nclust==loops$nclust[i],
                  region==loops$region[i],
                  season==loops$season[i],
                  var==loops$var[i]) %>% 
    dplyr::filter(!is.na(val))
  
  mod.i <- brm(count ~ years + val + (1|id.route) + (1|id.bird),
                  family=negbinomial(),
                  data = dat.i, 
                  warmup = 1000, 
                  iter   = 5000, 
                  chains = 3, 
                  inits  = "random",
                  cores  = 6,
                  prior=priors)
  
  mod.df <- summary(mod.i, prob = 0.80)$fixed %>% 
    mutate(nclust=loops$nclust[i],
           region=loops$region[i],
           season=loops$season[i],
           var=loops$var[i]) %>% 
    rbind(mod.df)

  print(paste0("Finished loop ", i, " of ", nrow(loops)))
}

#6. Select covariates----

#7. Full models----
