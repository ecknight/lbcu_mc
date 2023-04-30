library(tidyverse)
library(brms)
library(ggridges)
library(usdm)
library(lme4)

options(scipen=99999)

#Spring vars - t, t+1
#Breed vars - t+1
#Fall vars - t+1
#Winter vars - t+1

#A. TEST FOR DIFFERENCES####

#1. Get environmental attribute data----
raw <- read.csv("Data/MovementBehaviours.csv") %>% 
  dplyr::filter(!var %in% c("arrive2", "depart2")) %>% 
  group_by(nclust, region, id, year, season, var) %>% 
  mutate(row = row_number()) %>% 
  ungroup() %>% 
  mutate(year = year+1)

raw <- read.csv("Data/LBCU_environvars.csv") %>% 
  mutate(region = case_when(nclust=="3" & group==1 ~ "central",
                            nclust=="3" & group==2 ~ "west",
                            nclust=="3" & group==3 ~ "east",
                            nclust=="manual" & group==1 ~ "central",
                            nclust=="manual" & group==2 ~ "west",
                            nclust=="manual" & group==3 ~ "east - inland",
                            nclust=="manual" & group==4 ~ "east - coastal"))

#2. Add t-1 for spring migration chars----
springt <- raw %>% 
  dplyr::filter(season=="springmig") %>% 
  mutate(year = year-1,
         season = paste0(season, "t"))

covs <- rbind(raw, springt)

#3. Make long-----
covs.long <- covs %>% 
  pivot_longer(built:wetland, names_to="var", values_to="val") %>% 
  mutate(seasonvar = paste0(var, "_", season))

#4. Visualize----
ggplot(covs.long) +
  geom_density_ridges(aes(x=val, y = nclust, fill=region), alpha = 0.5) +
  facet_grid(season ~ var, scales="free")

ggplot(covs.long) +
  geom_boxplot(aes(y=val, x = nclust, fill=region), alpha = 0.5) +
  facet_grid(var ~ season, scales="free")


#5. Set up loop----
loop <- expand.grid(seasonvar = unique(covs.long$seasonvar),
                    nclust = unique(covs.long$nclust))

out <- data.frame()
for(i in 1:nrow(loop)){
  
  covs.i <- covs.long %>% 
    dplyr::filter(seasonvar==loop$seasonvar[i],
                  nclust==loop$nclust[i])
  
  test.i <- kruskal.test(covs.i$val, covs.i$group)
  
  out <- rbind(out,
               cbind(loop[i,],
                     data.frame(statistic = test.i$statistic,
                                p = test.i$p.value)))
  
}

#6. Look at results----
#use p < 0.01 threshold based on visual inspection of plots
out.use <- out %>% 
  separate(seasonvar, into=c("var", "season")) %>% 
  dplyr::filter(p < 0.001)

table(out.use$var, out.use$season, out.use$nclust)

#B. TEST FOR RELATIONSHIP WITH BBS DATA####

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

#3. Get environmental attributes for each region----
raw.reg <- read.csv("Data/LBCU_environvars_region.csv")  %>% 
  mutate(region = case_when(nclust=="3" & group==1 ~ "central",
                            nclust=="3" & group==2 ~ "west",
                            nclust=="3" & group==3 ~ "east",
                            nclust=="manual" & group==1 ~ "central",
                            nclust=="manual" & group==2 ~ "west",
                            nclust=="manual" & group==3 ~ "east - inland",
                            nclust=="manual" & group==4 ~ "east - coastal"))

#4. Add t-1 for spring migration attrs----
springt.reg <- raw.reg %>% 
  dplyr::filter(season=="springmig") %>% 
  mutate(year = year-1,
         season = paste0(season, "t"))
  
#5. Link attribute data to bbs data----
dat <- rbind(raw.reg, springt.reg) %>%  
  inner_join(bbs)

#6. Detrend----
#https://github.com/crushing05/WOTH_retro/blob/master/WOTH_data_prep.R
dat <- plyr::ddply(dat, c("nclust", "region", "season"), mutate, built = resid(lm(built~year)))
dat <- plyr::ddply(dat, c("nclust", "region", "season"), mutate, change = resid(lm(change~year)))
dat <- plyr::ddply(dat, c("nclust", "region", "season"), mutate, crop = resid(lm(crop~year)))
dat <- plyr::ddply(dat, c("nclust", "region", "season"), mutate, drought = resid(lm(drought~year)))
dat <- plyr::ddply(dat, c("nclust", "region", "season"), mutate, grass = resid(lm(grass~year)))
dat <- plyr::ddply(dat, c("nclust", "region", "season"), mutate, wetland = resid(lm(wetland~year)))

#7. Scale covariates---
dat.s <- dat %>% 
  mutate(built = as.numeric(scale(built)),
         change = as.numeric(scale(change)),
         crop = as.numeric(scale(crop)),
         drought = as.numeric(scale(drought)),
         grass = as.numeric(scale(grass)),
         wetland = as.numeric(scale(wetland)),
         year = as.numeric(scale(year))) %>% 
  dplyr::select(nclust, region, id.route, year, season, count, built, change, crop, drought, grass, wetland)

#8. Pivot wider by season----
dat.w <- dat.s %>% 
  pivot_wider(values_from=built:wetland, names_from="season")

#9. Check for covariance----
covs.use1 <- dat %>% 
  dplyr::filter(nclust=="manual", group==4) %>% 
  dplyr::select(built, change, crop, drought, grass, seasonality, wetland)

covs.cor <- data.frame(cor(covs.use1, use="na.or.complete"))
covs.cor

vif(covs.use1)

#let's take out seasonality
covs.use2 <- covs.use1 %>% 
  dplyr::select(-seasonality)

vif(covs.use2)

#8. Set up priors-----
priors <- c(prior(normal(0,10), class = "Intercept"),
            prior(normal(0,10), class = "b", coef ="years"),
            prior(normal(0,10), class = "b", coef ="val"))

#9. Set up loops----
loop <- dat.w %>%
  dplyr::select(nclust, region) %>% 
  unique()

mods <- list()
for(i in 1:nrow(loop)){
  
  dat.i <- dat.w %>% 
    dplyr::filter(nclust==loop$nclust[i],
                  region==loop$region[i])
  
  if(loop$nclust[i]=="3"){
    mods[[i]] <- brm(count ~ year +
                   # built_breed +
                   # built_fallmig +
                   # built_winter +
                   # change_breed +
                   # change_springmig +
                   # change_springmigt +
                   # change_winter +
                   crop_breed +
                   crop_fallmig +
                   crop_springmig +
                   crop_springmigt +
                   crop_winter +
                   # drought_breed +
                   # drought_winter +
                   # grass_breed +
                   # grass_fallmig +
                   # grass_winter +
                   # wetland_springmig +
                   # wetland_springmigt +
                   (1|id.route),
                 family=negbinomial(),
                 data = dat.i, 
                 warmup = 1000, 
                 iter   = 5000, 
                 chains = 3, 
                 inits  = "random",
                 cores  = 6)
    
    mods[[i]] <- glmer.nb(count ~ year + 
                        built_breed +
                        built_fallmig +
                        built_winter +
                        change_breed +
                        change_springmig +
                        change_springmigt +
                        change_winter +
                        crop_breed +
                        crop_fallmig +
                        crop_springmig +
                        crop_springmigt +
                        crop_winter +
                        drought_breed +
                        drought_winter +
                        grass_breed +
                        grass_fallmig +
                        grass_winter +
                        wetland_springmig +
                        wetland_springmigt +
                        (1|id.route), data=dat.i)
  }
  
  if(loop$nclust[i]=="manual"){
    mods[[i]] <- brm(count ~ year +
                   built_breed +
                   built_fallmig +
                   built_winter +
                   change_springmig +
                   change_springmigt +
                   change_winter +
                   crop_breed +
                   crop_fallmig +
                   crop_springmig +
                   crop_springmigt +
                   crop_winter +
                   drought_breed +
                   drought_winter +
                   grass_breed +
                   grass_springmig +
                   grass_springmigt +
                   grass_fallmig +
                   grass_winter +
                   wetland_fallmig +
                   wetland_springmig +
                   wetland_springmigt +
                   wetland_winter + 
                   (1|id.route),
                 family=negbinomial(),
                 data = dat.i, 
                 warmup = 1000, 
                 iter   = 5000, 
                 chains = 3, 
                 inits  = "random",
                 cores  = 6)
  }

  
  
}