library(tidyverse)
library(lubridate)
library(lme4)
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
            prior(normal(0,10), class = "b", coef ="var"))

mod.df <- data.frame()
for(i in 1:nrow(loops)){
  
  dat.i <- dat %>% 
    dplyr::filter(nclust==loops$nclust[i],
                  region==loops$region[i],
                  season==loops$season[i],
                  var==loops$var[i]) %>% 
    dplyr::filter(!is.na(val))
  
  #Get an estimate of theta
  fit <- glmer.nb(count ~ years  + (1|id.route) + (1|id.bird), data=dat.i)
  theta <- lme4:::getNBdisp(fit)
  
  fit1 <- try(glmer.nb(count ~ years + fallmig_arrive + (1|id.route) + (1|id.bird), data=dat.i))
  fit2 <- try(glmer.nb(count ~ years + fallmig_depart + (1|id.route) + (1|id.bird), data=dat.i))
  fit4 <- try(glmer.nb(count ~ years + fallmig_duration + (1|id.route) + (1|id.bird), data=dat.i))
  fit5 <- try(glmer.nb(count ~ years + fallmig_rate + (1|id.route) + (1|id.bird), data=dat.i))
  fit6 <- try(glmer.nb(count ~ years + fallmig_HRarea + (1|id.route) + (1|id.bird), data=dat.i))
  fit7 <- try(glmer.nb(count ~ years + fallmig_Stopovers + (1|id.route) + (1|id.bird), data=dat.i))
  fit8 <- try(glmer.nb(count ~ years + springmig_arrive + (1|id.route) + (1|id.bird), data=dat.i))
  fit9 <- try(glmer.nb(count ~ years + springmig_depart + (1|id.route) + (1|id.bird), data=dat.i))
  fit11 <- try(glmer.nb(count ~ years + springmig_duration + (1|id.route) + (1|id.bird), data=dat.i))
  fit12 <- try(glmer.nb(count ~ years + springmig_rate + (1|id.route) + (1|id.bird), data=dat.i))
  fit13 <- try(glmer.nb(count ~ years + springmig_HRarea + (1|id.route) + (1|id.bird), data=dat.i))
  fit14 <- try(glmer.nb(count ~ years + springmig_Stopovers + (1|id.route) + (1|id.bird), data=dat.i))
  fit15 <- try(glmer.nb(count ~ years + springmigt_arrive + (1|id.route) + (1|id.bird), data=dat.i))
  fit16 <- try(glmer.nb(count ~ years + springmigt_depart + (1|id.route) + (1|id.bird), data=dat.i))
  fit18 <- try(glmer.nb(count ~ years + springmigt_duration + (1|id.route) + (1|id.bird), data=dat.i))
  fit19 <- try(glmer.nb(count ~ years + springmigt_rate + (1|id.route) + (1|id.bird), data=dat.i))
  fit20 <- try(glmer.nb(count ~ years + springmigt_HRarea + (1|id.route) + (1|id.bird), data=dat.i))
  fit21 <- try(glmer.nb(count ~ years + springmigt_Stopovers + (1|id.route) + (1|id.bird), data=dat.i))
  fit22 <- try(glmer.nb(count ~ years + breed_HRarea + (1|id.route) + (1|id.bird), data=dat.i))
  fit23 <- try(glmer.nb(count ~ years + winter_HRarea + (1|id.route) + (1|id.bird), data=dat.i))
  fit24 <- try(glmer.nb(count ~ years + winter_WinterHRs + (1|id.route) + (1|id.bird), data=dat.i))
  
  mods <- list(fit1, fit2, fit4, fit5, fit7, fit8, fit9, fit11, fit12, fit13, fit14, fit15, fit16, fit18, fit19, fit21, fit22, fit24)
  
  for(j in 1:length(mods)){
    
    if(class(mods[[j]])=="glmerMod"){
      mod.df <- data.frame(summary(mods[[j]])$coefficients) %>% 
        mutate(nclust=loops$nclust[i],
               region=loops$region[i],
               cov = row.names(summary(mods[[j]])$coefficients)) %>% 
        dplyr::filter(!cov %in% c("(Intercept)", "years")) %>% 
        rbind(mod.df)
    }
  }
  print(paste0("Finished loop ", i, " of ", nrow(loops)))
}

#6. Select covariates----

#7. Full models----
