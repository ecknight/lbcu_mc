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
covs <- read.csv("Data/LBCU_environvars.csv") %>% 
  mutate(region = case_when(nclust=="3" & group==1 ~ "Central",
                            nclust=="3" & group==2 ~ "West",
                            nclust=="3" & group==3 ~ "East",
                            nclust=="manual" & group==1 ~ "West",
                            nclust=="manual" & group==2 ~ "Central",
                            nclust=="manual" & group==3 ~ "Inland east",
                            nclust=="manual" & group==4 ~ "Coastal east"))

#2. Make long-----
covs.long <- covs %>% 
  pivot_longer(built:wetland, names_to="var", values_to="val") %>% 
  mutate(seasonvar = paste0(var, "_", season))

#3. Visualize----
ggplot(covs.long) +
  geom_density_ridges(aes(x=val, y = nclust, fill=region), alpha = 0.5) +
  facet_grid(season ~ var, scales="free")

ggplot(covs.long) +
  geom_boxplot(aes(y=val, x = nclust, fill=region), alpha = 0.5) +
  facet_grid(var ~ season, scales="free")

#4. Set up loop----
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

#5. Look at results----
#use p < 0.01 threshold based on visual inspection of plots
out.use <- out %>% 
  separate(seasonvar, into=c("var", "season")) %>% 
  dplyr::filter(p < 0.01)

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
raw <- read.csv("Data/LBCU_environvars_region.csv")  %>% 
  mutate(region = case_when(nclust=="3" & group==1 ~ "central",
                            nclust=="3" & group==2 ~ "west",
                            nclust=="3" & group==3 ~ "east",
                            nclust=="manual" & group==1 ~ "central",
                            nclust=="manual" & group==2 ~ "west",
                            nclust=="manual" & group==3 ~ "east - inland",
                            nclust=="manual" & group==4 ~ "east - coastal")) %>% 
  dplyr::select(nclust, region, season, year, built, change, crop, drought, grass, wetland, seasonality)
  
#4. Link attribute data to bbs data----
dat <- inner_join(bbs, raw)

#5. Detrend----
#https://github.com/crushing05/WOTH_retro/blob/master/WOTH_data_prep.R
dat <- plyr::ddply(dat, c("nclust", "region", "season"), mutate, built = resid(lm(built~year)))
dat <- plyr::ddply(dat, c("nclust", "region", "season"), mutate, change = resid(lm(change~year)))
dat <- plyr::ddply(dat, c("nclust", "region", "season"), mutate, crop = resid(lm(crop~year)))
dat <- plyr::ddply(dat, c("nclust", "region", "season"), mutate, drought = resid(lm(drought~year)))
dat <- plyr::ddply(dat, c("nclust", "region", "season"), mutate, grass = resid(lm(grass~year)))
dat <- plyr::ddply(dat, c("nclust", "region", "season"), mutate, wetland = resid(lm(wetland~year)))
dat <- plyr::ddply(dat, c("nclust", "region", "season"), mutate, seasonality = resid(lm(seasonality~year)))

#6. Pivot wider by season----
dat.w <- dat %>% 
  pivot_wider(values_from=built:seasonality, names_from="season")

#7. Set up loops----
loop <- dat.w %>%
  dplyr::select(nclust, region) %>% 
  unique()

mods <- list()
for(i in 1:nrow(loop)){
  
  #8. Filter to cluster & region of interest----
  dat.i <- dat.w %>% 
    dplyr::filter(nclust==loop$nclust[i],
                  region==loop$region[i])
  
  #9. Scale covariates----
  dat.s <- dat.i %>% 
    mutate(year = year - min(year) + 1) %>% 
    mutate(across(built_breed:seasonality_winter, scale)) %>% 
    mutate(across(built_breed:seasonality_winter, as.numeric))
  
  #10. Get list of covs that are different----
  out.i <- dplyr::filter(out.use, nclust==loop$nclust[i]) %>% 
    mutate(varseason = paste0(var, "_", season))
  
  #9. Check for covariance----
  covs.use <- dat.s %>% 
    dplyr::select(out.i$varseason) %>% 
    data.frame()
  
  covs.vif <- vifstep(covs.use, th=10)
  
  #10. Make formula----
  f <- as.formula(paste("count ~ year + ", paste(covs.vif@results$Variables, collapse= "+"), " + (1|id.route)"))
  
  mods[[i]] <- brm(f,
                   family=negbinomial(),
                   data = dat.s, 
                   warmup = 1000, 
                   iter   = 5000, 
                   chains = 3, 
                   inits  = "random",
                   cores  = 6)
  
}
summary(mods[[1]])
summary(mods[[2]])
summary(mods[[3]])
summary(mods[[4]])
summary(mods[[5]])
summary(mods[[6]])
summary(mods[[7]])

betas <- data.frame()
for(i in 1:length(mods)){
  
  betas <- summary(mods[[i]])$fixed %>% 
    rename(l = 'l-95% CI', u = 'u-95% CI') %>% 
    mutate(effect = case_when(l < 0 & u < 0 ~ "negative",
                              l > 0 & u > 0 ~ "positive",
                              l < 0 & u > 0 ~ "none"),
           nclust=loop$nclust[i],
           region=loop$region[i],
           seasonvar=row.names(summary(mods[[i]])$fixed)) %>% 
    separate(seasonvar, into=c("season", "var")) %>% 
    rbind(betas)
  
}

betas.cov <- dplyr::filter(betas, !is.na(var),
                           abs(Estimate) < 10)
betas.year <- dplyr::filter(betas, season=="year")

ggplot(betas.year) +
  geom_errorbar(aes(x=region, ymin=l, ymax=u)) +
  facet_wrap(~nclust, scales="free")

ggplot(betas.cov %>% dplyr::filter(effect!="none")) +
  geom_hline(aes(yintercept=0)) +
  geom_errorbar(aes(x=var, ymin=l, ymax=u, colour=season)) +
  facet_wrap(region~nclust, scales="free")
