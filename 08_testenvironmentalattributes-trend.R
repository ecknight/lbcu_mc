library(tidyverse)
library(brms)
library(ggridges)

options(scipen=99999)

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
  mutate(seasonvar = paste0(season, "_", var))

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
  separate(seasonvar, into=c("season", "var")) %>% 
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

#3. Link attribute data to bbs data----
#Spring vars - t, t+1
#Breed vars - t+1
#Fall vars - t+1
#Winter vars - t+1

raw <- read.csv("Data/LBCU_environvars.csv") %>% 
  mutate(conv = covcrop + covbuilt,
         region = case_when(kdecluster==1 ~ "central", 
                            kdecluster==2 ~ "west",
                            kdecluster==3 ~ "east")) %>% 
  dplyr::select(-kdecluster, -season, -year, -cluster) %>% 
  separate(id, into=c("kdecluster", "season", "birdid", "year", "cluster"), remove=FALSE)

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

#5. Test for differences----