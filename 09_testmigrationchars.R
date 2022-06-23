library(tidyverse)
library(vegan)
library(lme4)
library(MuMIn)
library(sf)

options(scipen=99999)

n <- 2

#1. Get dominant cluster id for each bird----
clust <- read.csv("Data/LBCUKDEClusters.csv") %>% 
  dplyr::filter(nclust==n,
                !is.na(X)) %>% 
  group_by(id) %>% 
  summarize(kdecluster = round(mean(kdecluster))) %>% 
  ungroup()

#2. Read in daily estimate data---
raw <- read.csv("Data/LBCU_FilteredData_Segmented.csv") %>% 
  dplyr::filter(!id %in% c(46768277, 33088, 129945787)) %>% 
  left_join(clust) %>% 
  mutate(region = case_when(kdecluster==1 ~ "east", 
                            kdecluster==2 ~ "west")) %>% 
  dplyr::filter(!is.na(region)) 

#3. Wrangle migration duration----
dat.dur <- raw %>% 
  dplyr::filter(segment %in% c("depart", "arrive")) %>% 
  mutate(season = case_when(season=="breed" ~ "springmig",
                            season=="winter" ~ "fallmig",
                            !is.na(season) ~ season)) %>% 
  pivot_wider(id_cols=c(region, id, season, year), names_from=segment, values_from=doy) %>% 
  group_by(season) %>% 
  mutate(duration = arrive - depart, 
         arrive = arrive - min(arrive, na.rm=TRUE),
         depart = depart - min(depart, na.rm=TRUE)) %>% 
  dplyr::filter(!is.na(duration)) %>% 
  ungroup() %>% 
  dplyr::filter(duration < 70)

write.csv(dat.dur, "Data/MigrationTiming.csv", row.names = FALSE)
#Remove two outliers that are due to gaps in transmission

#4. Wrangle number of stopovers/home ranges----
dat.mig <- raw %>% 
  dplyr::filter(season %in% c("springmig", "fallmig")) %>% 
  group_by(region, id, year, season) %>% 
  summarize(n = max(stopovercluster, na.rm=TRUE)) %>% 
  ungroup() %>% 
  mutate(n = ifelse(is.infinite(n), 0, n))

dat.wint <- raw %>% 
  dplyr::filter(season=="winter", winter!="wintermig") %>% 
  mutate(year = as.numeric(ifelse(doy < 130, year-1, year)),
         cluster = as.numeric(str_sub(winter, -1, -1))) %>% 
  group_by(season, region, id, year) %>% 
  summarize(n = max(cluster)) %>% 
  ungroup()

dat.stop <- dat.mig %>% 
  rbind(dat.wint)

#5. Wrangle home range size----
dat.hr <- read_sf("gis/shp/kde_individual.shp") %>% 
  mutate(kdecluster = as.numeric(kdecluster),
         id = as.numeric(id),
         year = as.numeric(year)) %>% 
  left_join(clust) %>% 
  data.frame() %>% 
  dplyr::select(-geometry) %>% 
  mutate(region = case_when(kdecluster==1 ~ "east", 
                            kdecluster==2 ~ "west")) %>% 
  dplyr::filter(area < 10000) %>% 
  mutate(area.s = scale(area))

#6. Visualize
#Departure
ggplot(dat.dur) +
  geom_boxplot(aes(x=region, y=depart)) +
  facet_wrap(~season, scales="free")

#Arrival
ggplot(dat.dur) +
  geom_boxplot(aes(x=region, y=arrive)) +
  facet_wrap(~season, scales="free")

#Duration
ggplot(dat.dur) +
  geom_boxplot(aes(x=region, y=duration)) +
  facet_wrap(~season)

#Number of stopovers
ggplot(dat.mig) +
  geom_histogram(aes(x=n)) +
  facet_grid(region~season)

#Area
ggplot(dat.hr) +
  geom_boxplot(aes(x=region, y=area.s)) +
  facet_wrap(~season)

ggplot(dat.hr) +
  geom_histogram(aes(x=area)) +
  facet_grid(region~season)

#7. Test----

#Departure
m1 <- lmer(depart ~ region*season + (1|year) + (1|id), data=dat.dur, na.action="na.fail", REML=FALSE)
dredge(m1)
summary(m1)
plot(m1)
p <- data.frame(pred = predict(m1)) %>% 
  cbind(dat.dur)

ggplot(dat.dur) +
#  geom_boxplot(aes(x=region, y=depart)) +
  geom_boxplot(aes(x=region, y=pred), data=p, fill=NA, colour="red") +
  facet_wrap(~season, scales="free")

#Arrival
m2 <- lmer(arrive ~ region*season + (1|year) + (1|id), data=dat.dur, na.action="na.fail", REML=FALSE)
dredge(m2)
m2 <- lmer(arrive ~ region + season + (1|year) + (1|id), data=dat.dur, na.action="na.fail", REML=FALSE)
summary(m2)
plot(m2)

#Duration
m3 <- lmer(duration ~ region*season + (1|year) + (1|id), data=dat.dur, na.action="na.fail", REML=FALSE)
dredge(m3)
summary(m3)
plot(m3)

#Stopovers
m4 <- glmer(n ~ region*season + (1|year) + (1|id), data=dat.mig, na.action="na.fail", family="poisson")
dredge(m4)
m5 <- glmer(n ~ region + season + (1|year) + (1|id), data=dat.mig, na.action="na.fail", family=poisson)
dredge(m5)
summary(m5)
plot(m5)

#Home range
dat.hr.w <- dat.hr %>% 
  dplyr::filter(season=="winter")

m6 <- glmer(area ~ region + (1|id) + (1|year), data=dat.hr.w, na.action="na.fail", family=Gamma)
dredge(m6)
plot(m6)
summary(m6)

#Need to think about this....




#8. Playing with other things----

#Carry over effects
ggplot(dat.dur) +
  geom_point(aes(x=depart, y=duration, colour=region)) +
  facet_wrap(region~season, scales="free")

ggplot(dat.dur) +
  geom_point(aes(x=arrive, y=duration, colour=region)) +
  facet_wrap(region~season, scales="free")

ggplot(dat.dur) +
  geom_point(aes(x=depart, y=arrive, colour=region)) +
  facet_wrap(region~season, scales="free")

m3 <- lmer(duration ~ depart*region + depart*season + (1|year/id), na.action="na.fail", data=dat.dur)
dredge(m3)

m4 <- lmer(duration ~ arrive*region + arrive*season + (1|year/id), na.action="na.fail", data=dat.dur)
dredge(m4)

#Link to environmental factors
raw.cov <- read.csv("Data/LBCU_environvars.csv") %>% 
  mutate(region = case_when(kdecluster==1 ~ "east", 
                            kdecluster==2 ~ "west"),
         conv = covcrop + covbuilt) %>% 
  group_by(region, id, year, season) %>% 
  summarize(crops = mean(crops),
            flooded_vegetation = mean(flooded_vegetation),
            grass = mean(grass),
            occurrence = mean(occurrence),
            pdsi = mean(pdsi),
            recurrence = mean(recurrence),
            seasonality = mean(seasonality),
            conv = mean(conv)) %>% 
  ungroup()

dat.cov <- dat.mig %>% 
  left_join(raw.cov) %>% 
  left_join(raw.cov %>% 
              dplyr::filter(season=="winter") %>% 
              rename(crops.w = crops, flooded_vegetation.w = flooded_vegetation, grass.w = grass, occurrence.w = occurrence, pdsi.w = pdsi, recurrence.w = recurrence, seasonality.w = seasonality, conv.w = conv) %>% 
              dplyr::select(-season))
  
ggplot(dat.cov) +
  geom_smooth(aes(x=recurrence.w, y=duration, colour=season)) +
  geom_point(aes(x=recurrence.w, y=duration, colour=season)) +
  facet_grid(region ~ season, scales="free")