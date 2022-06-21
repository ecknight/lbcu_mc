library(tidyverse)
library(vegan)
library(lme4)
library(MuMIn)

options(scipen=99999)

n <- 2

#A. Home range area----

#1. Get dominant cluster id for each bird----
clust <- read.csv("Data/LBCUKDEClusters.csv") %>% 
  dplyr::filter(nclust==n,
                !is.na(X)) %>% 
  group_by(id) %>% 
  summarize(kdecluster = round(mean(kdecluster))) %>% 
  ungroup()

#2. Read in  home range data-----
dat <- read.csv("Data/LBCU_environvars.csv") %>% 
  mutate(region = case_when(kdecluster==1 ~ "east", 
                            kdecluster==2 ~ "west")) %>% 
  dplyr::filter(area< 10000)

#3. Read in daily estimate data---
raw <- read.csv("Data/LBCU_FilteredData_Segmented.csv") %>% 
  dplyr::filter(!id %in% c(46768277, 33088, 129945787)) %>% 
  left_join(clust) %>% 
  mutate(region = case_when(kdecluster==1 ~ "east", 
                            kdecluster==2 ~ "west")) %>% 
  dplyr::filter(!is.na(region))

#4. Wrangle arrival, departure, duration----
arr <- raw %>% 
  dplyr::filter(segment=="arrive")

dep <- raw %>% 
  dplyr::filter(segment=="depart")

duration <- raw %>% 
  dplyr::filter(segment %in% c("depart", "arrive")) %>% 
  mutate(season = case_when(season=="breed" ~ "springmig",
                            season=="winter" ~ "fallmig",
                            !is.na(season) ~ season)) %>% 
  pivot_wider(id_cols=c(region, id, season, year), names_from=segment, values_from=doy) %>% 
  mutate(duration = arrive - depart) %>% 
  dplyr::filter(!is.na(duration))

#5. Wrangle number of stopovers----
dat.stop <- dat %>% 
  dplyr::filter(season %in% c("springmig", "fallmig")) %>% 
  group_by(region, id, year, season) %>% 
  summarize(n = max(cluster)) %>% 
  ungroup() %>% 
  full_join(raw %>% 
              dplyr::filter(season %in% c("springmig", "fallmig")) %>% 
              dplyr::select(region, id, year, season) %>% 
              unique()) %>% 
  mutate(n=ifelse(is.na(n), 0, n))

#6. Put together----
dat.mig <- dat.stop %>% 
  full_join(duration) %>% 
  dplyr::filter(!is.na(duration))

#6. Visualize
#Home range
ggplot(dat) +
  geom_boxplot(aes(x=season, y=log(area), fill=region)) +
  ylab("log of core home range area (km2)")

#Departure
ggplot(dep) +
  geom_boxplot(aes(x=season, y=doy, fill=region))

#Arrival
ggplot(arr) +
  geom_boxplot(aes(x=season, y=doy, fill=region))

#Number of stopovers
ggplot(dat.stop) +
  geom_histogram(aes(x=n)) +
  xlab("Number of stopovers") +
  facet_grid(region~season)

ggplot(dat.stop) +
  geom_boxplot(aes(x=season, y=n, fill=region))

m2 <- glmer(n ~ region*season + (1|id), data=dat.mig, na.action="na.fail", family="poisson")
dredge(m2)

m2 <- glmer(n ~ region + season + (1|id), data=dat.mig, na.action="na.fail", family="poisson")
summary(m2)
plot(m2)

#Duration
ggplot(duration) +
  geom_boxplot(aes(x=season, y=duration, fill=region))

m1 <- lmer(duration ~ region*season + (1|id), data=dat.mig, na.action="na.fail", REML=FALSE)
dredge(m1)
summary(m1)
plot(m1)

#Colinearity
ggplot(dat.mig) +
  geom_smooth(aes(x=duration, y=n, colour=season)) +
  geom_jitter(aes(x=duration, y=n, colour=season))