library(tidyverse)
library(vegan)
library(lme4)
library(MuMIn)
library(sf)
library(adehabitatLT)
library(data.table)

options(scipen=99999)

n <- 3

#1. Get migration characteristic data----
dat <- read.csv("Data/MovementBehaviours.csv")

#2. Visualize----
#Departure
ggplot(dat %>% dplyr::filter(var=="depart")) +
  geom_boxplot(aes(x=region, y=val)) +
  facet_wrap(~season, scales="free")

#Arrival
ggplot(dat %>% dplyr::filter(var=="arrive")) +
  geom_boxplot(aes(x=region, y=val)) +
  facet_wrap(~season, scales="free")

#Duration
ggplot(dat %>% dplyr::filter(var=="duration")) +
  geom_boxplot(aes(x=region, y=val)) +
  facet_wrap(~season, scales="free")

#Distance
ggplot(dat %>% dplyr::filter(var=="dist")) +
  geom_boxplot(aes(x=region, y=log(val))) +
  facet_wrap(~season, scales="free")

#Rate
ggplot(dat %>% dplyr::filter(var=="rate")) +
  geom_boxplot(aes(x=region, y=log(val))) +
  facet_wrap(~season, scales="free")

#Number of stopovers
ggplot(dat %>% dplyr::filter(var=="Stopovers")) +
  geom_histogram(aes(x=val)) +
  facet_grid(region~season)

#Area
ggplot(dat %>% dplyr::filter(var=="HRarea")) +
  geom_boxplot(aes(x=region, y=log(val))) +
  facet_wrap(~season,scales="free")

#Wintering home ranges
ggplot(dat %>% dplyr::filter(var=="WinterHRs")) +
  geom_histogram(aes(x=val)) +
  facet_wrap(region~season)

#3. Test----
aic.list <- list()

#Departure
dat.dep <- dat %>% dplyr::filter(var=="depart")
m1 <- lmer(val ~ region*season + (1|id), data=dat.dep, na.action="na.fail", REML=FALSE)
aic.list[[1]] <- dredge(m1)
aic.list[[1]]$mod = "departure"
aic.list[[1]]
summary(m1)
plot(m1)
pdep <- data.frame(pred = predict(m1)) %>% 
  cbind(dat.dep)

ggplot(dat.dep) +
  geom_boxplot(aes(x=region, y=val)) +
  geom_boxplot(aes(x=region, y=pred), data=pdep, fill=NA, colour="red") +
  facet_wrap(~season, scales="free")

#Arrival
m2 <- lmer(arrive ~ region*season + (1|id), data=dat.dur, na.action="na.fail", REML=FALSE)
aic.list[[2]] <- dredge(m2)
aic.list[[2]]$mod = "arrival"
aic.list[[2]]
summary(m2)
plot(m2)

parr <- data.frame(pred = predict(m2)) %>% 
  cbind(dat.dur)

ggplot(dat.dur) +
  geom_boxplot(aes(x=region, y=arrive)) +
  geom_boxplot(aes(x=region, y=pred), data=parr, fill=NA, colour="red") +
  facet_wrap(~season, scales="free")

#Duration
m3 <- lmer(duration ~ region*season + (1|id), data=dat.dur, na.action="na.fail", REML=FALSE)
aic.list[[3]] <- dredge(m3)
aic.list[[3]]$mod = "duration"
aic.list[[3]]
summary(m3)
plot(m3)

pdur <- data.frame(pred = predict(m3)) %>% 
  cbind(dat.dur)

ggplot(dat.dur) +
  geom_boxplot(aes(x=region, y=duration)) +
  geom_boxplot(aes(x=region, y=pred), data=pdur, fill=NA, colour="red") +
  facet_wrap(~season, scales="free")

#Distance
m4 <- lmer(dist ~ region*season + (1|id), data=dat.dist, na.action="na.fail", REML=FALSE)
aic.list[[4]] <- dredge(m4)
aic.list[[4]]$mod = "distance"
aic.list[[4]]
summary(m4)
plot(m4)

pdist <- data.frame(pred = predict(m4)) %>% 
  cbind(dat.dist)

ggplot(dat.dist) +
  geom_boxplot(aes(x=region, y=dist)) +
  geom_boxplot(aes(x=region, y=pred), data=pdist, fill=NA, colour="red") +
  facet_wrap(~season, scales="free")

ggplot(dat.dist) +
  geom_boxplot(aes(x=season, y=dist)) +
  geom_boxplot(aes(x=season, y=pred), data=pdist, fill=NA, colour="red") +
  facet_wrap(~region, scales="free")

#Rate
m5 <- lmer(rate ~ region*season + (1|id), data=dat.rate, na.action="na.fail", REML=FALSE)
aic.list[[5]] <- dredge(m5)
aic.list[[5]]$mod = "rate"
aic.list[[5]]
summary(m5)
plot(m5)

prate <- data.frame(pred = predict(m5)) %>% 
  cbind(dat.rate)

ggplot(dat.rate) +
  geom_boxplot(aes(x=region, y=rate)) +
  geom_boxplot(aes(x=region, y=pred), data=prate, fill=NA, colour="red") +
  facet_wrap(~season, scales="free")

#Stopovers
m6 <- glmer(n ~ region*season + (1|id), data=dat.mig, na.action="na.fail", family="poisson")
aic.list[[6]] <- dredge(m6)
aic.list[[6]]$mod = "stopover"
aic.list[[6]]
m6 <- glmer(n ~ region + season + (1|id), data=dat.mig, na.action="na.fail", family="poisson")
summary(m6)
plot(m6)

pstop <- data.frame(pred = predict(m6, type="response")) %>% 
  cbind(dat.mig)

ggplot(dat.mig) +
  geom_boxplot(aes(x=region, y=n)) +
  geom_boxplot(aes(x=region, y=pred), data=pstop, fill=NA, colour="red") +
  facet_wrap(~season, scales="free")

#Home range area
m7 <- lmer(log(area) ~ region*season + (1|bird), data=dat.hr, na.action="na.fail", REML = FALSE)
aic.list[[7]] <- dredge(m7)
aic.list[[7]]$mod = "area"
aic.list[[7]]
plot(m7)
summary(m7)

phr <- data.frame(pred = predict(m7)) %>% 
  cbind(dat.hr)

ggplot(dat.hr) +
  geom_boxplot(aes(x=region, y=log(area))) +
  geom_boxplot(aes(x=region, y=pred), data=phr, fill=NA, colour="red") +
  facet_wrap(~season, scales="free")

#Number of home ranges
m8 <- glmer(n ~ region + (1|bird), data=dat.wint, na.action="na.fail", family="poisson")
aic.list[[8]] <- dredge(m8)
aic.list[[8]]$mod = "winterhrs"
aic.list[[8]]
plot(m8)
summary(m8)

#Stopover duration
m9 <- lmer(n ~ region*season + (1|id), data=dat.stop, na.action="na.fail", REML=FALSE)
aic.list[[9]] <- dredge(m9)
aic.list[[9]]
m9 <- lmer(n ~ season + (1|id), data=dat.stop, na.action="na.fail", REML=FALSE)
plot(m9)
summary(m9)

pstop <- data.frame(pred = predict(m9, type="response")) %>% 
  cbind(dat.stop)

ggplot(dat.stop) +
  geom_boxplot(aes(x=region, y=n)) +
  geom_boxplot(aes(x=region, y=pred), data=pstop, fill=NA, colour="red") +
  facet_wrap(~season, scales="free")

#11. Collapse AIC results----
aic.out <- rbindlist(aic.list, fill=TRUE) %>% 
  arrange(mod, -weight) %>% 
  mutate(weight = round(weight, 3))
