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
dat.dep <- dat %>% dplyr::filter(var=="depart2")
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
dat.arr <- dat %>% dplyr::filter(var=="arrive2")
m2 <- lmer(val ~ region*season + (1|id), data=dat.arr, na.action="na.fail", REML=FALSE)
aic.list[[2]] <- dredge(m2)
aic.list[[2]]$mod = "arrival"
aic.list[[2]]
summary(m2)
plot(m2)

parr <- data.frame(pred = predict(m2)) %>% 
  cbind(dat.arr)

ggplot(dat.arr) +
  geom_boxplot(aes(x=region, y=val)) +
  geom_boxplot(aes(x=region, y=pred), data=parr, fill=NA, colour="red") +
  facet_wrap(~season, scales="free")

#Duration
dat.dur <- dat %>% dplyr::filter(var=="duration")
m3 <- lmer(val ~ region*season + (1|id), data=dat.dur, na.action="na.fail", REML=FALSE)
aic.list[[3]] <- dredge(m3)
aic.list[[3]]$mod = "duration"
aic.list[[3]]
summary(m3)
plot(m3)

pdur <- data.frame(pred = predict(m3)) %>% 
  cbind(dat.dur)

ggplot(dat.dur) +
  geom_boxplot(aes(x=region, y=val)) +
  geom_boxplot(aes(x=region, y=pred), data=pdur, fill=NA, colour="red") +
  facet_wrap(~season, scales="free")

#Distance - maybe needs a different distribution (gamma?)
dat.dist <- dat %>% dplyr::filter(var=="dist", !is.na(val))
m4 <- lmer(val ~ region*season + (1|id), data=dat.dist, na.action="na.fail", REML=FALSE)
aic.list[[4]] <- dredge(m4)
aic.list[[4]]$mod = "distance"
aic.list[[4]]
m4 <- lmer(val ~ region + (1|id), data=dat.dist, na.action="na.fail", REML=FALSE)
summary(m4)
plot(m4)

pdist <- data.frame(pred = predict(m4)) %>% 
  cbind(dat.dist)

ggplot(dat.dist) +
  geom_boxplot(aes(x=region, y=log(val))) +
  geom_boxplot(aes(x=region, y=log(pred)), data=pdist, fill=NA, colour="red") +
  facet_wrap(~season, scales="free")

#Rate - maybe needs a different distribution (gamma?)
dat.rate <- dat %>% dplyr::filter(var=="rate", !is.na(val))
m5 <- lmer(val ~ region*season + (1|id), data=dat.rate, na.action="na.fail", REML=FALSE)
aic.list[[5]] <- dredge(m5)
aic.list[[5]]$mod = "rate"
aic.list[[5]]
m5 <- lmer(val ~ season + (1|id), data=dat.rate, na.action="na.fail", REML=FALSE)
summary(m5)
plot(m5)

prate <- data.frame(pred = predict(m5)) %>% 
  cbind(dat.rate)

ggplot(dat.rate) +
  geom_boxplot(aes(x=region, y=val)) +
  geom_boxplot(aes(x=region, y=pred), data=prate, fill=NA, colour="red") +
  facet_wrap(~season, scales="free")

#Stopovers
dat.stop <- dat %>% dplyr::filter(var=="Stopovers")
m6 <- glmer(val ~ region*season + (1|id), data=dat.stop, na.action="na.fail", family="poisson")
aic.list[[6]] <- dredge(m6)
aic.list[[6]]$mod = "stopover"
aic.list[[6]]
m6 <- glmer(val ~ region + season + (1|id), data=dat.stop, na.action="na.fail", family="poisson")
summary(m6)
plot(m6)

pstop <- data.frame(pred = predict(m6, type="response")) %>% 
  cbind(dat.stop)

ggplot(dat.stop) +
  geom_boxplot(aes(x=region, y=val)) +
  geom_boxplot(aes(x=region, y=pred), data=pstop, fill=NA, colour="red") +
  facet_wrap(~season, scales="free")

#Home range area
dat.area <- dat %>% dplyr::filter(var=="HRarea")
m7 <- lmer(log(val) ~ region*season + (1|id), data=dat.area, na.action="na.fail", REML = FALSE)
aic.list[[7]] <- dredge(m7)
aic.list[[7]]$mod = "area"
aic.list[[7]]
plot(m7)
summary(m7)

phr <- data.frame(pred = predict(m7)) %>% 
  cbind(dat.area)

ggplot(dat.area) +
  geom_boxplot(aes(x=region, y=log(val))) +
  geom_boxplot(aes(x=region, y=pred), data=phr, fill=NA, colour="red") +
  facet_wrap(~season, scales="free")

#Number of home ranges
dat.hr <- dat %>% dplyr::filter(var=="WinterHRs")
m8 <- glmer(val ~ region + (1|id), data=dat.hr, na.action="na.fail", family="poisson")
aic.list[[8]] <- dredge(m8)
aic.list[[8]]$mod = "winterhrs"
aic.list[[8]]
plot(m8)
summary(m8)

#Stopover duration
dat.stop <- dat %>% dplyr::filter(var=="stopoverduration")
m9 <- lmer(val ~ region*season + (1|id), data=dat.stop, na.action="na.fail", REML=FALSE)
aic.list[[9]] <- dredge(m9)
aic.list[[9]]$mod = "stopoverduration"
aic.list[[9]]
m9 <- lmer(val ~ season + (1|id), data=dat.stop, na.action="na.fail", REML=FALSE)
plot(m9)
summary(m9)

pstop <- data.frame(pred = predict(m9, type="response")) %>% 
  cbind(dat.stop)

ggplot(dat.stop) +
  geom_boxplot(aes(x=region, y=val)) +
  geom_boxplot(aes(x=region, y=pred), data=pstop, fill=NA, colour="red") +
  facet_wrap(~season, scales="free")

#11. Collapse AIC results----
aic.out <- rbindlist(aic.list, fill=TRUE) %>% 
  arrange(mod, -weight) %>% 
  data.frame() %>% 
  mutate(weight = round(weight, 2),
         logLik = round(logLik, 2),
         AICc = round(AICc, 2),
         delta = round(delta, 2),
         model = case_when(is.na(region.season) & region=="+" & season=="+" ~ "group + season",
                           !is.na(region.season) ~ "group * season",
                           is.na(region.season) & is.na(region) & season=="+" ~ "season",
                           is.na(region.season) & is.na(season) & region=="+" ~ "group",
                           is.na(region.season) & is.na(season) & is.na(region) ~ "1")) %>% 
  dplyr::select(mod, model, df, logLik, AICc, delta, weight)

#12. Save----
write.csv(aic.out, "Results/MigrationCharsAIC.csv", row.names = FALSE)
