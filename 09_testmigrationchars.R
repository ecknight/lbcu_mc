library(tidyverse)
library(vegan)
library(lme4)
library(MuMIn)
library(sf)
library(adehabitatLT)
library(data.table)

options(scipen=99999)

setwd("G:/My Drive/SMBC")

#1. Get migration characteristic data----
dat <- read.csv("Data/MovementBehaviours.csv") |> 
  mutate(group = as.factor(group))

#2. Visualize----
#Departure
ggplot(dat %>% dplyr::filter(var=="depart")) +
  geom_boxplot(aes(x=group, y=val)) +
  facet_wrap(nclust~season, scales="free")

#Arrival
ggplot(dat %>% dplyr::filter(var=="arrive")) +
  geom_boxplot(aes(x=group, y=val)) +
  facet_wrap(nclust~season, scales="free")

#Rate
ggplot(dat %>% dplyr::filter(var=="rate")) +
  geom_boxplot(aes(x=group, y=log(val))) +
  facet_wrap(nclust~season, scales="free")

#Number of stopovers
ggplot(dat %>% dplyr::filter(var=="Stopovers")) +
  geom_histogram(aes(x=val, fill=group)) +
  facet_grid(nclust~season)

#Area
ggplot(dat %>% dplyr::filter(var=="HRarea")) +
  geom_boxplot(aes(x=group, y=log(val))) +
  facet_wrap(nclust~season,scales="free")

#Wintering home ranges
ggplot(dat %>% dplyr::filter(var=="WinterHRs")) +
  geom_histogram(aes(x=val, fill=group)) +
  facet_wrap(nclust~season)

#3. Set up loop----
clusts <- unique(dat$nclust)

out.list <- list()
for(i in 1:length(clusts)){
  
  #3. Test----
  aic.list <- list()
  
  #Departure
  dat.dep <- dat %>% dplyr::filter(var=="depart2", nclust==clusts[i])
  m1 <- lmer(val ~ group*season + (1|id), data=dat.dep, na.action="na.fail", REML=FALSE)
  aic.list[[1]] <- dredge(m1)
  aic.list[[1]]$mod = "departure"
  aic.list[[1]]
  summary(m1)
  plot(m1)
  pdep <- data.frame(pred = predict(m1)) %>% 
    cbind(dat.dep)
  
  ggplot(dat.dep) +
    geom_boxplot(aes(x=group, y=val)) +
    geom_boxplot(aes(x=group, y=pred), data=pdep, fill=NA, colour="red") +
    facet_wrap(~season, scales="free")
  
  #Arrival
  dat.arr <- dat %>% dplyr::filter(var=="arrive2", nclust==clusts[i])
  m2 <- lmer(val ~ group*season + (1|id), data=dat.arr, na.action="na.fail", REML=FALSE)
  aic.list[[2]] <- dredge(m2)
  aic.list[[2]]$mod = "arrival"
  aic.list[[2]]
  summary(m2)
  plot(m2)
  
  parr <- data.frame(pred = predict(m2)) %>% 
    cbind(dat.arr)
  
  ggplot(dat.arr) +
    geom_boxplot(aes(x=group, y=val)) +
    geom_boxplot(aes(x=group, y=pred), data=parr, fill=NA, colour="red") +
    facet_wrap(~season, scales="free")

  #Rate - maybe needs a different distribution (gamma?)
  dat.rate <- dat %>% dplyr::filter(var=="rate", !is.na(val), nclust==clusts[i])
  m5 <- lmer(val ~ group*season + (1|id), data=dat.rate, na.action="na.fail", REML=FALSE)
  aic.list[[5]] <- dredge(m5)
  aic.list[[5]]$mod = "rate"
  aic.list[[5]]
  m5 <- lmer(val ~ season + (1|id), data=dat.rate, na.action="na.fail", REML=FALSE)
  summary(m5)
  plot(m5)

  prate <- data.frame(pred = predict(m5)) %>%
    cbind(dat.rate)

  ggplot(dat.rate) +
    geom_boxplot(aes(x=group, y=val)) +
    geom_boxplot(aes(x=group, y=pred), data=prate, fill=NA, colour="red") +
    facet_wrap(~season, scales="free")
  
  #Stopovers
  dat.stop <- dat %>% dplyr::filter(var=="Stopovers", nclust==clusts[i])
  m6 <- glmer(val ~ group*season + (1|id), data=dat.stop, na.action="na.fail", family="poisson")
  aic.list[[6]] <- dredge(m6)
  aic.list[[6]]$mod = "stopover"
  aic.list[[6]]
  m6 <- glmer(val ~ group + season + (1|id), data=dat.stop, na.action="na.fail", family="poisson")
  summary(m6)
  plot(m6)
  
  pstop <- data.frame(pred = predict(m6, type="response")) %>% 
    cbind(dat.stop)
  
  ggplot(dat.stop) +
    geom_boxplot(aes(x=group, y=val)) +
    geom_boxplot(aes(x=group, y=pred), data=pstop, fill=NA, colour="red") +
    facet_wrap(~season, scales="free")
  
  #Home range area
  dat.area <- dat %>% dplyr::filter(var=="HRarea", nclust==clusts[i])
  m7 <- lmer(log(val) ~ group*season + (1|id), data=dat.area, na.action="na.fail", REML = FALSE)
  aic.list[[7]] <- dredge(m7)
  aic.list[[7]]$mod = "area"
  aic.list[[7]]
  plot(m7)
  summary(m7)
  
  phr <- data.frame(pred = predict(m7)) %>% 
    cbind(dat.area)
  
  ggplot(dat.area) +
    geom_boxplot(aes(x=group, y=log(val))) +
    geom_boxplot(aes(x=group, y=pred), data=phr, fill=NA, colour="red") +
    facet_wrap(~season, scales="free")
  
  #Number of home ranges
  dat.hr <- dat %>% dplyr::filter(var=="WinterHRs", nclust==clusts[i])
  m8 <- glmer(val ~ group + (1|id), data=dat.hr, na.action="na.fail", family="poisson")
  aic.list[[8]] <- dredge(m8)
  aic.list[[8]]$mod = "winterhrs"
  aic.list[[8]]
  plot(m8)
  summary(m8)
  
  #Stopover duration
  dat.stop <- dat %>% dplyr::filter(var=="stopoverduration", nclust==clusts[i])
  m9 <- lmer(val ~ group*season + (1|id), data=dat.stop, na.action="na.fail", REML=FALSE)
  aic.list[[9]] <- dredge(m9)
  aic.list[[9]]$mod = "stopoverduration"
  aic.list[[9]]
  m9 <- lmer(val ~ season + (1|id), data=dat.stop, na.action="na.fail", REML=FALSE)
  plot(m9)
  summary(m9)
  
  pstop <- data.frame(pred = predict(m9, type="response")) %>% 
    cbind(dat.stop)
  
  ggplot(dat.stop) +
    geom_boxplot(aes(x=group, y=val)) +
    geom_boxplot(aes(x=group, y=pred), data=pstop, fill=NA, colour="red") +
    facet_wrap(~season, scales="free")
  
  #11. Collapse AIC results----
  out.list[[i]] <- rbindlist(aic.list, fill=TRUE) %>% 
    arrange(mod, -weight) %>% 
    data.frame() %>% 
    mutate(weight = round(weight, 2),
           logLik = round(logLik, 2),
           AICc = round(AICc, 2),
           delta = round(delta, 2),
           nclust = clusts[i],
           model = case_when(is.na(group.season) & group=="+" & season=="+" ~ "group + season",
                             !is.na(group.season) ~ "group * season",
                             is.na(group.season) & is.na(group) & season=="+" ~ "season",
                             is.na(group.season) & is.na(season) & group=="+" ~ "group",
                             is.na(group.season) & is.na(season) & is.na(group) ~ "1")) %>% 
    dplyr::select(nclust, mod, model, df, logLik, AICc, delta, weight)
  
  
}

aic.out <- do.call(rbind, out.list) |> 
  arrange(mod, nclust)

#12. Save----
write.csv(aic.out, "Results/MigrationCharsAIC.csv", row.names = FALSE)
