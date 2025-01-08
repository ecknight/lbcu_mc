library(tidyverse)
library(lme4)
library(MuMIn)

options(scipen=99999)

#1. Get migration characteristic data----
dat <- read.csv("Data/MovementBehaviours.csv") |> 
  mutate(group = as.factor(group))

#2. Visualize----
#Departure
ggplot(dat |> dplyr::filter(var=="depart")) +
  geom_boxplot(aes(x=group, y=val)) +
  facet_wrap(nclust~season, scales="free")

#Arrival
ggplot(dat |> dplyr::filter(var=="arrive")) +
  geom_boxplot(aes(x=group, y=val)) +
  facet_wrap(nclust~season, scales="free")

#Rate
ggplot(dat |> dplyr::filter(var=="rate")) +
  geom_boxplot(aes(x=group, y=log(val))) +
  facet_wrap(nclust~season, scales="free")

#Number of stopovers
ggplot(dat |> dplyr::filter(var=="Stopovers")) +
  geom_histogram(aes(x=val, fill=group)) +
  facet_grid(nclust~season)

#Area
ggplot(dat |> dplyr::filter(var=="HRarea")) +
  geom_boxplot(aes(x=group, y=log(val))) +
  facet_wrap(nclust~season,scales="free")

#Wintering home ranges
ggplot(dat |> dplyr::filter(var=="WinterHRs")) +
  geom_histogram(aes(x=val, fill=group)) +
  facet_wrap(nclust~season)

#3. Set up loop----
clusts <- unique(dat$nclust)

out.list <- list()
for(i in 1:length(clusts)){
  
  #3. Test----
  aic.list <- list()
  
  #Departure
  dat.dep <- dat |> dplyr::filter(var=="depart2", nclust==clusts[i])
  m1 <- lmer(val ~ group*season + (1|id), data=dat.dep, na.action="na.fail", REML=FALSE)
  aic.list[[1]] <- dredge(m1)
  aic.list[[1]]$mod = "departure"
  
  #Arrival
  dat.arr <- dat |> dplyr::filter(var=="arrive2", nclust==clusts[i])
  m2 <- lmer(val ~ group*season + (1|id), data=dat.arr, na.action="na.fail", REML=FALSE)
  aic.list[[2]] <- dredge(m2)
  aic.list[[2]]$mod = "arrival"

  #Rate - maybe needs a different distribution (gamma?)
  dat.rate <- dat |> dplyr::filter(var=="rate", !is.na(val), nclust==clusts[i])
  m5 <- lmer(val ~ group*season + (1|id), data=dat.rate, na.action="na.fail", REML=FALSE)
  aic.list[[5]] <- dredge(m5)
  aic.list[[5]]$mod = "rate"
  
  #Stopovers
  dat.stop <- dat |> dplyr::filter(var=="Stopovers", nclust==clusts[i])
  m6 <- glmer(val ~ group*season + (1|id), data=dat.stop, na.action="na.fail", family="poisson")
  aic.list[[6]] <- dredge(m6)
  aic.list[[6]]$mod = "stopover"
  
  #Home range area
  dat.area <- dat |> dplyr::filter(var=="HRarea", nclust==clusts[i])
  m7 <- lmer(log(val) ~ group*season + (1|id), data=dat.area, na.action="na.fail", REML = FALSE)
  aic.list[[7]] <- dredge(m7)
  aic.list[[7]]$mod = "area"
  
  #Number of home ranges
  dat.hr <- dat |> dplyr::filter(var=="WinterHRs", nclust==clusts[i])
  m8 <- glmer(val ~ group + (1|id), data=dat.hr, na.action="na.fail", family="poisson")
  aic.list[[8]] <- dredge(m8)
  aic.list[[8]]$mod = "winterhrs"
  
  #Stopover duration
  dat.stop <- dat |> dplyr::filter(var=="stopoverduration", nclust==clusts[i])
  m9 <- lmer(val ~ group*season + (1|id), data=dat.stop, na.action="na.fail", REML=FALSE)
  aic.list[[9]] <- dredge(m9)
  aic.list[[9]]$mod = "stopoverduration"
  
  #11. Collapse AIC results----
  out.list[[i]] <- rbindlist(aic.list, fill=TRUE) |> 
    arrange(mod, -weight) |> 
    data.frame() |> 
    mutate(weight = round(weight, 2),
           logLik = round(logLik, 2),
           AICc = round(AICc, 2),
           delta = round(delta, 2),
           nclust = clusts[i],
           model = case_when(is.na(group.season) & group=="+" & season=="+" ~ "group + season",
                             !is.na(group.season) ~ "group * season",
                             is.na(group.season) & is.na(group) & season=="+" ~ "season",
                             is.na(group.season) & is.na(season) & group=="+" ~ "group",
                             is.na(group.season) & is.na(season) & is.na(group) ~ "1")) |> 
    dplyr::select(nclust, mod, model, df, logLik, AICc, delta, weight)
  
  
}

aic.out <- do.call(rbind, out.list) |> 
  arrange(mod, nclust)

#12. Save----
write.csv(aic.out, "Results/MigrationCharsAIC.csv", row.names = FALSE)
