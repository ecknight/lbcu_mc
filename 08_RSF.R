library(tidyverse)
library(usdm)
library(lme4)
library(MuMIn)
library(ggridges)
library(mgcv)

options(scipen=99999)

#1. Load data----
dat <- read.csv("Data/LBCU_environvars_RSF.csv") %>% 
  mutate(response = ifelse(type=="used", 1, 0),
         change.raw = change,
         drought.raw = drought,
         seasonality.raw = seasonality) %>% 
  mutate(change = (change.raw - min(change.raw))/(max(change.raw) - min(change.raw)),
         drought = (drought.raw - min(drought.raw))/(max(drought.raw) - min(drought.raw)),
         seasonality = (seasonality.raw - min(seasonality.raw))/(max(seasonality.raw) - min(seasonality.raw))) %>% 
  mutate(Region = case_when(group==1 ~ "Central",
                            group==2 ~ "West",
                            group==3 ~ "East"),
         Region = factor(Region, levels=c("Central", "West", "East")))

#2. Test VIF----
loop <- dat %>% 
  dplyr::select(season) %>% 
  unique()

v.list <- list()
for(i in 1:nrow(loop)){
  
  cov.i <- dat %>% 
    dplyr::filter(season==loop$season[i]) %>% 
    dplyr::select(change, crop, grass, seasonality)
  
  v.list[[i]] <- vif(cov.i) %>% 
    mutate(season=loop$season[i])
  
}

v <- do.call(rbind, v.list) %>% 
  arrange(-VIF)
#all good

#3. Visualize----
#Check for sufficient distribution
#Check for polynomials

#wetland change
ggplot(dat) +
  geom_histogram(aes(x=change, colour=Region)) +
  facet_wrap(Region~season, scales="free")

ggplot(dat) +
#  geom_jitter(aes(x=change, y=response, colour=Region)) +
  geom_smooth(aes(x=change, y=response, colour=Region)) +
  facet_wrap(~season)
#maybe for winter

#crop
ggplot(dat) +
  geom_histogram(aes(x=crop, colour=Region)) +
  facet_wrap(Region~season, scales="free")

ggplot(dat) +
#  geom_jitter(aes(x=crop, y=response, colour=Region)) +
  geom_smooth(aes(x=crop, y=response, colour=Region)) +
  facet_wrap(~season)

#grass
ggplot(dat) +
  geom_histogram(aes(x=grass, colour=Region)) +
  facet_wrap(Region~season, scales="free")

ggplot(dat) +
#  geom_jitter(aes(x=grass, y=response, colour=Region)) +
  geom_smooth(aes(x=grass, y=response, colour=Region)) +
  facet_wrap(~season)

#seasonality
ggplot(dat) +
  geom_histogram(aes(x=seasonality, colour=Region)) +
  facet_wrap(Region~season, scales="free")

ggplot(dat) +
#  geom_jitter(aes(x=seasonality, y=response, colour=Region)) +
  geom_smooth(aes(x=seasonality, y=response, colour=Region)) +
  facet_wrap(~season)

#4. Set up loop----
loop <- dat %>% 
  dplyr::select(season) %>% 
  unique()

m.list <- list()
s.list <- list()
d.list <- list()
for(i in 1:nrow(loop)){
  
  #5. Subset data----
  dat.i <- dat %>% 
    dplyr::filter(season==loop$season[i])
  
  #6. Model----
  m.list[[i]] <- glm(response ~ crop*Region + grass*Region + seasonality*Region + change*Region, family = "binomial", data=dat.i, na.action = "na.fail")
  
  #7. Save summary----
  s.list[[i]] <- summary(m.list[[i]])$coefficients %>% 
    data.frame() %>% 
    mutate(season = loop$season[i],
           cov = row.names(summary(m.list[[i]])$coefficients)) 
  
  #8. Dredge----
  d.list[[i]] <- data.frame(dredge(m.list[[i]])) %>% 
    mutate(season = loop$season[i])
  
}

#9. Wrangle output----
summary <- do.call(rbind, s.list) %>% 
  dplyr::filter(cov!="(Intercept)") %>% 
  rename(p = 'Pr...z..', se = 'Std..Error') %>% 
  mutate(sig = ifelse(p < 0.05, 1, 0))

dredged <- do.call(rbind, d.list) %>% 
  dplyr::filter(delta < 2) %>%
  group_by(season) %>% 
  mutate(mindf = min(df)) %>% 
  mutate(pick = ifelse(df==mindf, 1, 0)) %>% 
  dplyr::filter(pick==1) %>%
  arrange(df) %>% 
  dplyr::filter(row_number()==1) %>% 
  ungroup() %>% 
  arrange(season)

#10. Run final models----

m.final <- list()

#breed
dat.breed <- dat %>% 
  dplyr::filter(season=="breed")

m.final[[1]] <- glm(response ~ crop + grass*Region + change*Region + seasonality*Region, family = "binomial", data=dat.breed, na.action = "na.fail")

#fall
dat.fall <- dat %>% 
  dplyr::filter(season=="fallmig")

m.final[[2]]  <- glm(response ~ crop*Region + grass, family = "binomial", data=dat.fall, na.action = "na.fail")

#winter
dat.winter <- dat %>% 
  dplyr::filter(season=="winter")

m.final[[3]] <- glm(response ~ crop + grass*Region + change*Region + seasonality*Region, family = "binomial", data=dat.winter, na.action = "na.fail")

#spring
dat.spring <- dat %>% 
  dplyr::filter(season=="springmig")

m.final[[4]]  <- glm(response ~ crop*Region + seasonality*Region, family = "binomial", data=dat.spring, na.action = "na.fail")

#11. Predictions----
loop <- unique(dat$season)

pred.list <- list()
for(i in 1:length(loop)){
  
  dat.i <- dat %>% 
    dplyr::filter(season==loop[i])
  
  newdat.crop <- expand.grid(Region=unique(dat.spring$Region), crop=seq(0, 1, 0.01), grass=0, seasonality=0, change=0)
  
  newdat.grass <- expand.grid(Region=unique(dat.spring$Region), crop=0, grass=seq(0, 1, 0.01), seasonality=0, change=0)
  
  newdat.seasonality <- expand.grid(Region=unique(dat.spring$Region), crop=0, grass=0, seasonality=seq(0, 1, 0.01), change=0)
  
  newdat.change = expand.grid(Region=unique(dat.spring$Region), crop=0, grass=0, seasonality=0, change=seq(0, 1, 0.01))
  
  pred.crop <- predict(m.final[[i]], newdat.crop, se=TRUE, type="response")  %>% 
    data.frame() %>% 
    cbind(newdat.crop %>% 
            dplyr::select(Region, crop)) %>% 
    rename(val=crop) %>% 
    mutate(up=fit+1.96*se.fit, low=fit-1.96*se.fit,
           cov="crop")
  
  pred.grass <- predict(m.final[[i]], newdat.grass, se=TRUE, type="response")  %>% 
    data.frame() %>% 
    cbind(newdat.grass %>% 
            dplyr::select(Region, grass)) %>% 
    rename(val=grass) %>% 
    mutate(up=fit+1.96*se.fit, low=fit-1.96*se.fit,
           cov="grass")
  
  pred.seasonality <- predict(m.final[[i]], newdat.seasonality, se=TRUE, type="response")  %>% 
    data.frame() %>% 
    cbind(newdat.seasonality %>% 
            dplyr::select(Region, seasonality)) %>% 
    rename(val=seasonality) %>% 
    mutate(up=fit+1.96*se.fit, low=fit-1.96*se.fit,
           cov="seasonality")
  
  pred.change <- predict(m.final[[i]], newdat.change, se=TRUE, type="response")  %>% 
    data.frame() %>% 
    cbind(newdat.change %>% 
            dplyr::select(Region, change)) %>% 
    rename(val=change) %>% 
    mutate(up=fit+1.96*se.fit, low=fit-1.96*se.fit,
           cov="change")
  
  pred.list[[i]] <- rbind(pred.crop, pred.grass, pred.seasonality, pred.change) %>% 
    mutate(season=loop[i])
  
}

pred <- do.call(rbind, pred.list) %>% 
  dplyr::filter(!(season=="fallmig" & cov %in% c("change", "seasonality")),
                !(season=="springmig" & cov %in% c("grass", "change")))

write.csv(pred, "Results/RSFPredictions.csv", row.names = FALSE)


#12. Plot----
ggplot(pred) +
  geom_ribbon(aes(x=val, ymin=low, ymax=up, group=Region), alpha=0.3) +
  geom_line(aes(x=val, y=fit, colour=Region)) +
  facet_grid(cov~season, scales="free")