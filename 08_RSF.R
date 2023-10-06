library(tidyverse)
library(usdm)
library(lme4)
library(MuMIn)
library(ggridges)
library(mgcv)

options(scipen=99999)

#1. Load data----
dat <- read.csv("Data/LBCU_environvars_RSF.csv") %>% 
  mutate(response = ifelse(type=="used", 1, 0)) %>% 
  mutate(change = (change - min(change))/(max(change) - min(change)),
         drought = (drought - min(drought))/(max(drought) - min(drought)),
         seasonality = (seasonality - min(seasonality))/(max(seasonality) - min(seasonality))) %>% 
  mutate(Region = case_when(group==1 ~ "Central",
                            group==2 ~ "West",
                            group==3 ~ "East"),
         Region = factor(Region, levels=c("Central", "West", "East")))

#2. Test VIF----
loop <- dat %>% 
  dplyr::select(Region, season) %>% 
  unique()

v.list <- list()
for(i in 1:nrow(loop)){
  
  cov.i <- dat %>% 
    dplyr::filter(Region==loop$Region[i],
                  season==loop$season[i]) %>% 
    dplyr::select(change, crop, drought, grass, seasonality)
  
  v.list[[i]] <- vif(cov.i) %>% 
    mutate(season=loop$season[i],
           Region=loop$Region[i])
  
}

v <- do.call(rbind, v.list) %>% 
  arrange(-VIF)
#spring mig for east population is only issue

cor(dat %>% 
      dplyr::filter(Region=="East",
                    season=="springmig") %>% 
      dplyr::select(built, change, crop, drought, grass, seasonality, wetland))
#crop and grass is -0.61. Let's leave it


#3. Visualize----
#Check for sufficient distribution
#Check for polynomials
ggplot(dat) +
  geom_histogram(aes(x=built, colour=Region)) +
  facet_wrap(Region~season, scales="free")

ggplot(dat) +
  geom_jitter(aes(x=built, y=response, colour=Region)) +
  geom_smooth(aes(x=built, y=response, colour=Region)) +
  facet_wrap(~season)

ggplot(dat) +
  geom_histogram(aes(x=change, colour=Region)) +
  facet_wrap(Region~season, scales="free")

ggplot(dat) +
#  geom_jitter(aes(x=change, y=response, colour=Region)) +
  geom_smooth(aes(x=change, y=response, colour=Region)) +
  facet_wrap(~season)
#maybe for winter

ggplot(dat) +
  geom_histogram(aes(x=crop, colour=Region)) +
  facet_wrap(Region~season, scales="free")

ggplot(dat) +
#  geom_jitter(aes(x=crop, y=response, colour=Region)) +
  geom_smooth(aes(x=crop, y=response, colour=Region)) +
  facet_wrap(~season)

ggplot(dat) +
  geom_histogram(aes(x=drought, colour=Region)) +
  facet_wrap(Region~season, scales="free")

ggplot(dat) +
#  geom_jitter(aes(x=drought, y=response, colour=Region)) +
  geom_smooth(aes(x=drought, y=response, colour=Region)) +
  facet_wrap(~season)

ggplot(dat) +
  geom_histogram(aes(x=grass, colour=Region)) +
  facet_wrap(Region~season, scales="free")

ggplot(dat) +
#  geom_jitter(aes(x=grass, y=response, colour=Region)) +
  geom_smooth(aes(x=grass, y=response, colour=Region)) +
  facet_wrap(~season)

ggplot(dat) +
  geom_histogram(aes(x=seasonality, colour=Region)) +
  facet_wrap(Region~season, scales="free")

ggplot(dat) +
#  geom_jitter(aes(x=seasonality, y=response, colour=Region)) +
  geom_smooth(aes(x=seasonality, y=response, colour=Region)) +
  facet_wrap(~season)

ggplot(dat) +
  geom_histogram(aes(x=wetland, colour=Region)) +
  facet_wrap(Region~season, scales="free")

ggplot(dat) +
#  geom_jitter(aes(x=wetland, y=response, colour=Region)) +
  geom_smooth(aes(x=wetland, y=response, colour=Region)) +
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
  group_by(season, Region) %>% 
  mutate(mindf = min(df)) %>% 
  ungroup() %>% 
  mutate(pick = ifelse(df==mindf, 1, 0)) %>% 
  dplyr::filter(pick==1) %>% 
  arrange(season, Region)

#10. Visualize selection----
#By pvalue
ggplot(summary %>% dplyr::filter(sig==1)) + 
  geom_errorbar(aes(x=season, ymin=Estimate-se, ymax=Estimate+se, colour=Region),
                position="dodge", width=1) +
  #  geom_point(aes(x=season, y=Estimate, colour=Region), position = position_dodge(width=1)) +
  geom_hline(aes(yintercept=0)) +
  facet_wrap(~cov, scales="free_y")

#By AIC
dredged.plot <- dredged %>% 
  pivot_longer(change:seasonality, names_to="cov", values_to = "intercept") %>% 
  dplyr::filter(!is.na(intercept)) %>% 
  left_join(summary)

ggplot(dredged.plot) + 
  geom_errorbar(aes(x=season, ymin=Estimate-se, ymax=Estimate+se, colour=Region),
                position="dodge", width=1) +
  #  geom_point(aes(x=season, y=Estimate, colour=Region), position = position_dodge(width=1)) +
  geom_hline(aes(yintercept=0)) +
  facet_wrap(~cov, scales="free_y")

#11. Visualize used & available----
dat.long <- dat %>% 
  dplyr::select(-built) %>% 
  pivot_longer(change:wetland, names_to="cov", values_to="val")

ggplot(dat.long %>% dplyr::filter(type=="avail")) +
  geom_density_ridges(aes(x=val, fill=factor(Region), y=season), alpha = 0.5) +
  facet_wrap(~cov, scales="free")

#12. Test for differences in availability----
dat.avail <- dat %>% 
  dplyr::filter(response==0)

ma <- lm(seasonality ~ Region*season, data=dat.avail, na.action="na.fail")
dredge(ma)

ma <- lm(crop ~ Region*season, data=dat.avail, na.action="na.fail")
dredge(ma)

ma <- lm(change ~ Region*season, data=dat.avail, na.action="na.fail")
dredge(ma)

ma <- lm(drought ~ Region*season, data=dat.avail, na.action="na.fail")
dredge(ma)

ma <- lm(grass ~ Region*season, data=dat.avail, na.action="na.fail")
dredge(ma)

#13. Test for functional response across groups----

#13a. Crop----
sum.crop <- summary %>% 
  dplyr::filter(cov=="crop") %>% 
  left_join(dat %>% 
              group_by(season, Region) %>% 
              summarize(avail = mean(crop)) %>% 
              ungroup())
m.crop <- lm(Estimate ~ log(avail), data=sum.crop)
summary(m.crop)

ggplot(sum.crop) +
  geom_point(aes(x=log(avail), y=Estimate, colour=season))

#13b. Drought----
sum.drought <- summary %>% 
  dplyr::filter(cov=="drought") %>% 
  left_join(dat %>% 
              group_by(season, Region) %>% 
              summarize(avail = mean(drought)) %>% 
              ungroup())
m.drought <- lm(Estimate ~ log(avail), data=sum.drought)
summary(m.drought)

ggplot(sum.drought) +
  geom_point(aes(x=log(avail), y=Estimate, colour=season))

#13c. Change----
sum.change <- summary %>% 
  dplyr::filter(cov=="change") %>% 
  left_join(dat %>% 
              group_by(season, Region) %>% 
              summarize(avail = mean(change)) %>% 
              ungroup())
m.change <- lm(Estimate ~ log(avail), data=sum.change)
summary(m.change)

ggplot(sum.change) +
  geom_point(aes(x=log(avail), y=Estimate, colour=season))

#13d. Grass----
sum.grass <- summary %>% 
  dplyr::filter(cov=="grass") %>% 
  left_join(dat %>% 
              group_by(season, Region) %>% 
              summarize(avail = mean(grass)) %>% 
              ungroup())
m.grass <- lm(Estimate ~ log(avail), data=sum.grass)
summary(m.grass)

ggplot(sum.grass) +
  geom_point(aes(x=log(avail), y=Estimate, colour=season, pch=Region))

#13e. Seasonality----
sum.seasonality <- summary %>% 
  dplyr::filter(cov=="seasonality") %>% 
  left_join(dat %>% 
              group_by(season, Region) %>% 
              summarize(avail = mean(seasonality)) %>% 
              ungroup())
m.seasonality <- lm(Estimate ~ log(avail), data=sum.seasonality)
summary(m.seasonality)

ggplot(sum.seasonality) +
  geom_point(aes(x=log(avail), y=Estimate, colour=season, pch=Region))

#14. Run final models----
#breed
dat.breed <- dat %>% 
  dplyr::filter(season=="breed")

m.breed <- glm(response ~ crop + grass*Region + seasonality*Region + change*Region, family = "binomial", data=dat.breed, na.action = "na.fail")
summary(m.breed)

#fall
dat.fall <- dat %>% 
  dplyr::filter(season=="fallmig")

m.fall <- glm(response ~ crop*Region + grass, family = "binomial", data=dat.fall, na.action = "na.fail")
summary(m.fall)

#winter
dat.winter <- dat %>% 
  dplyr::filter(season=="winter")

m.winter <- glm(response ~ crop + grass*Region + seasonality*Region + change*Region, family = "binomial", data=dat.winter, na.action = "na.fail")
summary(m.winter)

#spring
dat.spring <- dat %>% 
  dplyr::filter(season=="springmig")

m.spring <- glm(response ~ crop*Region + seasonality*Region, family = "binomial", data=dat.spring, na.action = "na.fail")
summary(m.spring)

#15. Predictions----
newdat <- expand.grid(Region=unique(dat$Region), crop=seq(0, 1, 0.01), grass=mean(dat.spring$grass), seasonality=mean(dat.spring$seasonality), change=mean(dat.spring$change))

pred.spring.crop <- data.frame(pred=predict(m.spring, newdat, se=TRUE)) %>% 
  rename(pred=pred.fit, se=pred.se.fit) %>% 
  dplyr::select(pred, se) %>% 
  cbind(newdat) %>% 
  mutate(up=pred+1.96*se, low=pred-1.96*se)

ggplot(pred.spring.crop) +
  geom_ribbon(aes(x=crop, ymin=low, ymax=up, group=Region), alpha=0.2) +
  geom_line(aes(x=crop, y=pred, colour=Region))
