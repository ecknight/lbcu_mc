library(tidyverse)
library(usdm)
library(lme4)
library(MuMIn)
library(ggridges)

options(scipen=99999)

#1. Load data----
dat <- read.csv("Data/LBCU_environvars_RSF.csv") %>% 
  mutate(response = ifelse(type=="used", 1, 0)) %>% 
  mutate(change = (change - min(change))/(max(change) - min(change)),
         drought = (drought - min(drought))/(max(drought) - min(drought)),
         seasonality = (seasonality - min(seasonality))/(max(seasonality) - min(seasonality))) %>% 
  mutate(Region = case_when(group==1 ~ "Central",
                            group==2 ~ "West",
                            group==3 ~ "East"))

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
  dplyr::select(Region, season) %>% 
  unique()

m.list <- list()
s.list <- list()
d.list <- list()
for(i in 1:nrow(loop)){
  
  #5. Subset data----
  dat.i <- dat %>% 
    dplyr::filter(Region==loop$Region[i],
                  season==loop$season[i]) %>% 
    mutate(change2=change^2)
  
  #6. Model----
  m.list[[i]] <- glmer(response ~ crop + grass + seasonality + change + drought + wetland + (1|bird), family = "binomial", data=dat.i, na.action = "na.fail")
  
  #7. Save summary----
  s.list[[i]] <- summary(m.list[[i]])$coefficients %>% 
    data.frame() %>% 
    mutate(Region = loop$Region[i],
           season = loop$season[i],
           cov = row.names(summary(m.list[[i]])$coefficients)) 
  
  #8. Dredge----
  d.list[[i]] <- dredge(m.list[[i]]) %>% 
    data.frame() %>% 
    mutate(Region = loop$Region[i],
           season = loop$season[i])
  
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

#10. Visualize used & available----
dat.long <- dat %>% 
  dplyr::select(-built) %>% 
  pivot_longer(change:wetland, names_to="cov", values_to="val")

ggplot(dat.long %>% dplyr::filter(type=="avail")) +
  geom_density_ridges(aes(x=val, fill=factor(Region), y=season), alpha = 0.5) +
  facet_wrap(~cov, scales="free")

#11. Test for differences in availability----
dat.avail <- dat %>% 
  dplyr::filter(response==0)

ma <- lm(seasonality ~ Region*season, data=dat.avail, na.action="na.fail")
dredge(ma)
plot(ma)

#12. Test for functional response across groups----

#5. Subset data----
dat.i <- dat %>% 
  dplyr::filter(Region==loop$Region[i],
                season==loop$season[i]) %>% 
  mutate(change2=change^2)

#6. Model----
m.list[[i]] <- glmer(response ~ crop + grass + seasonality + change + drought +  (1|bird), family = "binomial", data=dat.i, na.action = "na.fail")
