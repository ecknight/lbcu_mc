library(tidyverse)
library(usdm)
library(lme4)
library(MuMIn)

options(scipen=99999)

#1. Load data----
dat <- read.csv("Data/LBCU_environvars_RSF.csv") %>% 
  mutate(response = ifelse(type=="used", 1, 0)) %>% 
  mutate(change = (change - min(change))/(max(change) - min(change)),
         drought = (drought - min(drought))/(max(drought) - min(drought)),
         seasonality = (seasonality - min(seasonality))/(max(seasonality) - min(seasonality)))

#2. Test VIF----
cov <- dat %>% 
  dplyr::select(built, change, crop, drought, grass, seasonality, wetland)

vif(cov)
#all good! amazing

#3. Set up loop----
loop <- dat %>% 
  dplyr::select(group, season) %>% 
  unique()

m.list <- list()
s.list <- list()
d.list <- list()
for(i in 1:nrow(loop)){
  
  #4. Subset data----
  dat.i <- dat %>% 
    dplyr::filter(group==loop$group[i],
                  season==loop$season[i])
  #5. Model----
  m.list[[i]] <- glm(response ~ built + crop + drought + grass + seasonality + wetland, family = "binomial", data=dat.i, na.action = "na.fail")
  
  #6. Save summary----
  s.list[[i]] <- summary(m.list[[i]])$coefficients %>% 
    data.frame() %>% 
    mutate(group = loop$group[i],
           season = loop$season[i],
           cov = row.names(summary(m.list[[i]])$coefficients)) 
  
  #7. Dredge----
  d.list[[i]] <- dredge(m.list[[i]]) %>% 
    data.frame() %>% 
    mutate(group = loop$group[i],
           season = loop$season[i])
  
  
}

#8. Wrangle output----
summary <- do.call(rbind, s.list) %>% 
  dplyr::filter(cov!="(Intercept)") %>% 
  rename(p = 'Pr...z..', se = 'Std..Error') %>% 
  mutate(sig = ifelse(p < 0.05, 1, 0)) %>% 
  mutate(Region = case_when(group==1 ~ "Central",
                            group==2 ~ "West",
                            group==3 ~ "East"))

dredged <- do.call(rbind, d.list)

#9. Visualize----
ggplot(data %>% dplyr::filter(abs(Estimate) < 20)) + 
  geom_errorbar(aes(x=season, ymin=Estimate-se, ymax=Estimate+se, colour=Region, alpha=sig),
                position="dodge", width=1) +
  geom_point(aes(x=season, y=Estimate, colour=Region), position = position_dodge(width=1)) +
  geom_hline(aes(yintercept=0)) +
  facet_wrap(~cov, scales="free_y")
