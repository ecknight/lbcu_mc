library(tidyverse)
library(vegan)
library(lme4)
library(MuMIn)
library(rstatix)

options(scipen=99999)

#1. Read in data-----
dat <- read.csv("Data/LBCU_environvars.csv") %>% 
  mutate(region = case_when(kdecluster==1 ~ "east", 
                            kdecluster==2 ~ "west"),
         conv = covcrop + covbuilt)

#2. Visualize----
ggplot(dat) +
  geom_boxplot(aes(x=season, y=log(built), fill=region))

ggplot(dat) +
  geom_boxplot(aes(x=season, y=crops, fill=region))

ggplot(dat) +
  geom_boxplot(aes(x=season, y=grass, fill=region))

ggplot(dat) +
  geom_boxplot(aes(x=season, y=log(flooded_vegetation), fill=region))

ggplot(dat) +
  geom_boxplot(aes(x=season, y=change_norm, fill=region))

ggplot(dat) +
  geom_boxplot(aes(x=season, y=occurrence, fill=region))

ggplot(dat) +
  geom_boxplot(aes(x=season, y=recurrence, fill=region))

ggplot(dat) +
  geom_boxplot(aes(x=season, y=seasonality, fill=region))

ggplot(dat) +
  geom_boxplot(aes(x=season, y=log(covcrop+0.001), fill=region))

ggplot(dat) +
  geom_boxplot(aes(x=season, y=log(covbuilt+0.001), fill=region))

ggplot(dat) +
  geom_boxplot(aes(x=season, y=log(conv+0.001), fill=region))

ggplot(dat) +
  geom_boxplot(aes(x=season, y=pdsi, fill=region))

#3. Set up bootstraps & season loop----
boot <- 100

set.seed(1234)
season <- unique(dat$season)
out.list <- data.frame()
for(i in 1:length(season)){
  
  season.i <- season[i]
  
  for(j in 1:boot){
    
    #4. Select one point per bird and scale----
    dat.i <- dat %>% 
      dplyr::filter(season==season.i) %>% 
      group_by(id) %>% 
      sample_n(1) %>% 
      ungroup() %>% 
      mutate(crops.s = scale(crops),
             grass.s =scale(grass),
             built.s = scale(built),
             water.s = scale(flooded_vegetation),
             wchange.s = scale(change_norm),
             recur.s = scale(recurrence),
             conv.s = scale(conv),
             drought.s = scale(pdsi)) %>% 
      dplyr::select(region, crops.s, grass.s, built.s, water.s, wchange.s, recur.s, conv.s, drought.s)
    
    #5. NPMANOVA----
    vars.i <- dat.i %>% 
      dplyr::select(crops.s, grass.s, built.s, water.s, wchange.s, recur.s, conv.s, drought.s)
    
    a.i <- adonis(vars.i ~ dat.i$region, method="euclidean", parallel = 4)
    
    #6. Save info----
    out.list <- data.frame(f = a.i[["aov.tab"]]$F.Model[1],
                           r2 = a.i[["aov.tab"]]$R2[1],
                           p = a.i[["aov.tab"]]$`Pr(>F)`[1]) %>% 
      cbind(t(a.i$coefficients[2,])) %>% 
      cbind(t(data.frame(table(dat.i$region))$Freq)) %>% 
      rename(eastn = '1',
             westn = '2') %>% 
      mutate(boot=j,
             season = season.i) %>% 
      rbind(out.list)
  }

}

#7. Summarize----
out.sum <- out.list %>% 
  group_by(season) %>% 
  summarize(p = mean(p),
            r2 = mean(r2),
            crop = mean(crops.s),
            grass = mean(grass.s),
            built = mean(built.s),
            water = mean(water.s),
            wchange = mean(wchange.s),
            recur = mean(recur.s),
            conv = mean(conv.s),
            drought = mean(drought.s)) %>% 
  ungroup()
out.sum

out.covs <- out.list %>% 
  pivot_longer(crops.s:drought.s, names_to = "var", values_to="val")

ggplot(out.covs) +
  geom_hline(aes(yintercept=0), linetype = "dashed") +
  geom_boxplot(aes(x=var, y=val, fill=season))

#8. Test-----
dat.l <- dat %>% 
  mutate(crops.s = scale(crops),
         grass.s =scale(grass),
         built.s = scale(built),
         water.s = scale(flooded_vegetation),
         wchange.s = scale(change_norm),
         recur.s = scale(recurrence),
         conv.s = scale(conv),
         drought.s = scale(pdsi))

#Linear regression
nd <- dat %>% 
  dplyr::select(season, region, id, year) %>% 
  unique()

m1 <- lmer(drought.s ~ season*region + (1|id), data=dat.l, na.action="na.fail", REML=FALSE)
summary(m1)
dredge(m1)
plot(m1)
p1 <- data.frame(p = predict(m1, newdata=nd)) %>% 
  cbind(nd)

ggplot(p1) +
  geom_boxplot(aes(x=season, y=p, colour=region))

#Two-way anova
dat.a <- dat.l %>% 
 group_by(id, season) %>%
 sample_n(1) %>%
 ungroup() %>%
  mutate(crops.s = scale(crops),
         grass.s =scale(grass),
         built.s = scale(built),
         water.s = scale(flooded_vegetation),
         wchange.s = scale(change_norm),
         recur.s = scale(recurrence),
         conv.s = scale(conv),
         drought.s = scale(pdsi))

ggplot(dat.a) +
  geom_boxplot(aes(x=season, y=drought.s, colour=region))

m2 <- aov(crops.s ~ region + season + region:season, data = dat.a)
summary.aov(m2)
plot(m2)

TukeyHSD(m2, which = "region:season")
