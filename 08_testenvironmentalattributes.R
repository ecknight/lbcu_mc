library(tidyverse)
library(vegan)
library(lme4)
library(MuMIn)
library(rstatix)
library(ggridges)

options(scipen=99999)

#1. Read in data-----
dat <- read.csv("Data/LBCU_environvars.csv") %>% 
  mutate(region = case_when(kdecluster==1 ~ "east", 
                            kdecluster==2 ~ "west"),
         conv = covcrop + covbuilt)

#2. Visualize----
ggplot(dat) +
  geom_density_ridges(aes(y=season, x=built, fill=region), alpha=0.3)

ggplot(dat) +
  geom_density_ridges(aes(y=season, x=crops, fill=region), alpha=0.3)

ggplot(dat) +
  geom_density_ridges(aes(y=season, x=grass, fill=region), alpha=0.3)

ggplot(dat) +
  geom_density_ridges(aes(y=season, x=flooded_vegetation, fill=region), alpha=0.3)

ggplot(dat) +
  geom_density_ridges(aes(y=season, x=change_norm, fill=region), alpha=0.3)
  geom_hline(aes(yintercept=0), linetype="dashed")

ggplot(dat) +
  geom_density_ridges(aes(y=season, x=occurrence, fill=region), alpha=0.3)

ggplot(dat) +
  geom_density_ridges(aes(y=season, x=recurrence, fill=region), alpha=0.3)

ggplot(dat) +
  geom_density_ridges(aes(y=season, x=seasonality, fill=region), alpha=0.3)

ggplot(dat) +
  geom_density_ridges(aes(y=season, x=covcrop, fill=region), alpha=0.3)

ggplot(dat) +
  geom_density_ridges(aes(y=season, x=conv, fill=region), alpha=0.3)

ggplot(dat) +
  geom_density_ridges(aes(y=season, x=pdsi, fill=region), alpha=0.3)

#3. Set up bootstraps & season loop----
boot <- 100

set.seed(1234)
season <- unique(dat$season)
vars <- c("crops.s", "grass.s", "built.s", "water.s", "wchange.s", "recur.s", "conv.s", "drought.s")
npmanova.list <- data.frame()
t.list <- data.frame()
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
    
    npmanova.list <- data.frame(f = a.i[["aov.tab"]]$F.Model[1],
                           r2 = a.i[["aov.tab"]]$R2[1],
                           p = a.i[["aov.tab"]]$`Pr(>F)`[1]) %>% 
      cbind(t(a.i$coefficients[2,])) %>% 
      cbind(t(data.frame(table(dat.i$region))$Freq)) %>% 
      rename(eastn = '1',
             westn = '2') %>% 
      mutate(boot=j,
             season = season.i) %>% 
      rbind(npmanova.list)
    
    #6. T-test----
    for(k in 1:length(vars)){
      var.k <- vars[k]
      
      dat.k <- dat.i %>% 
        dplyr::select(region, var.k)
      colnames(dat.k) <- c("region", "var")
      dat.k$var <- dat.k$var[,1]
      
      m.i <- t.test(var ~ region, data=dat.k, var.equal=FALSE)
      
      t.list <- data.frame(t = m.i$statistic,
                               df = m.i$parameter,
                               p = m.i$p.value) %>% 
        mutate(boot=j,
               season = season.i,
               var = var.k) %>% 
        rbind(t.list)
      
    }
    
  }

}

write.csv(npmanova.list, "Data/NPMANOVA.csv", row.names = FALSE)

#7. Summarize NPMANOVA----
npmanova.sum <- npmanova.list %>% 
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
npmanova.sum

npmanova.covs <- npmanova.list %>% 
  pivot_longer(crops.s:drought.s, names_to = "var", values_to="val")

ggplot(npmanova.covs) +
  geom_density_ridges(aes(x=val, y=var, fill=season), alpha=0.3) +
  geom_vline(aes(xintercept=0), linetype = "dashed")


  geom_hline(aes(yintercept=0), linetype = "dashed") +
  geom_boxplot(aes(x=season, y=val, fill=var))

#8. Summarize t-test-----
t.sum <- t.list %>% 
  group_by(var, season) %>% 
  summarize(t = mean(t),
            df = mean(df),
            p = mean(p)) %>% 
  ungroup() %>% 
  dplyr::filter(p < 0.05)
t.sum

ggplot(t.list %>% 
         right_join(t.sum)) +
  geom_point(aes(x=var, y=p, colour=season))