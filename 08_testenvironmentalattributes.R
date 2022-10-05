library(tidyverse)
library(vegan)
library(lme4)
library(MuMIn)
library(rstatix)
library(ggridges)
library(data.table)

options(scipen=99999)

#1. Read in data-----
raw <- read.csv("Data/LBCU_environvars.csv") %>% 
  mutate(conv = covcrop + covbuilt,
         region = case_when(kdecluster==1 ~ "central", 
                            kdecluster==2 ~ "west",
                            kdecluster==3 ~ "east")) %>% 
  dplyr::select(-kdecluster, -season, -year, -cluster) %>% 
  separate(id, into=c("kdecluster", "season", "birdid", "year", "cluster"), remove=FALSE)

#2. Visualize----
ggplot(raw) +
  geom_density_ridges(aes(y=season, x=log(built+0.000001), fill=region), alpha=0.3)

ggplot(raw) +
  geom_density_ridges(aes(y=season, x=crops, fill=region), alpha=0.3)

ggplot(raw) +
  geom_density_ridges(aes(y=season, x=grass, fill=region), alpha=0.3)

ggplot(raw) +
  geom_density_ridges(aes(y=season, x=log(flooded_vegetation), fill=region), alpha=0.3)

ggplot(raw) +
  geom_density_ridges(aes(y=season, x=change_norm, fill=region), alpha=0.3) +
  geom_hline(aes(yintercept=0), linetype="dashed")

ggplot(raw) +
  geom_density_ridges(aes(y=season, x=occurrence, fill=region), alpha=0.3)

ggplot(raw) +
  geom_density_ridges(aes(y=season, x=recurrence, fill=region), alpha=0.3)

ggplot(raw) +
  geom_density_ridges(aes(y=season, x=seasonality, fill=region), alpha=0.3)

ggplot(raw) +
  geom_density_ridges(aes(y=season, x=covcrop, fill=region), alpha=0.3)

ggplot(raw) +
  geom_density_ridges(aes(y=season, x=log(conv+0.000001), fill=region), alpha=0.3)

ggplot(raw) +
  geom_density_ridges(aes(y=season, x=pdsi, fill=region), alpha=0.3)

#3. Find K for NMDS----
set.seed(1234)
seasons <- unique(raw$season)
stress <- data.frame()
for(i in 1:length(seasons)){
  
  dat.covs <- raw %>% 
    dplyr::filter(season==seasons[i]) %>% 
    group_by(birdid) %>% 
    sample_n(1) %>% 
    ungroup() %>% 
    dplyr::select(built, change_norm, crops, flooded_vegetation, seasonality, grass, pdsi, recurrence, conv) %>% 
    mutate_all(~((.-min(.)))) %>% 
    mutate_all(~scale(., center=FALSE))
  
  nm2 <- metaMDS(dat.covs, k=2, trymax=100, maxit=1000, trace=0)
  nm3 <- metaMDS(dat.covs, k=3, trymax=100, maxit=1000, trace=0)
  nm4 <- metaMDS(dat.covs, k=4, trymax=100, maxit=1000, trace=0)
  nm5 <- metaMDS(dat.covs, k=5, trymax=100, maxit=1000, trace=0)
  nm6 <- metaMDS(dat.covs, k=6, trymax=100, maxit=1000, trace=0)
  nm7 <- metaMDS(dat.covs, k=6, trymax=100, maxit=1000, trace=0)
  nm8 <- metaMDS(dat.covs, k=6, trymax=100, maxit=1000, trace=0)
  nm9 <- metaMDS(dat.covs, k=6, trymax=100, maxit=1000, trace=0)
  nm10 <- metaMDS(dat.covs, k=6, trymax=100, maxit=1000, trace=0)
  
  stress <- data.frame(stress = c(nm2$stress, nm3$stress, nm4$stress, nm5$stress, nm6$stress, nm7$stress, nm8$stress, nm9$stress, nm10$stress),
                      k =c(2:10),
                      season = seasons[i]) %>% 
    rbind(stress)
  
  print(paste0("Finished season ", seasons[i]))
  
}

ggplot(stress)+
  geom_point(aes(x=k, y=stress)) +
  geom_line(aes(x=k, y=stress)) +
  facet_wrap(~season)

#UGH six is optimal, but 2 is actually not bad all things considered

#4. Set up bootstrap----
boot <- 100

set.seed(123)
season <- unique(raw$season)
stress.list <- list()
scores.list <- list()
covscores.list <- list()
mrpp.list <- list()
for(i in 1:length(season)){
  
  season.i <- season[i]
  
  stress.i <- data.frame()
  scores.i <- data.frame()
  covscores.i <- data.frame()
  mrpp.i <- data.frame()
  j <- 1
  while(j <= boot){
    
    #5. Select one point per bird and scale----
    dat.i <- raw %>% 
      dplyr::filter(season==season.i) %>% 
      group_by(birdid) %>% 
      sample_n(1) %>% 
      ungroup()
    
    covs.i <- dat.i %>% 
      dplyr::select(built, change_norm, crops, flooded_vegetation, grass, pdsi, recurrence, conv, seasonality) %>% 
      mutate_all(~((.-min(.)))) %>% 
      mutate_all(~scale(., center=FALSE))
    
    #6. Run NMDS----
    nmds.i <- metaMDS(covs.i, k=2, trace=0, trymax=100, maxit=1000)
    
    if(nmds.i$converged==TRUE){
      
      #save out diagnostics
      stress.i <- data.frame(stress = nmds.i$stress,
                             converge = nmds.i$converged) %>%
        mutate(boot = j) %>%
        rbind(stress.i)
      
      #7. Save out coords for each individual----
      scores.i <- scores(nmds.i, "sites") %>%
        data.frame() %>%
        mutate(boot = j,
               region = dat.i$region, 
               id = dat.i$birdid) %>%
        rbind(scores.i)
      
      #8. Save out coords for each covariate----
      covscores.i <- scores(nmds.i, "species") %>%
        data.frame() %>%
        mutate(cov=rownames(scores(nmds.i, "species"))) %>%
        mutate(boot = j) %>%
        rbind(covscores.i)
      
      #9. MRPP for group diffs----
      mrpp.mod.i <- mrpp(covs.i, dat.i$region, distance = "mahalanobis")
      
      mrpp.i <- data.frame(delta = mrpp.mod.i[["classdelta"]],
                           p = mrpp.mod.i[["Pvalue"]],
                           region = names(mrpp.mod.i[["classdelta"]])) %>% 
        mutate(boot = j) %>%
        rbind(mrpp.i)
      
      j <- j + 1
      
    }

    print(paste0("Finished bootstrap ", j, " of ", boot))
    
  }
  
  stress.list[[i]] <- stress.i %>% 
    mutate(season=season.i)
  scores.list[[i]] <- scores.i %>% 
    mutate(season=season.i)
  covscores.list[[i]] <- covscores.i %>% 
    mutate(season=season.i)
  mrpp.list[[i]] <- mrpp.i %>% 
    mutate(season=season.i)
  
  print(paste0("Finished season ", i, " of ", length(season)))
  
}

stress.out <- rbindlist(stress.list)
scores.out <- rbindlist(scores.list)
covscores.out <- rbindlist(covscores.list)
mrpp.out <- rbindlist(mrpp.list)

write.csv(stress.out, "Data/NMDSStress.csv", row.names = FALSE)
write.csv(scores.out, "Data/NMDSScores.csv", row.names = FALSE)
write.csv(covscores.out, "Data/NMDSCovscores.csv", row.names = FALSE)
write.csv(mrpp.out, "Data/MRPP.csv", row.names = FALSE)

#10. Summarize mrpp----
mrpp.p <- mrpp.out %>% 
  dplyr::select(season, p) %>% 
  # dplyr::select(season, p, boot) %>% 
  # unique() %>% 
  group_by(season) %>%
  summarize(p.mn = mean(p),
            p.sd = sd(p)) %>% 
  ungroup()
mrpp.p

mrpp.delta <- mrpp.out %>% 
  group_by(season, region) %>% 
  summarize(delta.mn = mean(delta),
            delta.sd = sd(delta),
            delta.up = quantile(delta, 0.975),
            delta.lw = quantile(delta, 0.025)) %>% 
  ungroup()
mrpp.delta

mrpp.delta$region <- factor(mrpp.delta$region, levels=c("west", "central", "east"))

ggplot(mrpp.delta) +
  geom_point(aes(x=region, y=delta.mn)) +
  geom_errorbar(aes(x=region, ymin = delta.lw, ymax = delta.up)) +
  facet_wrap(~season)

#11. NMDS stress----
stress.out <- read.csv("Data/NMDSStress.csv") %>% 
  dplyr::filter(converge==TRUE)

table(stress.out$season, stress.out$converge)
#ACKKKKK

stress.sum <- stress.out  %>% 
  group_by(season) %>% 
  summarize(stress.mn = mean(stress),
            stress.sd = sd(stress)) %>% 
  ungroup
stress.sum

#11. Summarize nmds----
scores.sum <- scores.out %>% 
  group_by(season, boot) %>% 
  mutate(id = row_number()) %>% 
  ungroup() %>% 
  group_by(season, region, id) %>% 
  summarize(x.mn = mean(NMDS1),
            x.up = quantile(NMDS1, 0.975),
            x.lw = quantile(NMDS1, 0.025),
            y.mn = mean(NMDS2),
            y.up = quantile(NMDS2, 0.975),
            y.lw = quantile(NMDS2, 0.025)) %>% 
  ungroup() %>% 
  dplyr::filter(!(season=="winter" & id==1))

covscores.sum <- covscores.out %>% 
  group_by(season, cov) %>% 
  summarize(x.mn = mean(NMDS1),
            x.up = quantile(NMDS1, 0.975),
            x.lw = quantile(NMDS1, 0.025),
            y.mn = mean(NMDS2),
            y.up = quantile(NMDS2, 0.975),
            y.lw = quantile(NMDS2, 0.025)) %>% 
  ungroup()

ggplot() +
  geom_point(data=scores.sum, aes(x=x.mn, y=y.mn, fill=region), pch=21, size=2) +
  geom_segment(data=covscores.sum, aes(x=0, y=0, xend=x.mn, yend=y.mn, colour=cov), arrow = arrow(length = unit(0.2, "cm"))) +
  facet_wrap(~season, scales="free")
