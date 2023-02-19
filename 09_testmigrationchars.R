library(tidyverse)
library(vegan)
library(lme4)
library(MuMIn)
library(sf)
library(adehabitatLT)
library(data.table)

options(scipen=99999)

n <- c("3", "manual")

#1. Get dominant cluster id for each bird----
clust <- read.csv("Data/LBCUKDEClusters.csv") %>%
  dplyr::filter(nclust %in% n,
                !is.na(X)) %>%
  group_by(nclust, id) %>%
  summarize(group = round(mean(group))) %>%
  ungroup()

#2. Read in daily estimate data----
raw <- read.csv("Data/LBCU_FilteredData_Segmented.csv") %>% 
              dplyr::filter(!id %in% c(46768277, 33088, 129945787, 46770723, 46769927, 86872)) %>% 
  inner_join(clust) %>% 
  mutate(region = case_when(nclust=="3" & group==1 ~ "central",
                            nclust=="3" & group==2 ~ "west",
                            nclust=="3" & group==3 ~ "east",
                            nclust=="manual" & group==1 ~ "central",
                            nclust=="manual" & group==2 ~ "west",
                            nclust=="manual" & group==3 ~ "east - inland",
                            nclust=="manual" & group==4 ~ "east - coastal")) %>%
  dplyr::filter(!is.na(region)) %>%
  arrange(id, date)

ggplot(raw) + 
  geom_point(aes(x=X, y=Y, colour=region)) +
  facet_grid(nclust~season)

#3. Wrangle migration duration----
dat.dur <- raw %>% 
  dplyr::filter(segment %in% c("depart", "arrive")) %>% 
  mutate(season = case_when(season=="breed" ~ "springmig",
                            season=="winter" ~ "fallmig",
                            !is.na(season) ~ season)) %>% 
  pivot_wider(id_cols=c(nclust, region, group, id, season, year), names_from=segment, values_from=doy) %>% 
  group_by(season) %>% 
  mutate(duration = arrive - depart, 
         arrive2 = arrive - min(arrive, na.rm=TRUE) + 1,
         depart2 = depart - min(depart, na.rm=TRUE) + 1) %>% 
  dplyr::filter(!is.na(duration)) %>% 
  ungroup() %>% 
  dplyr::filter(duration < 70)
#Remove two outliers that are due to gaps in transmission

write.csv(dat.dur, "Data/MigrationTiming.csv", row.names = FALSE)

#4. Wrangle migration distance----

#Calculate distance
traj <- as.ltraj(xy=raw[,c("X", "Y")],
                 id=paste0(raw$nclust, "-", raw$id),
                 date=raw$date,
                 typeII=FALSE,
                 proj4string = CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"))
names(traj) <- unique(raw$id)

dat.traj <- rbindlist(traj, idcol=TRUE) %>% 
  rename(id='.id', X=x, Y=y) %>% 
  separate(id, into=c("nclust", "id")) %>% 
  mutate(id=as.numeric(id)) %>% 
  arrange(nclust, id, date) %>% 
  dplyr::select(dist) %>% 
  cbind(raw) %>% 
  mutate(dist = dist/1000)

#migrations that weren't completed (from data cleaning script) plus last individual who didn't have data for the first leg of migration
dist.remove <- data.frame(legid=c("1418934379-2021-1spring",
                                  "99900-2021-2fall",
                                  "1146533212-2021-1spring",
                                  "77639184-2016-2fall",
                                  "479364105-2019-1spring",
                                  "46768108-2016-2fall",
                                  "188150741-2021-1spring",
                                  "172070515-2019-2fall",
                                  "1418896449-2021-1spring",
                                  "1418878943-2021-1spring",
                                  "172070319-2016-2fall",
                                  "281981414-2018-1spring",
                                  "290347351-2019-2fall",
                                  "877974163-2020-1spring",
                                  "46768277-2014-2fall",
                                  "145698291-2021-1spring",
                                  "77638376-2015-2fall",
                                  "99902-2021-2fall",
                                  "890834800-2020-1spring",
                                  "1953212690-2011-1spring",
                                  "290352179-2020-1spring")) %>% 
  separate(legid, into=c("id", "year", "season")) %>% 
  mutate(id=as.numeric(id),
         year=as.numeric(year),
         season=ifelse(season=="1spring", "springmig", "fallmig"))

dat.dist <- dat.traj %>% 
  dplyr::filter(season %in% c("fallmig", "springmig") |
                  lead(segment)=="depart" |
                  segment=="arrive") %>% 
  mutate(season = case_when(segment=="stationary" & season=="winter" ~ "springmig",
                            segment=="stationary" & season=="breed" ~ "fallmig",
                            segment=="arrive" & season=="winter" ~ "fallmig",
                            segment=="arrive" & season=="breed" ~ "springmig",
                            !is.na(season) ~ season)) %>% 
  anti_join(dist.remove) %>% 
  group_by(nclust, region, group, id, year, season) %>% 
  summarize(dist = sum(dist)) %>% 
  ungroup()

hist(dat.dist$dist)

write.csv(dat.dist, "Data/MigrationDistance.csv", row.names = FALSE)
  
#5. Wrangle migration rate----
dat.rate <- dat.dist %>% 
  inner_join(dat.dur) %>% 
  mutate(rate = dist/duration)

write.csv(dat.rate, "Data/MigrationRate.csv", row.names = FALSE)

#6. Wrangle number of stopovers/home ranges----
dat.mig <- raw %>% 
  dplyr::filter(season %in% c("springmig", "fallmig")) %>% 
  group_by(nclust, region, group, id, year, season) %>% 
  summarize(n = max(stopovercluster, na.rm=TRUE)) %>% 
  ungroup() %>% 
  mutate(n = ifelse(is.infinite(n), 0, n))

write.csv(dat.mig, "Data/MigrationStopovers.csv", row.names = FALSE)

dat.wint <- raw %>% 
  dplyr::filter(season=="winter", winter!="wintermig") %>% 
  mutate(year = as.numeric(ifelse(doy < 130, year-1, year)),
         cluster = as.numeric(str_sub(winter, -1, -1))) %>% 
  group_by(nclust, season, region, group, id, year) %>% 
  summarize(n = max(cluster)) %>% 
  ungroup()

#7. Wrangle stopover duration---
dat.stop <- raw %>% 
  dplyr::filter(season %in% c("springmig", "fallmig"),
                stopover==1) %>% 
  group_by(nclust, region, id, year, season, stopovercluster) %>% 
  summarize(n = n()) %>% 
  ungroup()

write.csv(dat.stop, "Data/MigrationStopoverLength.csv", row.names = FALSE)

#7. Wrangle home range size----
dat.hr <- read_sf("gis/shp/kde_individual.shp") %>% 
  st_make_valid() %>% 
  st_centroid() %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  dplyr::select(X, Y) %>% 
  cbind(read_sf("gis/shp/kde_individual.shp")) %>% 
  mutate(group = as.numeric(group),
         bird = as.numeric(bird),
         year = as.numeric(year)) %>% 
  left_join(clust %>% 
              rename(bird=id)) %>% 
  data.frame() %>% 
  dplyr::select(-geometry) %>% 
  mutate(region = case_when(nclust=="3" & group==1 ~ "central",
                            nclust=="3" & group==2 ~ "west",
                            nclust=="3" & group==3 ~ "east",
                            nclust=="manual" & group==1 ~ "central",
                            nclust=="manual" & group==2 ~ "west",
                            nclust=="manual" & group==3 ~ "east - inland",
                            nclust=="manual" & group==4 ~ "east - coastal")) %>%
  dplyr::filter(area < 20000) %>% 
  mutate(area.s = scale(area))

write.csv(dat.hr, "Data/WinterHRSize.csv", row.names = FALSE)

#8. Wrangle number of wintering home ranges----
dat.wint <- dat.hr %>% 
  dplyr::filter(season=="winter") %>% 
  group_by(nclust, season, region, group, bird, year) %>% 
  summarize(n=n(),
            X = mean(X),
            Y = mean(Y)) %>% 
  ungroup()

write.csv(dat.wint, "Data/WinterHRs.csv", row.names=FALSE)

#9. Put all the data together for plotting----
dat.all <- dat.mig %>% 
  dplyr::select(-group) %>% 
  mutate(var = "Stopovers") %>% 
  rename(val=n) %>% 
  rbind(dat.hr %>% 
          dplyr::select(nclust, season, bird, year, area, region) %>% 
          rename(id=bird,
                 val=area) %>% 
          mutate(var = "HRarea")) %>% 
  rbind(dat.wint %>% 
          dplyr::select(nclust, bird, year, n, region) %>% 
          rename(id=bird,
                 val=n) %>% 
          mutate(var = "WinterHRs",
                 season="winter")) %>% 
  rbind(dat.rate %>% 
          dplyr::select(-group) %>% 
          pivot_longer(dist:rate, values_to = "val", names_to="var"))

write.csv(dat.all, "Data/MovementBehaviours.csv", row.names = FALSE)
dat.all <- read.csv("Data/MovementBehaviours.csv")

#9. Visualize----
#Departure
ggplot(dat.dur) +
  geom_boxplot(aes(x=region, y=depart)) +
  facet_grid(nclust~season, scales="free")

ggplot(dat.dur) +
  ggridges::geom_density_ridges(aes(y=season, x=depart, fill=region), alpha=0.3)

#Arrival
ggplot(dat.dur) +
  geom_boxplot(aes(x=region, y=arrive)) +
  facet_grid(nclust~season, scales="free")

ggplot(dat.dur) +
  ggridges::geom_density_ridges(aes(y=season, x=arrive, fill=region), alpha=0.3)

#Duration
ggplot(dat.dur) +
  geom_boxplot(aes(x=region, y=duration)) +
  facet_grid(nclust~season)

ggplot(dat.dur) +
  ggridges::geom_density_ridges(aes(y=season, x=duration, fill=region), alpha=0.3)

#Distance
ggplot(dat.dist) +
  geom_boxplot(aes(x=region, y=log(dist))) +
  facet_grid(nclust~season)

ggplot(dat.dist) +
  ggridges::geom_density_ridges(aes(y=season, x=dist, fill=region), alpha=0.3)

#Rate
ggplot(dat.rate) +
  geom_boxplot(aes(x=region, y=log(rate))) +
  facet_grid(nclust~season)

ggplot(dat.rate) +
  ggridges::geom_density_ridges(aes(y=season, x=rate, fill=region), alpha=0.3)

#Number of stopovers
ggplot(dat.mig) +
  geom_histogram(aes(x=n, fill=region)) +
  facet_grid(nclust~season)

ggplot(dat.mig) +
  ggridges::geom_density_ridges(aes(y=season, x=n, fill=region), alpha=0.3)

#Stopover duration
ggplot(dat.stop) +
  geom_boxplot(aes(x=region, y=n)) +
  facet_grid(nclust~season)

ggplot(dat.stop) +
  ggridges::geom_density_ridges(aes(y=season, x=n, fill=region), alpha=0.3)

#Area
ggplot(dat.hr) +
  geom_boxplot(aes(x=region, y=log(area))) +
  facet_grid(nclust~season,scales="free")

ggplot(dat.hr) +
  ggridges::geom_density_ridges(aes(y=season, x=log(area), fill=region), alpha=0.3)

#Wintering home ranges
ggplot(dat.wint) +
  geom_boxplot(aes(x=region, y=n)) +
  facet_grid(nclust~season,scales="free")

ggplot(dat.wint) +
  ggridges::geom_density_ridges(aes(y=season, x=n, fill=region), alpha=0.3)

ggplot(dat.wint) +
  geom_histogram(aes(x=n)) +
  facet_wrap(~region)

#10. Test----
n <- c("3", "manual")

aic.out <- data.frame()
for(i in 1:length(n)){
  
  aic.list <- list()
  
  dat.i <- dat.all %>% 
    dplyr::filter(nclust==n[i])
  
  #Departure
  dat.j <- dplyr::filter(dat.i, var=="depart") %>% 
    dplyr::filter(!is.na(val))
  
  m1 <- lmer(val ~ region*season + (1|id), data=dat.j, na.action="na.fail", REML=FALSE)
  aic.list[[1]] <- dredge(m1)
  aic.list[[1]]$mod = "departure"
  aic.list[[1]]
  # summary(m1)
  # plot(m1)
  # pdep <- data.frame(pred = predict(m1)) %>% 
  #   cbind(dat.dur)
  # 
  # ggplot(dat.dur) +
  #   geom_boxplot(aes(x=region, y=depart)) +
  #   geom_boxplot(aes(x=region, y=pred), data=pdep, fill=NA, colour="red") +
  #   facet_wrap(~season, scales="free")
  
  #Arrival
  dat.j <- dplyr::filter(dat.i, var=="arrive") %>% 
    dplyr::filter(!is.na(val))
  
  m2 <- lmer(val ~ region*season + (1|id), data=dat.j, na.action="na.fail", REML=FALSE)
  aic.list[[2]] <- dredge(m2)
  aic.list[[2]]$mod = "arrival"
  aic.list[[2]]
  # summary(m2)
  # plot(m2)
  # 
  # parr <- data.frame(pred = predict(m2)) %>% 
  #   cbind(dat.dur)
  # 
  # ggplot(dat.dur) +
  #   geom_boxplot(aes(x=region, y=arrive)) +
  #   geom_boxplot(aes(x=region, y=pred), data=parr, fill=NA, colour="red") +
  #   facet_wrap(~season, scales="free")
  
  #Duration
  dat.j <- dplyr::filter(dat.i, var=="duration") %>% 
    dplyr::filter(!is.na(val))
  
  m3 <-lmer(val ~ region*season + (1|id), data=dat.j, na.action="na.fail", REML=FALSE)
  aic.list[[3]] <- dredge(m3)
  aic.list[[3]]$mod = "duration"
  aic.list[[3]]
  # summary(m3)
  # plot(m3)
  # 
  # pdur <- data.frame(pred = predict(m3)) %>% 
  #   cbind(dat.dur)
  # 
  # ggplot(dat.dur) +
  #   geom_boxplot(aes(x=region, y=duration)) +
  #   geom_boxplot(aes(x=region, y=pred), data=pdur, fill=NA, colour="red") +
  #   facet_wrap(~season, scales="free")
  
  #Distance
  dat.j <- dplyr::filter(dat.i, var=="dist") %>% 
    dplyr::filter(!is.na(val))
  
  m4 <- lmer(val ~ region*season + (1|id), data=dat.j, na.action="na.fail", REML=FALSE)
  aic.list[[4]] <- dredge(m4)
  aic.list[[4]]$mod = "distance"
  aic.list[[4]]
  # summary(m4)
  # plot(m4)
  # 
  # pdist <- data.frame(pred = predict(m4)) %>% 
  #   cbind(dat.dist)
  # 
  # ggplot(dat.dist) +
  #   geom_boxplot(aes(x=region, y=dist)) +
  #   geom_boxplot(aes(x=region, y=pred), data=pdist, fill=NA, colour="red") +
  #   facet_wrap(~season, scales="free")
  # 
  # ggplot(dat.dist) +
  #   geom_boxplot(aes(x=season, y=dist)) +
  #   geom_boxplot(aes(x=season, y=pred), data=pdist, fill=NA, colour="red") +
  #   facet_wrap(~region, scales="free")
  
  #Rate
  dat.j <- dplyr::filter(dat.i, var=="rate") %>% 
    dplyr::filter(!is.na(val))
  
  m5 <- lmer(val ~ region*season + (1|id), data=dat.j, na.action="na.fail", REML=FALSE)
  aic.list[[5]] <- dredge(m5)
  aic.list[[5]]$mod = "rate"
  aic.list[[5]]
  # summary(m5)
  # plot(m5)
  # 
  # prate <- data.frame(pred = predict(m5)) %>% 
  #   cbind(dat.rate)
  # 
  # ggplot(dat.rate) +
  #   geom_boxplot(aes(x=region, y=rate)) +
  #   geom_boxplot(aes(x=region, y=pred), data=prate, fill=NA, colour="red") +
  #   facet_wrap(~season, scales="free")
  
  #Stopovers
  dat.j <- dplyr::filter(dat.i, var=="Stopovers") %>% 
    dplyr::filter(!is.na(val))
  
  m6 <- glmer(val ~ region*season + (1|id), dat.j, na.action="na.fail", family="poisson")
  aic.list[[6]] <- dredge(m6)
  aic.list[[6]]$mod = "stopover"
  aic.list[[6]]
  # m6 <- glmer(n ~ region + season + (1|id), data=dat.mig, na.action="na.fail", family="poisson")
  # summary(m6)
  # plot(m6)
  # 
  # pstop <- data.frame(pred = predict(m6, type="response")) %>% 
  #   cbind(dat.mig)
  # 
  # ggplot(dat.mig) +
  #   geom_boxplot(aes(x=region, y=n)) +
  #   geom_boxplot(aes(x=region, y=pred), data=pstop, fill=NA, colour="red") +
  #   facet_wrap(~season, scales="free")
  
  #Home range area
  dat.j <- dplyr::filter(dat.i, var=="HRarea") %>% 
    dplyr::filter(!is.na(val))
  
  m7 <- lmer(log(val) ~ region*season + (1|id), data=dat.j, na.action="na.fail", REML = FALSE)
  aic.list[[7]] <- dredge(m7)
  aic.list[[7]]$mod = "area"
  aic.list[[7]]
  # plot(m7)
  # summary(m7)
  # 
  # phr <- data.frame(pred = predict(m7)) %>% 
  #   cbind(dat.hr)
  # 
  # ggplot(dat.hr) +
  #   geom_boxplot(aes(x=region, y=log(area))) +
  #   geom_boxplot(aes(x=region, y=pred), data=phr, fill=NA, colour="red") +
  #   facet_wrap(~season, scales="free")
  
  #Number of home ranges
  dat.j <- dplyr::filter(dat.i, var=="WinterHRs") %>% 
    dplyr::filter(!is.na(val))
  
  m8 <- glmer(val ~ region + (1|id), data=dat.j, na.action="na.fail", family="poisson")
  aic.list[[8]] <- dredge(m8)
  aic.list[[8]]$mod = "winterhrs"
  aic.list[[8]]
  # plot(m8)
  # summary(m8)
  
  # #Stopover duration
  # m9 <- lmer(n ~ region*season + (1|id), data=dat.stop, na.action="na.fail", REML=FALSE)
  # aic.list[[9]] <- dredge(m9)
  # aic.list[[9]]
  # m9 <- lmer(n ~ season + (1|id), data=dat.stop, na.action="na.fail", REML=FALSE)
  # plot(m9)
  # summary(m9)
  # 
  # pstop <- data.frame(pred = predict(m9, type="response")) %>% 
  #   cbind(dat.stop)
  # 
  # ggplot(dat.stop) +
  #   geom_boxplot(aes(x=region, y=n)) +
  #   geom_boxplot(aes(x=region, y=pred), data=pstop, fill=NA, colour="red") +
  #   facet_wrap(~season, scales="free")
  
  #11. Collapse AIC results----
  aic.out <- rbindlist(aic.list, fill=TRUE) %>% 
    arrange(mod, -weight) %>% 
    mutate(weight = round(weight, 3),
           nclust = n[i]) %>% 
    rbind(aic.out)
  
}

aic.0 <- aic.out %>% 
  dplyr::filter(delta==0)

#12. Play with trend piece----
trend <- read.csv("Data/LBCU_trend_gamye.csv")

trend.sum <- trend %>% 
  dplyr::select(nclust, knncluster, Trend, 'Trend_Q0.025', 'Trend_Q0.975') %>% 
  unique() %>% 
  rename(up = 'Trend_Q0.975', down = 'Trend_Q0.025', group=knncluster) %>% 
  mutate(region = case_when(nclust=="3" & group==1 ~ "central",
                            nclust=="3" & group==2 ~ "west",
                            nclust=="3" & group==3 ~ "east",
                            nclust=="manual" & group==1 ~ "central",
                            nclust=="manual" & group==2 ~ "west",
                            nclust=="manual" & group==3 ~ "east - inland",
                            nclust=="manual" & group==4 ~ "east - coastal"))

trend.all <- dat.all %>% 
  left_join(trend.sum)

ggplot(trend.all) +
  geom_point(aes(x=val, y=Trend, colour=season)) +
  geom_smooth(aes(x=val, y=Trend, colour=season), method="lm") +
  facet_wrap(var~nclust, scales="free")
