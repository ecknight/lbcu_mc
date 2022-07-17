library(tidyverse)
library(ggmap)
library(sf)
library(gridExtra)
library(mapdata)
library(maptools)
library(ggpubr)
library(ggridges)

my.theme <- theme_classic() +
  theme(text=element_text(size=12, family="Arial"),
        axis.text.x=element_text(size=12),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        plot.title=element_text(size=12, hjust = 0.5))

map.theme <- theme_nothing() +
  theme(text=element_text(size=12, family="Arial"),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1))

nam <- map_data("usa") %>% 
  rbind(map_data("worldHires", "Canada")) %>% 
  rbind(map_data("worldHires", "Mexico")) %>% 
  mutate(group = paste0(region,"-", group))

#1. Study area figure----
#Get data and wrangle for deployment locations----
dep <- read.csv("/Users/ellyknight/Documents/SMBC/Analysis/lbcu_exploration/Data/LBCUCleanedData.csv") %>% 
  group_by(id) %>% 
  dplyr::filter(row_number()==1) %>% 
  ungroup()

plot.sa <- ggplot() +
  geom_polygon(data=nam, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  coord_sf(xlim=c(-170, -30), expand = FALSE, crs=4326) +
  geom_point(data=dep, aes(x=long, y=lat, fill=study), pch=21, colour="grey85", alpha = 0.5, size=5) +
  xlab("") +
  ylab("") +
  map.theme +
  theme(legend.position = "right")

ggsave(plot.sa, filename="Figs/Fig1StudyArea.jpeg", width=8, height=6)

#2. Clustering figure----

#2a. Wrangle----
clust.raw <- read.csv("Data/LBCUKDEClusters.csv")

clust <- clust.raw %>% 
  group_by(id, season, nclust) %>% 
  summarize(X = mean(X, na.rm=FALSE), 
            Y = mean(Y, na.rm=FALSE),
            kdecluster = round(mean(kdecluster))) %>% 
  ungroup()

clust.ll <- clust  %>% 
  dplyr::filter(!is.na(X)) %>% 
  st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
  st_transform(crs=4326) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  rename(lat = Y, lon = X) %>% 
  cbind(clust %>% 
          dplyr::filter(!is.na(X)))

clust.ll$season <- factor(clust.ll$season, levels=c("breed", "fallmig", "winter", "springmig"),
                             labels=c("Breeding", "Postbreeding\nmigration", "Nonbreeding", "Prebreeding\nmigration"))
clust.ll$nclust <- factor(clust.ll$nclust, levels=c(2,3,4,5,6), labels=c("2 regions", "3 regions", "4 regions", "5 regions", "6 regions"))
  
#2b. Plot----
plot.clust <- ggplot(clust.ll) +
  geom_polygon(data = nam, aes(x=long, y = lat, group = group), fill = "gray75", color="gray90", size=0.3) +
  geom_point(aes(x=lon, y=lat, fill=factor(kdecluster)), size = 3, pch=21, colour="grey90") +
  coord_sf(xlim=c(min(clust.ll$lon)-5, max(clust.ll$lon)+5), ylim = c(min(clust.ll$lat)-5, max(clust.ll$lat)+5), expand = FALSE, crs=4326) +
  facet_grid(season ~ nclust) +
  xlab("") +
  ylab("") +
  map.theme +
  theme(legend.position = "bottom") +
  scale_fill_viridis_d(name="Region") +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  theme(panel.spacing = unit(0.2, "lines"))

ggsave(plot.clust, filename="Figs/Fig2Cluster.jpeg", width = 8, height = 10)

#3. Connectivity figure----
mn <- read.csv("Data/LBCUMCLocations.csv")

#3a. Maps---
plot.bw <- ggplot(mn %>% dplyr::filter(season%in% c("breed", "winter"))) +
  geom_polygon(data = nam, aes(x=long, y = lat, group = group), fill = "gray75", color="gray90", size=0.3) +
  geom_path(aes(x=lon, y=lat, group=id), colour="grey20") +
  geom_point(aes(x=lon, y=lat, fill=season), pch=21, colour="grey20", size=2.5, alpha = 0.7) +
  coord_sf(xlim=c(min(mn$lon)-5, max(mn$lon)+5), ylim = c(min(mn$lat)-5, max(mn$lat)+5), expand = FALSE, crs=4326) +
  scale_fill_manual(values=c("cyan4", "darkgoldenrod1")) +
  xlab("") +
  ylab("") +
  map.theme +
  ggtitle("Breeding to nonbreeding")
#plot.bw

plot.bf <- ggplot(mn %>% dplyr::filter(season %in% c("breed", "fallmig"))) +
  geom_polygon(data = nam, aes(x=long, y = lat, group = group), fill = "gray75", color="gray90", size=0.3) +
  geom_path(aes(x=lon, y=lat, group=id), colour="grey20") +
  geom_point(aes(x=lon, y=lat, fill=season), pch=21, colour="grey20", size=2.5, alpha = 0.7) +
  coord_sf(xlim=c(min(mn$lon)-5, max(mn$lon)+5), ylim = c(min(mn$lat)-5, max(mn$lat)+5), expand = FALSE, crs=4326) +
  scale_fill_manual(values=c("cyan4", "darkred")) +
  xlab("") +
  ylab("") +
  map.theme +
  ggtitle("Breeding to postbreeding migration")
#plot.bf

plot.bs <- ggplot(mn %>% dplyr::filter(season %in% c("breed", "springmig"))) +
  geom_polygon(data = nam, aes(x=long, y = lat, group = group), fill = "gray75", color="gray90", size=0.3) +
  geom_path(aes(x=lon, y=lat, group=id), colour="grey20") +
  geom_point(aes(x=lon, y=lat, fill=season), pch=21, colour="grey20", size=2.5, alpha = 0.7) +
  coord_sf(xlim=c(min(mn$lon)-5, max(mn$lon)+5), ylim = c(min(mn$lat)-5, max(mn$lat)+5), expand = FALSE, crs=4326) +
  scale_fill_manual(values=c("cyan4", "darkolivegreen4")) +
  xlab("") +
  ylab("") +
  map.theme +
  ggtitle("Breeding to prebreeding migration")
#plot.bs

#3b. Connectivity----
results <- read.csv("Data/LBCUMigConnectivity.csv")

sum <- results %>% 
  dplyr::select(MC, MClow, MChigh, nclust, season, boot) %>% 
  unique() %>% 
  group_by(season, nclust) %>% 
  summarize(MC = mean(MC),
            MClow = mean(MClow),
            MChigh = mean(MChigh)) %>% 
  ungroup() %>% 
  rbind(results %>% 
          dplyr::select(MC, MClow, MChigh, nclust, season, boot) %>% 
          unique() %>% 
          group_by(nclust) %>% 
          summarize(MC = mean(MC),
                    MClow = mean(MClow),
                    MChigh = mean(MChigh)) %>% 
          ungroup() %>% 
          mutate(season="Combined"))

mc.bf <- ggplot(sum %>% dplyr::filter(season=="fallmig")) +
  geom_point(aes(x=nclust, y=MC)) +
  geom_errorbar(aes(x=nclust, ymin=MClow, ymax=MChigh)) +
  my.theme +
  ylab("Strength of migratory connectivity") +
  xlab("") +
  ylim(c(0,1))

mc.bw <- ggplot(sum %>% dplyr::filter(season=="winter")) +
  geom_point(aes(x=nclust, y=MC)) +
  geom_errorbar(aes(x=nclust, ymin=MClow, ymax=MChigh)) +
  my.theme +
  ylab("") +
  xlab("Number of regions") +
  ylim(c(0,1))

mc.bs <- ggplot(sum %>% dplyr::filter(season=="springmig")) +
  geom_point(aes(x=nclust, y=MC)) +
  geom_errorbar(aes(x=nclust, ymin=MClow, ymax=MChigh)) +
  my.theme +
  ylab("") +
  xlab("") +
  ylim(c(0,1))

#3c. Put it together----
plot.mc <- grid.arrange(plot.bf, plot.bw, plot.bs, mc.bf, mc.bw, mc.bs,
             widths=c(2,2,2),
             heights=c(4,2),
             layout_matrix = rbind(c(1,2,3),
                                   c(4,5,6)))

ggsave(plot.mc, filename="Figs/Fig3MC.jpeg", width=10, height=8)

#4. Trend----

#4a. Read in data----
trend.list <- read.csv("Data/LBCU_trend_gamye.csv") %>% 
  st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
  st_transform(crs=4326) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  rename(lat = Y, lon = X) %>% 
  cbind(read.csv("Data/LBCU_trend_gamye.csv")) %>% 
  mutate(Region = ifelse(Region==1, "east", "west"))

trend <- trend.list %>% 
  dplyr::select(Region, Trend, 'Trend_Q0.025', 'Trend_Q0.975') %>% 
  unique() %>% 
  rename(up = 'Trend_Q0.975', down = 'Trend_Q0.025')

#4b. BBS Regions----
clust <- read.csv("Data/LBCUKDEClusters.csv")  %>% 
  dplyr::filter(season=="breed") %>% 
  st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
  st_transform(crs=4326) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  rename(lat = Y, lon = X) %>% 
  cbind(read.csv("Data/LBCUKDEClusters.csv") %>% 
          dplyr::filter(season=="breed")) %>% 
  dplyr::filter(season=="breed",
                nclust==2) %>% 
  dplyr::select(id, kdecluster, lat, lon) %>% 
  unique() %>% 
  left_join(trend %>% 
              dplyr::select(svmcluster, Trend) %>% 
              rename(kdecluster = svmcluster) %>% 
              unique()) %>% 
  mutate(Region = ifelse(kdecluster==1, "east", "west"))

plot.map <- ggplot(trend.list) +
  geom_polygon(data = nam, aes(x=long, y = lat, group = group), fill = "gray75", color="gray90", size=0.3) +
  geom_point(aes(x=lon, y=lat, colour=Trend)) +
  geom_point(data=clust, aes(x=lon, y=lat, colour=Trend), pch=21, fill="white", size=3) +
  scale_colour_gradient2(high="blue", low="red") +
  coord_sf(xlim=c(min(trend.list$lon)-5, max(trend.list$lon)+5), ylim = c(min(trend.list$lat)-5, max(trend.list$lat)+5), expand = FALSE, crs=3857) +
  map.theme +
  xlab("") +
  ylab("")
plot.map

#4c. Trend----
plot.bar <- ggplot(trend) +
  geom_hline(aes(yintercept=0), linetype="dashed", colour="grey30") +
  geom_point(aes(x=Region, y=Trend, colour=Trend), size=2) +
  geom_errorbar(aes(x=Region, ymin=down, ymax=up, colour=Trend)) +
  scale_colour_gradient2(high="blue", low="red") +
  ylab("Population trend") +
  my.theme +
  theme(legend.position="none")
plot.bar

#4d. Trajectories----
indices <- read.csv("Data/LBCU_indices_gamye.csv") %>% 
  mutate(Region = ifelse(Region==1, "east", "west")) %>% 
  left_join(trend %>% 
              dplyr::select(Region, Trend) %>% 
              unique())

plot.index <- ggplot(indices) +
  geom_ribbon(aes(x=Year, ymin=Index_q_0.025, ymax=Index_q_0.975, group=factor(Region)), alpha = 0.5) +
  geom_line(aes(x=Year, y=Index, group=Region, colour=Trend)) +
  scale_colour_gradient2(high="blue", low="red", name="Population\ntrend") +
  ylab("Relative abundance index") +
  my.theme +
  theme(legend.position="none")
plot.index

#4e. Legend----
plot.legend <- ggplot(indices) +
  geom_ribbon(aes(x=Year, ymin=Index_q_0.025, ymax=Index_q_0.975, group=factor(Region)), alpha = 0.5) +
  geom_line(aes(x=Year, y=Index, group=Region, colour=Trend)) +
  scale_colour_gradient2(high="blue", low="red", name="Population\ntrend") +
  ylab("Relative abundance index")

legend <- get_legend(plot.legend)

#4f. Put together----
ggsave(grid.arrange(plot.map, plot.bar, plot.index, legend,
                    widths = c(3,3,3,1), heights = c(3),
                    layout_matrix = rbind(c(1,2,3,4))), height=4, width=12, units='in', filename="Figs/Fig3Trend.jpeg")

#5. Behavioural attributes----

#5a. Get data----
dat.mig <- read.csv("Data/MigrationStopovers.csv") %>% 
  rename(Region = region)
dat.mig$season <- factor(dat.mig$season, levels=c("fallmig", "springmig"),
                         labels=c("Postbreeding\nmigration", "Prebreeding\nmigration"))
dat.mig$Region <- factor(dat.mig$Region, levels=c("west", "east"), labels=c("West", "East"))

dat.dur <- read.csv("Data/MigrationTiming.csv") %>% 
  rename(Region = region)
dat.dur$season <- factor(dat.dur$season, levels=c("fallmig", "springmig"),
                         labels=c("Postbreeding\nmigration", "Prebreeding\nmigration"))
dat.dur$Region <- factor(dat.dur$Region, levels=c("west", "east"), labels=c("West", "East"))

#5b. Departure----
plot.dep <- ggplot(dat.dur) +
  geom_density_ridges(aes(x=depart, y=season, fill=Region), alpha = 0.2) + 
  my.theme +
  xlab("Day of migration departure season") +
  ylab("") +
  theme(legend.position = "none")

#5c. Arrival----
plot.arr <- ggplot(dat.dur) +
  geom_density_ridges(aes(x=arrive, y=season, fill=Region), alpha = 0.2) + 
  my.theme +
  xlab("Day of migration arrival season") +
  ylab("") +
  theme(axis.text.y = element_blank())
#plot.arr

#5d. Duration----
plot.dur <- ggplot(dat.dur) +
  geom_density_ridges(aes(x=duration, y=season, fill=Region), alpha = 0.2) + 
  my.theme +
  xlab("Days of migration") +
  ylab("") +
  theme(legend.position = "none")
#plot.arr

#5e. Stopovers----
plot.stop <- ggplot(dat.mig) +
  geom_density_ridges(aes(x=n, y=season, fill=Region), alpha = 0.2, stat="binline") + 
  my.theme +
  xlab("Number of migration stopovers") +
  ylab("") +
  theme(axis.text.y = element_blank())
#plot.stop

#5f. Put together----
plot.behav <- grid.arrange(plot.dep, plot.arr, plot.dur, plot.stop,
             widths = c(4,4), heights = c(4,4),
             layout_matrix = rbind(c(1,2), c(3,4)))

ggsave(plot.behav, filename="Figs/Fig5Behave.jpeg", width = 10, height = 6)

#6. Environmental attributes----

#6a. Wrangle----
np <- read.csv("Data/NPMANOVA.csv") %>% 
  pivot_longer(crops.s:drought.s, names_to = "var", values_to="val") 
np$var <- factor(np$var, levels=c("water.s", "recur.s", "wchange.s", "grass.s", "conv.s", "crops.s", "built.s", "drought.s"),
                 labels=c("Wetland", "Water\nrecurrence", "Water\nchange", "Grassland", "Grassland\nconversion", "Cropland", "Urban", "Drought"))
np$season <- factor(np$season, levels=c("breed", "fallmig", "winter", "springmig"),
                    labels = c("Breed", "Postbreeding\nmigration", "Nonbreeding", "Prebreeding\nmigration"))

#6b. Plot----
plot.np <- ggplot(np) +
  geom_density_ridges(aes(x=val, y=var, fill=season), alpha=0.3) +
  geom_vline(aes(xintercept=0), linetype = "dashed") +
  my.theme +
  xlab("Strength of effect on region classification (west to east)") +
  ylab("") +
  scale_fill_viridis_d(name="Season")

ggsave(plot.np, filename="Figs/Fig6NPMANOVA.jpeg", width = 8, height=8)


#7. SUMMARY STATS####

#7a. Dataset----
raw.raw <- read.csv("Data/LBCU_FilteredData_Segmented.csv") %>% 
  dplyr::filter(!id %in% c(46768277, 33088, 129945787, 46770723, 46769927))

locs.raw <- read.csv("Data/LBCUMCLocations.csv")

locs.n <-locs.raw  %>% 
  dplyr::filter(season %in% c("breed", "winter")) %>% 
  group_by(id, season) %>% 
  summarize(n=n()) %>% 
  group_by(id) %>% 
  summarize(n=n()) %>% 
  dplyr::filter(n==2) %>% 
  ungroup()

locs <- locs.raw %>% 
  dplyr::filter(id %in% locs.n$id)

raw <- raw.raw %>% 
  dplyr::filter(id %in% locs.n$id)

#Number of locations
nrow(raw)

#Number of individuals
length(unique(locs$id))

#Years of locations
nrow(raw)/365

#Years of tracking per individual
raw %>% 
  group_by(id) %>% 
  summarize(n=n()) %>% 
  ungroup() %>% 
  summarize(mean = mean(n),
            sd = sd(n),
            mean.yr = mean/365,
            sd.yr = sd/365,
            max = max(n),
            min = min(n),
            max.yr = max/365,
            min.yr = min/365)

#Number of summarized locations per season
table(locs$season)

#7b. Clustering----
clust.raw <- read.csv("Data/LBCUKDEClusters.csv")

clust.sum <- clust.raw %>% 
  group_by(id, season, nclust, kdecluster) %>% 
  summarize(n=n()) %>% 
  ungroup() %>% 
  dplyr::filter(n < 100) %>% 
  group_by(season, nclust, n) %>% 
  summarize(count=n()) %>% 
  ungroup()
unique(clust.sum$n)
clust.sum %>% 
  dplyr::filter(n==1,
                season=="breed")

#7c. MC----
mc.raw <- read.csv("Data/LBCUMigConnectivity.csv")

#across numbers of clusters
mc.raw %>% 
  group_by(season) %>% 
  summarize(MC = mean(MC),
            MClow = mean(MClow),
            MChigh = mean(MChigh))

#for each number of clusters
mc.raw %>% 
  group_by(season, nclust) %>% 
  summarize(MC = mean(MC),
            MClow = mean(MClow),
            MChigh = mean(MChigh))

#across seasons
mc.raw %>% 
  group_by(nclust) %>% 
  summarize(MC = mean(MC),
            MClow = mean(MClow),
            MChigh = mean(MChigh))

#7d. Trend----
clust <- read.csv("Data/LBCUBBSClusters.csv") 

load("~/Library/Application Support/bbsBayes/bbs_raw_data.RData")

bird <- bbs_data[["bird"]] %>% 
  mutate(id = paste(countrynum, statenum, Route, sep="-")) %>% 
  dplyr::filter(id %in% clust$id,
                AOU==2640) %>% 
  mutate(rt.uni=paste(statenum, Route, sep="-"),
         rt.uni.y=paste(rt.uni, Year, sep="-"))

route <- bbs_data[["route"]] %>% 
  mutate(id = paste(countrynum, statenum, Route, sep="-")) %>% 
  inner_join(clust) %>% 
  mutate(strat_name=svmcluster,
         rt.uni=paste(statenum, Route, sep="-"),
         rt.uni.y=paste(rt.uni, Year, sep="-")) %>% 
  left_join(bird) %>% 
  mutate(SpeciesTotal = ifelse(is.na(SpeciesTotal), 0, SpeciesTotal))

#Number of routes
length(unique(route$id))

#Years
summary(route$Year)

#Number of LBCU
summary(route$SpeciesTotal)
sd(route$SpeciesTotal)

#Routes per cluster
table(clust$svmcluster)
