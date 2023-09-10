library(tidyverse)
library(ggmap)
library(sf)
library(gridExtra)
library(mapdata)
library(maptools)
library(ggpubr)
library(ggridges)
library(paletteer)
library(viridis)
library(ebirdst)
library(raster)

#TO DO: STANDARDIZE POINT BORDER COLOUR

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

#Get map data
# nam <- map_data("usa") %>% 
#   rbind(map_data("worldHires", "Canada")) %>% 
#   rbind(map_data("worldHires", "Mexico")) %>% 
#   mutate(group = paste0(region,"-", group),
#          subregion = ifelse(is.na(subregion), "Mainland", subregion)) %>%  
#   dplyr::filter(!str_sub(subregion, -6, -1)%in%c("Island", "island"),
#                 long < 0)

country <- map_data("world", region=c("Canada", 
                                      "USA",
                                      "Mexico"))

lake <- map_data("lakes")

#Get colours
seasons <- paletteer_d("nord::victory_bonds", n=4) 
seasons

groups <- viridis_pal(option="viridis")(n=5)

#1. Figure 1: Study area figure----
#Get data and wrangle for deployment locations----
dep <- read.csv("/Users/ellyknight/Documents/SMBC/Analysis/lbcu_exploration/Data/LBCUCleanedData.csv") %>% 
  group_by(id) %>% 
  dplyr::filter(row_number()==1) %>% 
  ungroup()

#make up study IDs
study <- data.frame(study = unique(dep$study)) %>% 
  arrange(study) %>% 
  mutate(studyid = c(6,4,9,7,2,3,5,1,8)) %>% 
  left_join(dep) %>% 
  mutate(longr = round(long/2)*2,
         latr = round(lat/2)*2) %>% 
  group_by(study, studyid, longr, latr) %>% 
  summarize(n=n()) %>% 
  ungroup()

#get eBird range
ebd <- load_raster("/Users/ellyknight/Library/Application Support/ebirdst/lobcur-ERD2019-STATUS-20200930-da308f90", product="abundance_seasonal")

#wrangle eBird range - this is slow
breed <- projectRaster(ebd$breeding, crs="+proj=longlat +datum=WGS84 +no_defs +type=crs")
winter <- projectRaster(ebd$nonbreeding, crs="+proj=longlat +datum=WGS84 +no_defs +type=crs")

breed.pt <- breed %>% 
  aggregate(10) %>% 
  rasterToPoints(., spatial = TRUE) %>% 
  data.frame() %>% 
  dplyr::filter(breeding > 0)

winter.pt <- winter %>% 
  aggregate(10) %>% 
  rasterToPoints(., spatial = TRUE) %>% 
  data.frame() %>% 
  dplyr::filter(nonbreeding > 0)

range.pt <- full_join(breed.pt, winter.pt) %>% 
  mutate(season = case_when(is.na(breeding) ~ "Nonbreeding",
                            is.na(nonbreeding) ~ "Breeding",
                            !is.na(breeding) & !is.na(nonbreeding) ~ "Year round"))

#plot
plot.sa <- ggplot() +
  geom_polygon(data=country, aes(x=long, y=lat, group=group), fill="gray80", colour = "gray90", size=0.3) +
  geom_polygon(data=lake, aes(x=long, y=lat, group=group), fill="white", colour = "white", size=0.3) +
  coord_sf(xlim=c(-129, -75), ylim=c(15, 57), expand = FALSE, crs=4326) +
  geom_tile(data = range.pt, aes(x = x, y = y, fill=season), alpha = 0.7) +
  geom_jitter(data=study, aes(x=longr, y=latr, size=n), pch=21, colour="grey90", fill="grey10", alpha = 0.7, width=1, height=1) +
  xlab("") +
  ylab("") +
  map.theme +
  theme(legend.position = "right") +
  scale_size(range = c(2, 7), name="Individuals\ntagged") +
  scale_fill_manual(values=c(as.character(seasons[c(1,4)]), "purple4"), name="Seasonal\nrange")
plot.sa

ggsave(plot.sa, filename="Figs/Fig1StudyArea.jpeg", width=8, height=6)

#2. Figure 2: All tracks figure----
raw <- read.csv("Data/LBCU_FilteredData_Segmented.csv")

#read in & wrangle data
loc <- read.csv("Data/LBCUMCLocations.csv") %>% 
  dplyr::filter(season=="winter") %>% 
  mutate(manual = case_when(X > -10960000 & distance < 100000 ~ 4,
                           lon < -105 & lon > -108 & distance < 10000 ~ 2,
                           lon < -108 & lon > -118 ~ 2,
                           lon < -118 ~ 1,
                           !is.na(lon) ~ 3)) %>% 
  dplyr::select(id, manual) %>% 
  unique() %>% 
  full_join(read.csv("Data/LBCUMCLocations.csv")) %>% 
  mutate(season = factor(season, levels=c("winter", "springmig", "breed", "fallmig"))) %>% 
  arrange(id, year, season, cluster)

#separate into spring and fall
loc.spring <- loc %>% 
  dplyr::filter(season!="fallmig") %>% 
  mutate(plot=paste0(id, "-", year))

loc.fall <- loc %>% 
  dplyr::filter(season!="springmig") %>% 
  mutate(plot=paste0(id, "-", year))

#define migratory divides
div <- data.frame(div=rep(c("1-2", "2-3", "3-4"), 2),
         type=c(rep("start", 3), rep("end", 3)),
         lat = c(55, 52, 54,
                 32, 22, 18),
         lon = c(-118, -113.2, -111,
                 -118, -111.2, -97),
         season = rep("season", 6))

#plot spring
plot.spring <- ggplot() +
  geom_polygon(data=country, aes(x=long, y=lat, group=group), fill="gray80", colour = "gray90", size=0.3) +
  geom_polygon(data=lake, aes(x=long, y=lat, group=group), fill="white", colour = "white", size=0.3) +
  coord_sf(xlim=c(-129, -75), ylim=c(15, 57), expand = FALSE, crs=4326) +
  geom_line(data=loc.spring, aes(x=lon, y=lat, group=plot), colour="grey30", alpha = 0.4) +
  geom_point(data=loc.spring, aes(x=lon, y=lat, fill=season), pch=21, colour="grey90", alpha = 0.7, size=3) +
  xlab("") +
  ylab("") +
  map.theme +
  geom_line(data=div, aes(x=lon, y=lat, group=div), linetype="dashed", size=1) +
  ggtitle("Prebreeding migration") +
  scale_fill_manual(values=as.character(seasons[c(4,3,1)]))
#plot.spring

#plot fall
plot.fall <- ggplot() +
  geom_polygon(data=country, aes(x=long, y=lat, group=group), fill="gray80", colour = "gray90", size=0.3) +
  geom_polygon(data=lake, aes(x=long, y=lat, group=group), fill="white", colour = "white", size=0.3) +
  coord_sf(xlim=c(-129, -75), ylim=c(15, 57), expand = FALSE, crs=4326) +
  geom_line(data=loc.fall, aes(x=lon, y=lat, group=plot), colour="grey30", alpha = 0.4) +
  geom_point(data=loc.fall, aes(x=lon, y=lat, fill=season), pch=21, colour="grey90", alpha = 0.7, size=3) +
  xlab("") +
  ylab("") +
  map.theme +
  geom_line(data=div, aes(x=lon, y=lat, group=div), linetype="dashed", size=1) +
  ggtitle("Postbreeding migration")  +
  scale_fill_manual(values=as.character(seasons[c(4,1,3)]))
#plot.fall

#make legend
plot.legend <- ggplot() +
  geom_point(data=sample_n(loc.fall, 10), aes(x=lon, y=lat, fill=season), pch=21, colour="grey90", alpha = 0.7, size=3) +
  geom_line(data=div, aes(x=lon, y=lat, group=div, linetype=season), size=1) +
  scale_fill_manual(values=as.character(seasons[c(4,1,3)]), name="Season", labels=c("Stationary\nnonbreeding",  "Breeding", "Migration\nstopover")) + 
  scale_linetype_manual(values=c("dashed"), labels=c("Hypothesized\nmigratory divide"), name="") +
  theme(legend.position = "bottom")
#plot.legend
legend <- cowplot::get_legend(plot.legend)

#put together and save
ggsave(grid.arrange(plot.spring, plot.fall, legend,
                    heights = c(4,1),
                    widths = c(4,4),
                    layout_matrix = rbind(c(1,2),
                                          c(3,3))), filename="Figs/Fig2Tracks.jpeg", width=8, height=6)

#3. Figure 3: Clustering figure----

#3a. Wrangle----
clust.raw <- read.csv("Data/LBCUKDEClusters.csv")

clust <- clust.raw %>% 
  group_by(id, season, nclust, group) %>%
  summarize(n=n()) %>% 
  group_by(id, season, nclust) %>% 
  dplyr::filter(n == max(n)) %>% 
  ungroup() %>% 
  left_join(clust.raw) %>% 
  group_by(id, season, nclust, group) %>%
  summarize(X=mean(X),
            Y=mean(Y)) 

clust.ll <- clust  %>% 
  dplyr::filter(!is.na(X)) %>% 
  st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
  st_transform(crs=4326) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  rename(lat = Y, lon = X) %>% 
  cbind(clust %>% 
          dplyr::filter(!is.na(X)))

#write.csv(dplyr::filter(clust.ll, nclust %in% c("3", "manual")), "Data/LBCUclusterlocs.csv", row.names=FALSE)

clust.ll$season <- factor(clust.ll$season, levels=c("breed", "fallmig", "winter", "springmig"),
                             labels=c("Breeding", "Postbreeding\nmigration\nstopover", "Nonbreeding", "Prebreeding\nmigration\nstopover"))
clust.ll$nclust <- factor(clust.ll$nclust, levels=c("2", "3", "4", "5"), labels=c("2 groups", "3 groups", "4 groups", "5 groups"))
  
#3b. Clustering plot----
groups <- viridis_pal(option="viridis")(n=5)

plot.clust <- ggplot(clust.ll) +
  geom_polygon(data=country, aes(x=long, y=lat, group=group), fill="gray80", colour = "gray90", size=0.3) +
  geom_polygon(data=lake, aes(x=long, y=lat, group=group), fill="white", colour = "white", size=0.3) +
  geom_point(aes(x=lon, y=lat, fill=factor(group)), size = 3, pch=21, colour="grey90", alpha=0.7) +
  coord_sf(xlim=c(min(clust.ll$lon)-5, max(clust.ll$lon)+5), ylim = c(min(clust.ll$lat)-5, max(clust.ll$lat)+5), expand = FALSE, crs=4326) +
  facet_grid(season ~ nclust) +
  xlab("") +
  ylab("") +
  map.theme +
  theme(legend.position = "bottom") +
  scale_fill_manual(values=groups[c(5,3,1,4,2)], name="Group") +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  theme(panel.spacing = unit(0.2, "lines"))
plot.clust

#3c. Connectivity wrangling----
results <- read.csv("Data/LBCUMigConnectivity.csv")

sum <- results %>% 
  dplyr::select(MC, MClow, MChigh, nclust, originseason, targetseason, boot) %>% 
  unique() %>% 
  group_by(originseason, targetseason, nclust) %>% 
  summarize(MC = mean(MC),
            MClow = mean(MClow),
            MChigh = mean(MChigh)) %>% 
  ungroup()

sum$targetseason <- factor(sum$targetseason, levels=c("breed", "fallmig", "winter", "springmig"),
                          labels=c("Breeding", "Postbreeding\nmigration\nstopover", "Nonbreeding", "Prebreeding\nmigration\nstopover"))
sum$originseason <- factor(sum$originseason, levels=c("breed","winter"),
                           labels=c("Breeding", "Nonbreeding"))

#3d. Connectivity plot----
plot.mc <- ggplot(sum) +
  geom_point(aes(x=nclust, y=MC)) +
  geom_errorbar(aes(x=nclust, ymin=MClow, ymax=MChigh)) +
  my.theme +
  ylab("Strength of migratory connectivity") +
  xlab("Number of groups") +
  ylim(c(0,1)) +
  facet_grid(targetseason ~ originseason) +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank()) +
  ggtitle("Reference season for migratory\nconnectivity estimation")
plot.mc

#3e. Put together----
ggsave(grid.arrange(plot.clust, plot.mc,
                    widths = c(3.5,6.5),
                    heights = c(1,20),
                    layout_matrix = rbind(c(2,NA),
                                          c(2,1))), filename="Figs/Fig3Cluster.jpeg", width=8, height=8)

#4. Trend----

#4a. Read in data----
trend.list <- read.csv("Data/LBCU_trend_gamye.csv") %>% 
  st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
  st_transform(crs=4326) %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  rename(lat = Y, lon = X) %>% 
  cbind(read.csv("Data/LBCU_trend_gamye.csv")) %>% 
  rename(group = knncluster) %>% 
  mutate(Region = case_when(group==1 ~ "Central",
                            group==2 ~ "West",
                            group==3 ~ "East"))

trend <- trend.list %>% 
  dplyr::select(Region, Trend, 'Trend_Q0.025', 'Trend_Q0.975') %>% 
  unique() %>% 
  rename(up = 'Trend_Q0.975', down = 'Trend_Q0.025')

trend$Region <- factor(trend$Region, levels=c("West", "Central", "East"))

#4b. BBS Regions----
clust <- read.csv("Data/LBCUKDEClusters.csv") %>% 
  dplyr::filter(season=="breed", nclust==3) %>% 
  dplyr::select(id, group, lat, lon) %>% 
  unique()  %>% 
  mutate(Region = case_when(group==1 ~ "Central",
                            group==2 ~ "West",
                            group==3 ~ "East")) %>% 
  left_join(trend %>% 
              dplyr::select(Region, Trend) %>% 
              unique())

plot.map <- ggplot(trend.list) +
  geom_polygon(data=country, aes(x=long, y=lat, group=group), fill="gray80", colour = "gray90", size=0.3) +
  geom_polygon(data=lake, aes(x=long, y=lat, group=group), fill="white", colour = "white", size=0.3) +
  geom_point(aes(x=lon, y=lat, colour=factor(Trend)), pch=21, fill="white", size=1.5) +
  geom_point(data=clust, aes(x=lon, y=lat, fill=factor(Trend)), colour="grey90", pch=21,  size=3, alpha = 0.7) +
  scale_colour_manual(values=groups[c(1,3,5)]) +
  scale_fill_manual(values=groups[c(1,3,5)]) +
  coord_sf(xlim=c(min(trend.list$lon)-5, max(trend.list$lon)+5), ylim = c(min(trend.list$lat)-5, max(trend.list$lat)+5), expand = FALSE, crs=3857) +
  map.theme +
  xlab("") +
  ylab("") +
#  geom_text(aes(x=-95, y=56, label="A"), size=10)
  theme(plot.margin = margin(t = 0.5, r = 1, b = 1, l = 0))

plot.map

#4c. Trend----
plot.bar <- ggplot(trend) +
  geom_hline(aes(yintercept=0), linetype="dashed", colour="grey30") +
  geom_errorbar(aes(x=Region, ymin=down, ymax=up)) +
  geom_point(aes(x=Region, y=Trend, fill=factor(Trend)), size=4, pch=21) +
  scale_fill_manual(values=groups[c(1,3,5)]) +
  ylab("Population trend") +
  my.theme +
  theme(legend.position="none",
        strip.background = element_blank(),
        strip.text.x = element_blank())
#  geom_text(aes(x=3.4, y=4.4, label="B"), size=10)
plot.bar

#4d. Trajectories----
indices <- read.csv("Data/LBCU_indices_gamye.csv") %>% 
  mutate(Region = case_when(Region==1 ~ "Central",
                            Region==2 ~ "West",
                            Region==3 ~ "East")) %>% 
  left_join(trend %>% 
              dplyr::select(Region, Trend) %>% 
              unique())

plot.index <- ggplot(indices) +
  geom_ribbon(aes(x=Year, ymin=Index_q_0.025, ymax=Index_q_0.975, group=factor(Region)), alpha = 0.2) +
  geom_line(aes(x=Year, y=Index, group=Region, colour=factor(Trend)), size=1) +
  geom_vline(aes(xintercept = 2007), linetype = "dashed") +
  scale_colour_manual(values=groups[c(1,3,5)], name="Population\ntrend") +
  ylab("Relative abundance index") +
  my.theme +
  theme(legend.position="none",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = margin(t=,1, r = 5))
#  geom_text(aes(x=2020, y=7.9, label="C"), size=10)
plot.index

#4e. Legend----
plot.legend <- ggplot(indices) +
  geom_point(data=trend.list, aes(x=lon, y=lat, colour=factor(Trend))) +
  geom_line(aes(x=Year, y=Index, group=Region, colour=factor(Trend))) +
  scale_colour_manual(values=groups[c(3,5,1)], name="Group", labels=c("West", "Central", "East")) +
  ylab("Relative abundance index") + 
  theme(legend.position = "bottom")

legend <- get_legend(plot.legend)

#4f. Put together----
ggsave(grid.arrange(plot.map, plot.bar, plot.index, legend,
                    widths = c(3,3,3), heights = c(2,0.3),
                    layout_matrix = rbind(c(1,2,3),
                                          c(NA,4,NA))),
       height=3.8, width=10, units='in', filename="Figs/Fig4Trend.jpeg")

#5. Behavioural attributes----

#5a. Get data----
dat.mig <- read.csv("Data/MigrationStopovers.csv") %>% 
  rename(Group = region) %>% 
  mutate(n = as.numeric(n)) %>% 
  mutate(n = case_when(Group=="west" ~ (n-0.2),
                       Group=="central" ~ n,
                       Group=="east" ~ (n + 0.2)))
dat.mig$season <- factor(dat.mig$season, levels=c("fallmig", "springmig"),
                         labels=c("Postbreeding\nmigration", "Prebreeding\nmigration"))
dat.mig$Group <- factor(dat.mig$Group, levels=c("west", "central", "east"), labels=c("West", "Central", "East")) 

dat.stop <- read.csv("Data/MigrationStopoverLength.csv") %>% 
  rename(Group = region) %>% 
  mutate(n = as.numeric(n)) %>% 
  mutate(n = case_when(Group=="west" ~ (n-0.2),
                       Group=="central" ~ n,
                       Group=="east" ~ (n + 0.2)))
dat.stop$season <- factor(dat.stop$season, levels=c("fallmig", "springmig"),
                         labels=c("Postbreeding\nmigration", "Prebreeding\nmigration"))
dat.stop$Group <- factor(dat.stop$Group, levels=c("west", "central", "east"), labels=c("West", "Central", "East")) 

dat.dur <- read.csv("Data/MigrationTiming.csv") %>% 
  rename(Group = region)
dat.dur$season <- factor(dat.dur$season, levels=c("fallmig", "springmig"),
                         labels=c("Postbreeding\nmigration", "Prebreeding\nmigration"))
dat.dur$Group <- factor(dat.dur$Group, levels=c("west", "central", "east"), labels=c("West", "Central", "East"))

dat.dist <- read.csv("Data/MigrationDistance.csv") %>% 
  rename(Group = region)
dat.dist$season <- factor(dat.dist$season, levels=c("fallmig", "springmig"),
                         labels=c("Postbreeding\nmigration", "Prebreeding\nmigration"))
dat.dist$Group <- factor(dat.dist$Group, levels=c("west", "central", "east"), labels=c("West", "Central", "East"))

dat.rate <- read.csv("Data/MigrationRate.csv") %>% 
  rename(Group = region)
dat.rate$season <- factor(dat.rate$season, levels=c("fallmig", "springmig"),
                          labels=c("Postbreeding\nmigration", "Prebreeding\nmigration"))
dat.rate$Group <- factor(dat.rate$Group, levels=c("west", "central", "east"), labels=c("West", "Central", "East"))

dat.wint <- read.csv("Data/WinterHRs.csv") %>% 
  rename(Group = region) %>% 
  mutate(n = as.numeric(n)) %>% 
  mutate(n = case_when(Group=="west" ~ (n-0.1),
                       Group=="central" ~ n,
                       Group=="east" ~ (n + 0.05)))
dat.wint$season <- factor(dat.wint$season, levels=c("fallmig", "springmig"),
                          labels=c("Postbreeding\nmigration", "Prebreeding\nmigration"))
dat.wint$Group <- factor(dat.wint$Group, levels=c("west", "central", "east"), labels=c("West", "Central", "East"))

dat.hr <- read.csv("Data/WinterHRSize.csv") %>% 
  rename(Group = region)
dat.hr$season <- factor(dat.hr$season, levels=c("breed", "fallmig", "winter", "springmig"),
                          labels=c("Breeding", "Postbreeding\nmigration", "Nonbreeding", "Prebreeding\nmigration"))
dat.hr$Group <- factor(dat.hr$Group, levels=c("west", "central", "east"), labels=c("West", "Central", "East"))

#5b. Departure----
plot.dep <- ggplot(dat.dur) +
  geom_density_ridges(aes(x=depart, y=season, fill=Group), alpha = 0.4) + 
  my.theme +
  xlab("Day of migration departure season") +
  ylab("") +
  scale_x_continuous(breaks=c(32, 91, 152, 213), labels=c("Feb", "Apr", "Jun", "Aug")) +
  scale_fill_manual(values=groups) +
  theme(legend.position = "none")
plot.dep

#5c. Arrival----
plot.arr <- ggplot(dat.dur) +
  geom_density_ridges(aes(x=arrive, y=season, fill=Group), alpha = 0.4) + 
  my.theme +
  xlab("Day of migration arrival season") +
  ylab("") +
  scale_x_continuous(breaks=c(91, 152, 213, 274), labels=c("Apr", "Jun", "Aug", "Oct")) +
  scale_fill_manual(values=groups) +
  theme(axis.text.y = element_blank()) +
  theme(legend.position = "none")
plot.arr

#5d. Duration----
plot.dur <- ggplot(dat.dur) +
  geom_density_ridges(aes(x=duration, y=season, fill=Group), alpha = 0.4) + 
  my.theme +
  xlab("Days of migration") +
  ylab("") +
  scale_fill_manual(values=groups) +
  theme(axis.text.y = element_blank()) +
  theme(legend.position = "none")
plot.dur

#5e. Stopovers----
plot.stop <- ggplot(dat.mig) +
  geom_density_ridges(aes(x=n, y=season, fill=Group), alpha = 0.4, stat="binline") + 
  my.theme +
  xlab("Number of migration stopovers") +
  ylab("") +
  scale_fill_manual(values=groups) +
  theme(axis.text.y = element_blank()) +
  theme(legend.position = "none")
plot.stop

#5f. Stopover duration----
plot.stopdur <- ggplot(dat.stop) +
  geom_density_ridges(aes(x=n, y=season, fill=Group), alpha = 0.4) + 
  my.theme +
  xlab("Total days of migration stopover") +
  ylab("") +
  scale_fill_manual(values=groups) +
  theme(axis.text.y = element_blank()) +
  theme(legend.position = "none")
plot.stopdur

#5f. Distance----
plot.dist <- ggplot(dat.dist) +
  geom_density_ridges(aes(x=dist, y=season, fill=Group), alpha = 0.4) + 
  my.theme +
  xlab("Migration distance (km)") +
  ylab("") +
  scale_fill_manual(values=groups) +
  theme(axis.text.y = element_blank()) +
  theme(legend.position = "none")
plot.dist

#5g. Rate----
plot.rate <- ggplot(dat.rate) +
  geom_density_ridges(aes(x=rate, y=season, fill=Group), alpha = 0.4) + 
  my.theme +
  xlab("Migration rate (km/day)") +
  ylab("") +
  scale_fill_manual(values=groups) +
  theme(legend.position = "none")
plot.rate

#5h. Range area----
plot.hr <- ggplot(dat.hr) +
  geom_density_ridges(aes(x=log(area), y=season, fill=Group), alpha = 0.4) + 
  my.theme +
  xlab("Natural log of use area") +
  ylab("") +
  scale_fill_manual(values=groups) +
  theme(legend.position = "none")
plot.hr

#5i. Number of wintering ranges----
plot.wint <- ggplot(dat.wint) +
  geom_density_ridges(aes(x=n, y=season, fill=Group), alpha = 0.4, stat="binline") + 
  my.theme +
  xlab("Number of wintering home ranges") +
  ylab("") +
  scale_fill_manual(values=groups) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(breaks=c(1,2,3), labels=c(1,2,3)) +
  theme(legend.position = "none")
plot.wint

#5j. Legend----
plot.legend <- ggplot(dat.wint) +
  geom_density_ridges(aes(x=n, y=season, fill=Group), alpha = 0.4, stat="binline") +
  scale_fill_manual(values=groups) +
  theme(legend.position = "right")
plot.legend
legend <- cowplot::get_legend(plot.legend)

#5h. Put together----
plot.behav <- grid.arrange(plot.dep, plot.arr, plot.dur, plot.dist, plot.rate, plot.stopdur, plot.stop, plot.hr, plot.wint, legend,
                           widths = c(4,3,3,3), heights=c(3,3,3),
                           layout_matrix=rbind(c(1,2,3,4),
                                               c(5,6,7,NA),
                                               c(8,9,10,NA)))

ggsave(plot.behav, filename="Figs/Fig5Behave.jpeg", width = 12, height = 9)

#6. Environmental attributes----

#6a. Wrangle mrpp deltas----
mrpp <- read.csv("Data/MRPP.csv") %>% 
  group_by(season, region) %>% 
  summarize(delta.mn = mean(delta),
            delta.sd = sd(delta),
            delta.up = quantile(delta, 0.975),
            delta.lw = quantile(delta, 0.025)) %>% 
  ungroup()

mrpp$region <- factor(mrpp$region, levels=c("west", "central", "east"), labels=c("West", "Central", "East"))
mrpp$season <- factor(mrpp$season, levels=c("breed", "fallmig", "winter", "springmig"),
                      labels=c("Breeding home range", "Postbreeding migration\nstopover", "Nonbreeding home range", "Prebreeding migration\nstopover"))

#6b. Wrangle nmds scores----
scores <- read.csv("Data/NMDSScores.csv") %>% 
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

scores$region <- factor(scores$region, levels=c("west", "central", "east"), labels=c("West", "Central", "East"))
scores$season <- factor(scores$season, levels=c("breed", "fallmig", "winter", "springmig"),
                        labels=c("Breeding home range", "Postbreeding migration\nstopover", "Nonbreeding home range", "Prebreeding migration\nstopover"))

covscores <- read.csv("Data/NMDSCovscores.csv") %>% 
  group_by(season, cov) %>% 
  summarize(x.mn = mean(NMDS1),
            x.up = quantile(NMDS1, 0.975),
            x.lw = quantile(NMDS1, 0.025),
            y.mn = mean(NMDS2),
            y.up = quantile(NMDS2, 0.975),
            y.lw = quantile(NMDS2, 0.025)) %>% 
  ungroup()

covscores$season <- factor(covscores$season, levels=c("breed", "fallmig", "winter", "springmig"),
                      labels=c("Breeding home range", "Postbreeding migration\nstopover", "Nonbreeding home range", "Prebreeding migration\nstopover"))
covscores$cov <- factor(covscores$cov, levels=c("flooded_vegetation", "seasonality", "recurrence", "change_norm", "grass", "conv", "crops", "built", "pdsi"),
                       labels=c("Wetland", "Water\nseasonality", "Water\nrecurrence", "Water\nchange", "Grassland", "Grassland\nconversion", "Crop", "Urban", "Drought"))

#6c. Plot----
plot.mrpp <- ggplot(mrpp) +
  geom_errorbar(aes(x=region, ymin=delta.lw, ymax=delta.up, colour=region)) +
  geom_point(aes(x=region, y=delta.mn, colour=region)) +
  facet_grid(season~.) +
  scale_colour_manual(values=groups) +
  my.theme +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank()) + 
  ylab("Within-group MRPP mean")
#plot.mrpp

plot.nmds <- ggplot() +
  geom_point(data=scores, aes(x=x.mn, y=y.mn, fill=region), pch=21, size=2) +
  geom_segment(data=covscores, aes(x=0, y=0, xend=x.mn, yend=y.mn, colour=cov), arrow = arrow(length = unit(0.2, "cm"))) +
  facet_wrap(.~season, scales="free", ncol=1, strip.position = "right") +
  scale_fill_manual(values=groups, name="Group") +
  scale_colour_paletteer_d("ggthemes_ptol::qualitative", 9, direction=1, name="Environmental\nattribute") +
  my.theme + 
  ylab("") +
  xlab("") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank())
#plot.nmds

plot.covs <- grid.arrange(plot.mrpp, plot.nmds, widths = c(2, 3), heights =c(4),
             layout_matrix=rbind(c(1,2)))

ggsave(plot.covs, filename="Figs/Fig6NMDSMRPP.jpeg", width = 8, height=10)


#8. SUMMARY STATS####

#8a. Dataset----
raw.raw <- read.csv("Data/LBCU_FilteredData_Segmented.csv") %>% 
  dplyr::filter(!id %in% c(46768277, 33088, 129945787, 46770723, 46769927, 86872))

locs.raw <- read.csv("Data/LBCUMCLocations.csv") %>% 
  dplyr::filter(!id %in% c(46768277, 33088, 129945787, 46770723, 46769927, 86872))

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

#Histogram of tracking years
years <- raw %>% 
  group_by(id) %>% 
  summarize(n=n()/365) %>% 
  ungroup()
  

#Number of summarized locations per season
table(locs$season)

#Number of individuals with spring & fall mig stopovers
locs %>% 
  dplyr::filter(season%in%c("springmig", "fallmig")) %>% 
  dplyr::select(season, id) %>% 
  unique() %>% 
  group_by(season) %>% 
  summarize(n=n()) %>% 
  ungroup()

#8b. Clustering----
clust.raw <- read.csv("Data/LBCUKDEClusters.csv")

clust.ind <- clust.raw %>% 
  group_by(id, season, nclust, group) %>% 
  summarize(n=n()) %>% 
  ungroup() %>% 
  dplyr::filter(n < 100) %>% 
  dplyr::select(-season) %>% 
  unique()

clust.sum <- clust.ind %>% 
  group_by(id, nclust) %>% 
  summarize(groups=n()) %>% 
  ungroup()

#number of individuals with # of clusters
table(clust.sum$nclust, clust.sum$groups)

#8c. MC----
mc.raw <- read.csv("Data/LBCUMigConnectivity.csv")

mc.raw %>% 
  summarize(MC = mean(MC), 
            MClow = mean(MClow),
            MChigh = mean(MChigh))

#by strategy
mc.raw %>% 
  mutate(strategy = ifelse(nclust=="manual", "manual", "kmeans")) %>% 
  group_by(strategy) %>% 
  summarize(MC = mean(MC), 
            MClow = mean(MClow),
            MChigh = mean(MChigh))

#across seasons
mc.raw %>% 
  group_by(originseason, targetseason) %>% 
  summarize(MC = mean(MC),
            MClow = mean(MClow),
            MChigh = mean(MChigh))

#for each number of clusters by season
mc.raw %>% 
  group_by(originseason, targetseason, nclust) %>% 
  summarize(MC = mean(MC),
            MClow = mean(MClow),
            MChigh = mean(MChigh)) %>% 
  ungroup() %>% 
  arrange(originseason, targetseason, -MC) %>% 
  View()

#across number of clusters
mc.raw %>% 
  group_by(nclust) %>% 
  summarize(MC = mean(MC),
            MClow = mean(MClow),
            MChigh = mean(MChigh))

#8d. Trend----
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
  dplyr::filter(Year >= 1970) %>% 
  mutate(strat_name=knncluster,
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
table(clust$knncluster)

#8e. Movement chars----
dat.mig <- read.csv("Data/MigrationStopovers.csv") 
dat.dur <- read.csv("Data/MigrationTiming.csv")
dat.dist <- read.csv("Data/MigrationDistance.csv")
dat.rate <- read.csv("Data/MigrationRate.csv")
dat.wint <- read.csv("Data/WinterHRs.csv")
dat.hr <- read.csv("Data/WinterHRSize.csv")

#Distance
dat.dist %>% 
  group_by(region) %>% 
  summarize(dist= mean(dist))

dat.dist %>% 
  group_by(region, season) %>% 
  summarize(dist= mean(dist))

dat.dist %>% 
  group_by(season) %>% 
  summarize(dist= mean(dist))

#Duration
dat.dur %>% 
  group_by(region) %>% 
  summarize(dist= mean(duration))

dat.dur %>% 
  group_by(season) %>% 
  summarize(dist= mean(duration))

dat.dur %>% 
  group_by(region, season) %>% 
  summarize(dist= mean(duration))

#Departure
dat.dur %>% 
  group_by(region) %>% 
  summarize(depart= mean(depart))

dat.dur %>% 
  group_by(season, region) %>% 
  summarize(depart= mean(depart))

#Arrival
dat.dur %>% 
  group_by(region) %>% 
  summarize(arrive= mean(arrive))

dat.dur %>% 
  group_by(season, region) %>% 
  summarize(arrive= mean(arrive))

#Stopovers
dat.mig %>% 
  group_by(region) %>% 
  summarize(n= mean(n))

dat.mig %>% 
  group_by(region, season) %>% 
  summarize(n= mean(n))

dat.mig %>% 
  group_by(season) %>% 
  summarize(n= mean(n))

#Rate
dat.rate %>% 
  group_by(region) %>% 
  summarize(rate= mean(rate))

dat.rate %>% 
  group_by(season, region) %>% 
  summarize(rate= mean(rate))

dat.rate %>% 
  group_by(season) %>% 
  summarize(rate= mean(rate))

#Areas
dat.hr %>% 
  group_by(region) %>% 
  summarize(area= mean(area))

dat.hr %>% 
  group_by(season, region) %>% 
  summarize(area= mean(area))

dat.hr %>% 
  group_by(season) %>% 
  dplyr::filter(region!="west") %>% 
  summarize(area= mean(area),
            sd = sd(area, na.rm=TRUE))

#Distance between wintering hrs
kd.wint <- read_sf("gis/shp/kde_individual.shp") %>% 
  dplyr::filter(season=="winter") %>% 
  st_make_valid() %>% 
  st_centroid() %>% 
  st_transform(crs=3857) %>% 
  mutate(id = paste0(bird, "-", year, "", season)) %>% 
  arrange(id, cluster)

kd.xy <- st_coordinates(kd.wint)

traj <- as.ltraj(xy=kd.xy[,c("X", "Y")],
                 id=kd.wint$id,
                 date=kd.wint$cluster,
                 typeII=FALSE,
                 proj4string = CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"))
names(traj) <- unique(kd.wint$id)

dat.traj <- rbindlist(traj, idcol=TRUE) %>% 
  rename(id='.id', X=x, Y=y) %>% 
  mutate(id=as.numeric(id)) %>% 
  arrange(id, date) %>% 
  dplyr::select(dist) %>% 
  cbind(kd.wint) %>% 
  dplyr::filter(!is.na(dist)) %>% 
  mutate(dist = dist/1000) %>% 
  mutate(region = case_when(kdecluster==2 ~ "West", 
                            kdecluster==1 ~ "Central",
                            kdecluster==3 ~ "East"))

ggplot(dat.traj) +
  geom_histogram(aes(x=dist, fill=region)) +
  xlab("Distance between\nsuccessive winter home ranges")

ggsave(filename="Figs/WinteringHomeRangeDistances.jpeg", width = 3, height = 3)

#8f. Included individuals----
locs <- read.csv("Data/LBCUMCLocations.csv")
dat <- read.csv("Data/LBCU_FilteredData_Segmented.csv")

locs.n <-locs.raw  %>% 
  dplyr::filter(season %in% c("breed", "winter")) %>% 
  group_by(id, season) %>% 
  summarize(n=n()) %>% 
  group_by(id) %>% 
  summarize(n=n()) %>% 
  dplyr::filter(n==2) %>% 
  ungroup()

individuals <- locs.n %>% 
  dplyr::select(id) %>% 
  unique() %>% 
  left_join(dat %>% 
              mutate(date = lubridate::ymd_hms(date)) %>% 
              group_by(study, id, sensor, sex) %>% 
              summarize(start = min(date),
                        end = max(date)) %>% 
              ungroup())

write.csv(individuals, "Data/LBCUIndividualsIncluded.csv", row.names = FALSE)

#APPENDICES####

#1. Appendix 1: Deployment metadata----
birds <- read.csv("Data/LBCUMCLocations.csv") %>% 
  group_by(study, id, year) %>% 
  summarize(days=sum(days)) %>% 
  group_by(study, id) %>% 
  summarize(startyear = min(year),
            endyear = max(year),
           days = sum(days)) %>% 
  ungroup() %>% 
  mutate(use = ifelse(id %in% c(46768277, 33088, 129945787, 46770723, 46769927, 86872), "No", "Yes"), 
         studyid = case_when(study=="BC" ~ 6,
                             study=="IW" ~ 4,
                             study=="Jay" ~ 9,
                             study=="MN" ~ 7,
                             study=="MX" ~ 2,
                             study=="NB" ~ 3,
                             study=="TX" ~ 5,
                             study=="USGS" ~ 1,
                             study=="WY" ~ 8)) %>% 
  dplyr::select(studyid, id, startyear, endyear, days, use) %>% 
  arrange(studyid, id)

write.csv(birds, "Data/LBCUIndividualMetadata.csv", row.names = FALSE)
