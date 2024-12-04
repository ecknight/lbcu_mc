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
library(terra)

setwd("G:/My Drive/SMBC")

my.theme <- theme_classic() +
  theme(text=element_text(size=12),
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
  theme(text=element_text(size=12),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12)) +
  theme(plot.margin = margin(t = 1, r = 1, b = 1, l = 1))

country <- map_data("world", region=c("Canada", 
                                      "USA",
                                      "Mexico"))

lake <- map_data("lakes")

#Get colours
seasons <- paletteer_d("nord::victory_bonds", n=4) 
seasons

groups <- palette.colors(palette = "Dark 2")[c(1:3, 6, 4)]

ggplot() + 
  geom_point(aes(x=c(1:length(groups)), y=c(1:length(groups)), colour=factor(c(1:length(groups)))), size=5, show.legend=FALSE) +
  scale_colour_manual(values=as.character(groups))

cols <- groups[c(2,1,2,1,4,3,3)]

#1. Figure 1: Study area figure----
#Get data and wrangle for first location per tag
dat <- read.csv("Data/LBCUCleanedData.csv") |> 
  group_by(study, id) |> 
  dplyr::filter(row_number()==1) |> 
  ungroup()
  
#Get study ids
study <- data.frame(study = unique(dat$study)) |> 
  arrange(study) |> 
  mutate(studyid = c(6,4,9,7,2,3,5,1,8))

#Get deployment data and match to nearest location for sample size per deployment location
dep.raw <- read.csv("Data/DeploymentLocations.csv")|> 
  st_as_sf(coords=c("long", "lat"), crs=4326, remove=FALSE)

dep <- dat |> 
  mutate(deprow = dat |> 
           st_as_sf(coords=c("long", "lat"), crs=4326) |> 
           st_nearest_feature(dep.raw)) |> 
  dplyr::select(-long, -lat) |> 
  left_join(read.csv("Data/DeploymentLocations.csv") |> 
              mutate(deprow = row_number())) |> 
  group_by(studyid, long, lat) |> 
  summarize(n=n()) |> 
  ungroup()

#get eBird range
# ebirdst_download_status(species = "lobcur",
#                         download_abundance = TRUE,
#                         download_ranges = FALSE,
#                         force=TRUE)

ebd <- load_raster("lobcur", product="abundance", period="seasonal", resolution = "3km")

#wrangle eBird range - this is slow
breed <- raster(project(ebd$breeding, "EPSG:4326"))
winter <- raster(project(ebd$nonbreeding, "EPSG:4326"))

breed.pt <- breed |> 
  aggregate(10) |> 
  rasterToPoints(., spatial = TRUE) |> 
  data.frame() |> 
  dplyr::filter(breeding > 0)

winter.pt <- winter |> 
  aggregate(10) |> 
  rasterToPoints(., spatial = TRUE) |> 
  data.frame() |> 
  dplyr::filter(nonbreeding > 0)

range.pt <- full_join(breed.pt, winter.pt) |> 
  mutate(season = case_when(is.na(breeding) ~ "Nonbreeding",
                            is.na(nonbreeding) ~ "Breeding",
                            !is.na(breeding) & !is.na(nonbreeding) ~ "Year round"))

#plot
plot.sa <- ggplot() +
  geom_polygon(data=country, aes(x=long, y=lat, group=group), fill="gray80", colour = "gray90", size=0.3) +
  geom_polygon(data=lake, aes(x=long, y=lat, group=group), fill="white", colour = "white", size=0.3) +
  coord_sf(xlim=c(-129, -75), ylim=c(15, 57), expand = FALSE, crs=4326) +
  geom_tile(data = range.pt, aes(x = x, y = y, fill=season), alpha = 1, show.legend = FALSE) +
  geom_point(data=dep, aes(x=long, y=lat, size=n), pch=21, colour="black", fill="white", alpha = 0.7) +
  xlab("") +
  ylab("") +
  map.theme +
  theme(legend.position = "none") +
  ggtitle("Tag deployments")  +
  scale_size(range = c(2, 7), name="Individuals\ntagged") +
  scale_fill_manual(values=c(as.character(seasons[c(1,4)]), "grey50"), name="Seasonal\nrange") + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size=12))
plot.sa

#Get tracking data
raw <- read.csv("Data/LBCU_FilteredData_Segmented.csv")

#read in & wrangle data
loc <- read.csv("Data/LBCUMCLocations.csv") |> 
  dplyr::filter(season=="winter") |> 
  mutate(manual = case_when(X > -10960000 & distance < 100000 ~ 4,
                           lon < -105 & lon > -108 & distance < 10000 ~ 2,
                           lon < -108 & lon > -118 ~ 2,
                           lon < -118 ~ 1,
                           !is.na(lon) ~ 3)) |> 
  dplyr::select(id, manual) |> 
  unique() |> 
  full_join(read.csv("Data/LBCUMCLocations.csv")) |> 
  mutate(season = factor(season, levels=c("winter", "springmig", "breed", "fallmig"))) |> 
  arrange(id, year, season, cluster)

#separate into spring and fall
loc.spring <- loc |> 
  dplyr::filter(season!="fallmig") |> 
  mutate(plot=paste0(id, "-", year))

loc.fall <- loc |> 
  dplyr::filter(season!="springmig") |> 
  mutate(plot=paste0(id, "-", year))

#plot spring
plot.spring <- ggplot() +
  geom_polygon(data=country, aes(x=long, y=lat, group=group), fill="gray90", colour = "gray70", linewidth=0.3) +
  geom_polygon(data=lake, aes(x=long, y=lat, group=group), fill="white", colour = "grey70", linewidth=0.3) +
  coord_sf(xlim=c(-129, -75), ylim=c(14.6, 57), expand = FALSE, crs=4326) +
  geom_line(data=loc.spring, aes(x=lon, y=lat, group=plot), colour="grey30", alpha = 0.2) +
  geom_point(data=loc.spring, aes(x=lon, y=lat, fill=season), pch=21, colour="grey90", alpha = 0.7, size=3) +
  map.theme +
  ggtitle("Prebreeding migration") +
  scale_fill_manual(values=as.character(seasons[c(4,3,1)])) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size=12))
plot.spring

#plot fall
plot.fall <- ggplot() +
  geom_polygon(data=country, aes(x=long, y=lat, group=group), fill="gray90", colour = "gray70", linewidth=0.3) +
  geom_polygon(data=lake, aes(x=long, y=lat, group=group), fill="white", colour = "grey70", linewidth=0.3) +
  coord_sf(xlim=c(-129, -75), ylim=c(14.6, 57), expand = FALSE, crs=4326) +
  geom_line(data=loc.fall, aes(x=lon, y=lat, group=plot), colour="grey30", alpha = 0.2) +
  geom_point(data=loc.fall, aes(x=lon, y=lat, fill=season), pch=21, colour="grey90", alpha = 0.7, size=3) +
  map.theme +
  ggtitle("Postbreeding migration")  +
  scale_fill_manual(values=as.character(seasons[c(4,1,3)])) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size=12))
plot.fall

#make legend

plot.legend <- ggplot(legend.dat) +
  geom_point(aes(x=lon, y=lat, fill=season), pch=21, colour="grey90", alpha = 0.7, size=3) +
  geom_point(data=dep, aes(x=long, y=lat, size=n), pch=21, colour="black", fill="white", alpha = 0.7) +
  scale_size(range = c(2, 7), name="Individuals\ntagged") +
  scale_fill_manual(values=as.character(c(seasons[c(4,1,3)], "grey50")), name="Season", labels=c("Stationary nonbreeding",  "Breeding", "Migration stopover", "Year round")) + 
  theme(legend.position = "right",
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))
plot.legend
legend <- get_legend(plot.legend)

#put together and save


ggsave(grid.arrange(plot.sa, legend, plot.spring, plot.fall,
                    heights = c(4,4),
                    widths = c(4,4),
                    layout_matrix = rbind(c(1,2),
                                          c(3,4))),
       filename="Figs/Fig1SA&Tracks.jpeg", width=6, height=6.5)

#2. Figure 2: Clustering figure----

#2a. Wrangle----
#Fix the group #s to match
clust <- read.csv("Data/LBCUKDEClusters.csv") |> 
  mutate(group = case_when(nclust=="flyway" & group==1 ~ 2,
                           nclust=="flyway" & group==2 ~ 1,
                           nclust=="expert" & group==1 ~ 2,
                           nclust=="expert" & group==2 ~ 1,
                           nclust=="expert" & group==3 ~ 4,
                           nclust=="expert" & group==4 ~ 3,
                           !is.na(group) ~ group)) |> 
  dplyr::select(-lat, -lon)

clust.ll <- clust  |> 
  dplyr::filter(!is.na(X)) |> 
  st_as_sf(coords=c("X", "Y"), crs=3857) |> 
  st_transform(crs=4326) |> 
  st_coordinates() |> 
  data.frame() |> 
  rename(lat = Y, lon = X) |> 
  cbind(clust |> 
          dplyr::filter(!is.na(X)))

#write.csv(dplyr::filter(clust.ll, nclust %in% c("3", "manual")), "Data/LBCUclusterlocs.csv", row.names=FALSE)

clust.ll$season <- factor(clust.ll$season, levels=c("breed", "fallmig", "winter", "springmig"),
                             labels=c("Breeding", "Postbreeding\nmigration\nstopover", "Nonbreeding", "Prebreeding\nmigration\nstopover"))
clust.ll$nclust <- factor(clust.ll$nclust, levels=c("flyway", "expert", "2", "3", "4", "5"), labels=c("Flyway groups", "Expert groups", "2 groups", "3 groups", "4 groups", "5 groups"))
  
#2b. Clustering plot----
plot.clust <- ggplot(clust.ll) +
  geom_polygon(data=country, aes(x=long, y=lat, group=group), fill="gray90", colour = "gray70", size=0.3) +
  # geom_polygon(data=lake, aes(x=long, y=lat, group=group), fill="white", colour = "gray70", size=0.3) +
  geom_point(aes(x=lon, y=lat, fill=factor(group)), size = 3, pch=21, colour="grey90", alpha=0.7) +
  coord_sf(xlim=c(min(clust.ll$lon)-5, max(clust.ll$lon)+5), ylim = c(min(clust.ll$lat)-5, max(clust.ll$lat)+5), expand = FALSE, crs=4326) +
  facet_grid(season ~ nclust) +
  xlab("") +
  ylab("") +
  map.theme +
  theme(legend.position = "bottom") +
  scale_fill_manual(values=groups, name="Group") +
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
  theme(panel.spacing = unit(0.2, "lines"),
        legend.margin=margin(-10,0,0,0))
plot.clust

#2c. Connectivity wrangling----
sum <- read.csv("Data/LBCUMigConnectivity.csv")

sum$targetseason <- factor(sum$targetseason, levels=c("breed", "fallmig", "winter", "springmig"),
                          labels=c("Breeding", "Postbreeding\nmigration\nstopover", "Nonbreeding", "Prebreeding\nmigration\nstopover"))
sum$originseason <- factor(sum$originseason, levels=c("breed","winter"),
                           labels=c("Breeding", "Nonbreeding"))
sum$nclust <- factor(sum$nclust, levels=c("flyway", "expert", "2", "3", "4", "5"), labels=c("Flyway", "Expert", "2", "3", "4", "5"))

#2d. Connectivity plot----
plot.mc <- ggplot(sum) +
  geom_point(aes(x=nclust, y=MC)) +
  geom_errorbar(aes(x=nclust, ymin=MClow, ymax=MChigh)) +
  my.theme +
  ylab("Strength of migratory connectivity") +
  xlab("Number of groups") +
  ylim(c(0,1)) +
  facet_grid(targetseason ~ originseason) +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank(),
        strip.text.x = element_text(size=10),
        axis.text.x=element_text(size=10, angle = 45, vjust = 1, hjust = 1),
        axis.text.y=element_text(size=10),
        legend.title=element_text(size=12)) +
  ggtitle("Reference season for migratory\nconnectivity estimation")
plot.mc

#2e. Put together----
ggsave(grid.arrange(plot.clust, plot.mc,
                    widths = c(3.5,6.5),
                    heights = c(1.5,20),
                    layout_matrix = rbind(c(2,NA),
                                          c(2,1))), filename="Figs/Fig2Cluster.jpeg", width=10, height=7)

#3. Trend----

#3a. Read in data----
trend.list <- read.csv("Data/LBCU_trend_gamye.csv") |> 
  st_as_sf(coords=c("X", "Y"), crs=3857) |> 
  st_transform(crs=4326) |> 
  st_coordinates() |> 
  data.frame() |> 
  rename(lat = Y, lon = X) |> 
  cbind(read.csv("Data/LBCU_trend_gamye.csv")) |> 
  rename(group = knncluster) |> 
  mutate(Region = case_when(nclust=="3" & group==1 ~ "Intermountain",
                            nclust=="3" & group==2 ~ "West",
                            nclust=="3" & group==3 ~ "Plains",
                            nclust=="flyway" & group==1 ~ "Pacific",
                            nclust=="flyway" & group==2 ~ "Midcontinent",
                            nclust=="expert" & group==1 ~ "West",
                            nclust=="expert" & group==2 ~ "Intermountain",
                            nclust=="expert" & group==3 ~ "Chihuahuan",
                            nclust=="expert" & group==4 ~ "Gulf"),
         Region = factor(Region, levels = c("Pacific", "West", "Intermountain", "Midcontinent", "Chihuahuan", "Gulf", "Plains")))  |> 
  mutate(nclust = factor(nclust, levels = c("flyway", "expert", "3"),
                         labels = c("Flyway", "Expert", "Clustering")))


trend <- trend.list |> 
  dplyr::select(nclust, Region, Trend, 'Trend_Q0.025', 'Trend_Q0.975') |> 
  unique() |> 
  rename(up = 'Trend_Q0.975', down = 'Trend_Q0.025') |> 
  mutate(wd = case_when(nclust=="Clustering" ~ 0.45,
                        nclust=="Expert" ~ 0.6,
                        nclust=="Flyway" ~ 0.3))

#trend$Region <- factor(trend$Region, levels=c("West", "Intermountain", "Plains"))

#3b. BBS Regions----
clust.raw <- read.csv("Data/LBCUKDEClusters.csv") |> 
  dplyr::filter(season=="breed", nclust %in% c("3", "expert", "flyway"))

clust <-  clust.raw |> 
  group_by(nclust, id, group) |>
  summarize(n=n()) |> 
  group_by(id) |> 
  dplyr::filter(n == max(n)) |> 
  ungroup() |> 
  left_join(clust.raw) |> 
  group_by(nclust, id, group) |>
  summarize(lat=mean(lat),
            lon=mean(lon))  |> 
  ungroup() |> 
  mutate(Region = case_when(nclust=="3" & group==1 ~ "Intermountain",
                            nclust=="3" & group==2 ~ "West",
                            nclust=="3" & group==3 ~ "Plains",
                            nclust=="flyway" & group==1 ~ "Pacific",
                            nclust=="flyway" & group==2 ~ "Midcontinent",
                            nclust=="expert" & group==1 ~ "West",
                            nclust=="expert" & group==2 ~ "Intermountain",
                            nclust=="expert" & group==3 ~ "Chihuahuan",
                            nclust=="expert" & group==4 ~ "Gulf"),
         Region = factor(Region, levels = c("Pacific", "West", "Intermountain", "Midcontinent", "Chihuahuan", "Gulf", "Plains")))  |> 
  mutate(nclust = factor(nclust, levels = c("flyway", "expert", "3"),
                         labels = c("Flyway", "Expert", "Clustering"))) |> 
  left_join(trend |> 
              dplyr::select(nclust, Region, Trend) |> 
              unique())

plot.map <- ggplot(trend.list) +
  geom_polygon(data=country, aes(x=long, y=lat, group=group), fill="gray90", colour = "gray70", size=0.3) +
  geom_polygon(data=lake, aes(x=long, y=lat, group=group), fill="white", colour = "gray70", size=0.3) +
  geom_point(aes(x=lon, y=lat, colour=Trend), pch=21, fill="white", size=1.5) +
  geom_point(data=clust, aes(x=lon, y=lat, fill=Trend), colour="grey90", pch=21,  size=3, alpha = 0.7) +
  scale_colour_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  coord_sf(xlim=c(min(trend.list$lon)-5, max(trend.list$lon)+5), ylim = c(min(trend.list$lat)-5, max(trend.list$lat)+5), expand = FALSE, crs=3857) +
  map.theme +
  xlab("") +
  ylab("") +
  theme(plot.margin = margin(t = 0.5, r = 0, b = -3, l = 0),
        strip.text.x = element_text(size=12)) + 
  facet_wrap(~nclust)
plot.map

#3c. Trend----
plot.bar <- ggplot(trend) +
  geom_hline(aes(yintercept=0), linetype="dashed", colour="grey30") +
  geom_errorbar(aes(x=Region, ymin=down, ymax=up, width = wd)) +
  geom_point(aes(x=Region, y=Trend, fill=Trend), size=4, pch=21) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  ylab("Population trend") +
  my.theme +
  xlab("") +
  theme(legend.position="bottom",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.y=element_text(vjust=-1.5)) + 
  facet_wrap(~nclust, scales = "free_x")
plot.bar

#3f. Put together----
ggsave(grid.arrange(plot.map, plot.bar,
                    widths = c(5), heights = c(2,3),
                    layout_matrix = rbind(c(1),
                                          c(2))),
       height=7, width=8, units='in', filename="Figs/Fig3Trend.jpeg")

#4. Behavioural attributes----

#migration departure and arrival date, migration duration, migration distance, migration rate, number of migration stopovers, migration stopover duration, number of wintering home ranges, and core stopover or home range size

#4a. Get data----
dat <- read.csv("Data/MovementBehaviours.csv") |> 
  mutate(season = factor(season, levels=c("breed", "fallmig", "winter", "springmig"),
                         labels=c("Breeding", "Postbreeding\nmigration", "Nonbreeding", "Prebreeding\nmigration"))) |> 
  mutate(Region = case_when(nclust=="3" & group==1 ~ "Intermountain",
                            nclust=="3" & group==2 ~ "West",
                            nclust=="3" & group==3 ~ "Plains",
                            nclust=="flyway" & group==1 ~ "Pacific",
                            nclust=="flyway" & group==2 ~ "Midcontinent",
                            nclust=="expert" & group==1 ~ "West",
                            nclust=="expert" & group==2 ~ "Intermountain",
                            nclust=="expert" & group==3 ~ "Chihuahuan",
                            nclust=="expert" & group==4 ~ "Gulf"),
         Region = factor(Region, levels = c("Pacific", "Midcontinent", "West", "Intermountain",  "Chihuahuan", "Gulf", "Plains"))) |> 
  mutate(nclust = factor(nclust, levels = c("flyway", "expert", "3"),
                         labels = c("Flyway", "Expert", "Clustering")))

#4b. Departure----
plot.dep <- ggplot(dat |> dplyr::filter(var=="depart")) +
  geom_density_ridges(aes(x=val, y=season, fill=Region), alpha = 0.5, scale=1) + 
  my.theme +
  xlab("Day of migration departure season") +
  ylab("") +
  scale_x_continuous(breaks=c(32, 91, 152, 213), labels=c("Feb", "Apr", "Jun", "Aug")) +
  scale_fill_manual(values=cols) +
  theme(legend.position = "none") +
  theme(plot.margin = margin(r=10, b=10),
        axis.title.x.bottom = element_text(vjust=5)) +
  facet_wrap(~nclust)
plot.dep

#4c. Arrival----
plot.arr <- ggplot(dat |> dplyr::filter(var=="arrive")) +
  geom_density_ridges(aes(x=val, y=season, fill=Region), alpha = 0.5, scale=1) + 
  my.theme +
  xlab("Day of migration arrival season") +
  ylab("") +
  scale_x_continuous(breaks=c(91, 152, 213, 274), labels=c("Apr", "Jun", "Aug", "Oct")) +
  scale_fill_manual(values=cols) +
  theme(legend.position = "none") +
  theme(plot.margin = margin(r=10, b=10),
        axis.title.x.bottom = element_text(vjust=5),
        strip.text = element_blank()) +
  facet_wrap(~nclust)
plot.arr

#4d. Rate----
plot.rate <- ggplot(dat |> dplyr::filter(var=="rate")) +
  geom_density_ridges(aes(x=val/100, y=season, fill=Region), alpha = 0.5, scale=1) + 
  my.theme +
  xlab("Migration rate (100 km/day)") +
  ylab("") +
  scale_fill_manual(values=cols) +
  theme(legend.position = "none") +
  theme(plot.margin = margin(r=10, b=10),
        axis.title.x.bottom = element_text(vjust=5),
        strip.text = element_blank()) +
  facet_wrap(~nclust)
plot.rate

#4e. Stopovers----
plot.stop <- ggplot(dat |> dplyr::filter(var=="Stopovers")) +
  geom_density_ridges(aes(x=val, y=season, fill=Region), alpha = 0.5, stat="binline", position=position_dodge(width=1), scale=1) + 
  my.theme +
  xlab("Number of migration stopovers") +
  ylab("") +
  scale_fill_manual(values=cols) +
  theme(legend.position = "none") +
  theme(plot.margin = margin(r=10, b=10),
        axis.title.x.bottom = element_text(vjust=5),
        strip.text = element_blank()) +
  facet_wrap(~nclust)
plot.stop

#4f. Stopover duration----
plot.stopdur <- ggplot(dat |> dplyr::filter(var=="stopoverduration")) +
  geom_density_ridges(aes(x=val, y=season, fill=Region), alpha = 0.5, scale=1) + 
  my.theme +
  xlab("Total days of migration stopover") +
  ylab("") +
  scale_fill_manual(values=cols) +
  theme(legend.position = "none") +
  theme(plot.margin = margin(r=10, b=10),
        axis.title.x.bottom = element_text(vjust=5),
        strip.text = element_blank()) +
  facet_wrap(~nclust)
plot.stopdur

#4g. Number of wintering ranges----
plot.wint <- ggplot(dat |> dplyr::filter(var=="WinterHRs")) +
  geom_density_ridges(aes(x=val, y=season, fill=Region), alpha = 0.5, stat="binline", position=position_dodge(width=0.8), scale=1) + 
  my.theme +
  xlab("Number of wintering home ranges") +
  ylab("") +
  scale_fill_manual(values=cols) +
  scale_x_continuous(breaks=c(1,2,3), labels=c(1,2,3)) +
  theme(legend.position = "none") +
  theme(plot.margin = margin(r=10, b=10),
        axis.title.x.bottom = element_text(vjust=5),
        strip.text = element_blank()) +
  facet_wrap(~nclust)
plot.wint

#4h. Range area----
plot.hr <- ggplot(dat |> dplyr::filter(var=="HRarea")) +
  geom_density_ridges(aes(x=log(val), y=season, fill=Region), alpha = 0.5, scale=1) + 
  my.theme +
  xlab(bquote(Natural~log~of~home~range~or~stopover~area~(km^2))) +
  ylab("") +
  scale_fill_manual(values=cols) +
  theme(legend.position = "none") +
  theme(plot.margin = margin(r=10, b=10),
        axis.title.x.bottom = element_text(vjust=5),
        strip.text = element_blank()) +
  facet_wrap(~nclust)
plot.hr

#4i. Legend----
plot.legend.flyway <- ggplot(dat |> dplyr::filter(nclust=="Flyway")) +
  geom_density_ridges(aes(x=val, y=season, fill=Region), alpha = 0.5, stat="binline", scale=1) +
  scale_fill_manual(values=cols[c(1,2)], name="Flyway\nregions") +
  theme(legend.position = "right",
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.key.size = unit(2,"line"),
        strip.text = element_blank())
plot.legend.flyway
legend.flyway <- get_legend(plot.legend.flyway)

plot.legend.expert <- ggplot(dat |> dplyr::filter(nclust=="Expert")) +
  geom_density_ridges(aes(x=val, y=season, fill=Region), alpha = 0.5, stat="binline", scale=1) +
  scale_fill_manual(values=cols[c(3:6)], name="Expert\nregions") +
  theme(legend.position = "right",
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.key.size = unit(2,"line"),
        strip.text = element_blank())
plot.legend.expert
legend.expert <- get_legend(plot.legend.expert)

plot.legend.cluster <- ggplot(dat |> dplyr::filter(nclust=="Clustering")) +
  geom_density_ridges(aes(x=val, y=season, fill=Region), alpha = 0.5, stat="binline", scale=1) +
  scale_fill_manual(values=cols[c(3,4,7)], name="Clustering\nregions") +
  theme(legend.position = "right",
        legend.text=element_text(size=12),
        legend.title=element_text(size=12),
        legend.key.size = unit(2,"line"),
        strip.text = element_blank())
plot.legend.cluster
legend.cluster <- get_legend(plot.legend.cluster)
 
#4h. Put together----
plot.behav <- grid.arrange(plot.dep, plot.arr, plot.hr, plot.wint,
                           legend.flyway, legend.expert, legend.cluster,
                           widths = c(6,2), heights=c(3.5,3,4,2),
                           layout_matrix=rbind(c(1,8),
                                               c(2,9),
                                               c(3,10),
                                               c(4,NA)))

ggsave(plot.behav, filename="Figs/Fig4Behave.jpeg", width = 10, height = 10)

#5. Habitat selection----

#5a. Read in predictions----
pred <- read.csv("Results/RSFPredictions.csv") |> 
  mutate(season = factor(season, levels=c("breed", "fallmig", "winter", "springmig"),
                         labels=c("Breeding", "Postbreeding\nmigration\nstopover", "Nonbreeding", "Prebreeding\nmigration\nstopover"))) |> 
  mutate(Region = case_when(nclust=="3" & group==1 ~ "Intermountain",
                            nclust=="3" & group==2 ~ "West",
                            nclust=="3" & group==3 ~ "Plains",
                            nclust=="flyway" & group==1 ~ "Pacific",
                            nclust=="flyway" & group==2 ~ "Midcontinent",
                            nclust=="expert" & group==1 ~ "West",
                            nclust=="expert" & group==2 ~ "Intermountain",
                            nclust=="expert" & group==3 ~ "Chihuahuan",
                            nclust=="expert" & group==4 ~ "Gulf",
                            group==0 ~  "No group\ndifference"),
         Region = factor(Region, levels = c("Pacific", "Midcontinent", "West", "Intermountain",  "Chihuahuan", "Gulf", "Plains", "No group\ndifference"))) |> 
  mutate(nclust = factor(nclust, levels = c("flyway", "expert", "3"),
                         labels = c("Flyway", "Expert", "Clustering"))) |> 
  mutate(cov = factor(cov, levels=c("crop", "grass", "wetland"),
                      labels=c("Cropland", "Grassland", "Wetland"))) |> 
  dplyr::filter(!is.na(Region))

#5b. Plot----
plot.covsflyway <- ggplot(pred |> dplyr::filter(nclust=="Flyway")) +
  geom_ribbon(aes(x=val, ymin=low, ymax=up, group=Region), fill="grey30", alpha = 0.1) +
  geom_line(aes(x=val, y=fit, colour=Region), linewidth=1) +
  facet_grid(season~cov, scales="free") +
  my.theme +
  scale_colour_manual(values=c(cols[c(1,2)], "grey50"), name="Flyway\nregions") +
  xlab("") +
  ylab("Marginal relative probability of selection") +
  ggtitle("Flyway") +
  theme(axis.text.x = element_text(size=9),
        axis.text.y = element_text(size=9),
        legend.text = element_text(size=12),
        strip.text.x = element_text(size=12),
        strip.text.y = element_blank(),
        legend.position = "bottom")
plot.covsflyway

plot.covsexpert <- ggplot(pred |> dplyr::filter(nclust=="Expert")) +
  geom_ribbon(aes(x=val, ymin=low, ymax=up, group=Region), fill="grey30", alpha = 0.1) +
  geom_line(aes(x=val, y=fit, colour=Region), linewidth=1) +
  facet_grid(season~cov, scales="free") +
  my.theme +
  scale_colour_manual(values=c(cols[c(3:6)], "grey50"), name="Expert\nregions") +
  xlab("Attribute value") +
  ylab("Marginal relative probability of selection") +
  ggtitle("Expert") +
  theme(axis.text.x = element_text(size=9),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.text = element_text(size=12),
        strip.text.x = element_text(size=12),
        strip.text.y = element_blank(),
        legend.position = "bottom")
plot.covsexpert

plot.covs3 <- ggplot(pred |> dplyr::filter(nclust=="Clustering")) +
  geom_ribbon(aes(x=val, ymin=low, ymax=up, group=Region), fill="grey30", alpha = 0.1) +
  geom_line(aes(x=val, y=fit, colour=Region), linewidth=1) +
  facet_grid(season~cov, scales="free") +
  my.theme +
  scale_colour_manual(values=c(cols[c(3,4,7)], "grey50"), name="Clustering\nregions") +
  xlab("") +
  ylab("Marginal relative probability of selection") +
  ggtitle("Clustering") +
  theme(axis.text.x = element_text(size=9),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.text = element_text(size=12),
        strip.text.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        legend.position = "bottom")
plot.covs3

plot.covs <- grid.arrange(plot.covsflyway, plot.covsexpert, plot.covs3, widths = c(6,5,6))

ggsave(plot.covs, filename="Figs/Fig5RSF.jpeg", width = 12, height=5)

#8. SUMMARY STATS####

#8a. Dataset----
raw.raw <- read.csv("Data/LBCU_FilteredData_Segmented.csv") |> 
  dplyr::filter(!id %in% c(46768277, 33088, 129945787, 46770723, 46769927, 86872))

locs.raw <- read.csv("Data/LBCUMCLocations.csv") |> 
  dplyr::filter(!id %in% c(46768277, 33088, 129945787, 46770723, 46769927, 86872))

locs.n <-locs.raw  |> 
  dplyr::filter(season %in% c("breed", "winter")) |> 
  group_by(id, season) |> 
  summarize(n=n()) |> 
  group_by(id) |> 
  summarize(n=n()) |> 
  dplyr::filter(n==2) |> 
  ungroup()

locs <- locs.raw |> 
  dplyr::filter(id %in% locs.n$id)

raw <- raw.raw |> 
  dplyr::filter(id %in% locs.n$id)

#Number of locations
nrow(raw)

#Number of individuals
length(unique(locs$id))

#Years of locations
nrow(raw)/365

#Years of tracking per individual
raw |> 
  group_by(id) |> 
  summarize(n=n()) |> 
  ungroup() |> 
  summarize(mean = mean(n),
            sd = sd(n),
            mean.yr = mean/365,
            sd.yr = sd/365,
            max = max(n),
            min = min(n),
            max.yr = max/365,
            min.yr = min/365)

#Histogram of tracking years
years <- raw |> 
  group_by(id) |> 
  summarize(n=n()/365) |> 
  ungroup()
  

#Number of summarized locations per season
table(locs$season)

#Number of individuals with spring & fall mig stopovers
locs |> 
  dplyr::filter(season%in%c("springmig", "fallmig")) |> 
  dplyr::select(season, id) |> 
  unique() |> 
  group_by(season) |> 
  summarize(n=n()) |> 
  ungroup()

#8b. Clustering----
clust.raw <- read.csv("Data/LBCUKDEClusters.csv")

table(clust.raw$nclust, clust.raw$group)

#8c. MC----
mc.raw <- read.csv("Data/LBCUMigConnectivity.csv")

mc.raw |> 
  summarize(MC = mean(MC), 
            MClow = mean(MClow),
            MChigh = mean(MChigh))

#by strategy
mc.raw |> 
  mutate(strategy = ifelse(nclust %in% c("expert", "flyway"), nclust, "kmeans")) |> 
  group_by(strategy) |> 
  summarize(MC = mean(MC), 
            MClow = mean(MClow),
            MChigh = mean(MChigh))

#across seasons
mc.raw |> 
  group_by(originseason, targetseason) |> 
  summarize(MC = mean(MC),
            MClow = mean(MClow),
            MChigh = mean(MChigh))

#for each number of clusters by season
mc.raw |> 
  group_by(originseason, targetseason, nclust) |> 
  summarize(MC = mean(MC),
            MClow = mean(MClow),
            MChigh = mean(MChigh)) |> 
  ungroup() |> 
  arrange(originseason, targetseason, -MC) |> 
  View()

#across number of clusters
mc.raw |> 
  group_by(nclust) |> 
  summarize(MC = mean(MC),
            MClow = mean(MClow),
            MChigh = mean(MChigh))

#8d. Trend----
clust <- read.csv("Data/LBCUBBSClusters.csv") 

load("~/Library/Application Support/bbsBayes/bbs_raw_data.RData")

bird <- bbs_data[["bird"]] |> 
  mutate(id = paste(countrynum, statenum, Route, sep="-")) |> 
  dplyr::filter(id %in% clust$id,
                AOU==2640) |> 
  mutate(rt.uni=paste(statenum, Route, sep="-"),
         rt.uni.y=paste(rt.uni, Year, sep="-"))

route <- bbs_data[["route"]] |> 
  mutate(id = paste(countrynum, statenum, Route, sep="-")) |> 
  inner_join(clust) |> 
  dplyr::filter(Year >= 1970) |> 
  mutate(strat_name=knncluster,
         rt.uni=paste(statenum, Route, sep="-"),
         rt.uni.y=paste(rt.uni, Year, sep="-")) |> 
  left_join(bird) |> 
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

aic.out <- read.csv("Results/MigrationCharsAIC.csv")

aic.choose <- aic.out |> 
  dplyr::filter(delta < 2) |> 
  group_by(nclust, mod) |> 
  mutate(mindf = min(df)) |> 
  ungroup() |> 
  dplyr::filter(df==mindf) |> 
  dplyr::select(nclust, mod, model, delta) |> 
  pivot_wider(names_from = nclust, values_from=c(model, delta))

dat.all <- read.csv("Data/MovementBehaviours.csv") |> 
  mutate(Region = case_when(nclust=="3" & group==1 ~ "Intermountain",
                            nclust=="3" & group==2 ~ "West",
                            nclust=="3" & group==3 ~ "Plains",
                            nclust=="flyway" & group==1 ~ "Pacific",
                            nclust=="flyway" & group==2 ~ "Midcontinent",
                            nclust=="expert" & group==1 ~ "West",
                            nclust=="expert" & group==2 ~ "Intermountain",
                            nclust=="expert" & group==3 ~ "Chihuahuan",
                            nclust=="expert" & group==4 ~ "Gulf"),
         Region = factor(Region, levels = c("Pacific", "Midcontinent", "West", "Intermountain",  "Chihuahuan", "Gulf", "Plains")))

#Departure
dat.all |>
  dplyr::filter(var=="depart2") |> 
  group_by(season, nclust, Region) |> 
  summarize(depart= mean(val))

#Arrival
dat.all |>
  dplyr::filter(var=="arrive2") |> 
  group_by(season, nclust, Region) |> 
  summarize(depart= mean(val))

#Stopovers
dat.all |>
  dplyr::filter(var=="Stopovers") |> 
  group_by(season, nclust, Region) |> 
  summarize(depart= mean(val))

#Stopover days
dat.all |>
  dplyr::filter(var=="stopoverduration") |> 
  group_by(season) |> 
  summarize(stopovers= mean(val))

#Rate
dat.all |>
  dplyr::filter(var=="rate") |> 
  group_by(season, nclust, Region) |> 
  summarize(depart= mean(val))

#Areas
dat.all |>
  dplyr::filter(var=="HRarea") |> 
  group_by(season, nclust, Region) |> 
  summarize(depart= mean(val)) |> 
  View()

#number of wintering home ranges
dat.all |>
  dplyr::filter(var=="WinterHRs") |> 
  group_by(season, nclust, Region) |> 
  summarize(depart= mean(val))


#Distance between wintering hrs
kd.wint <- read_sf("gis/shp/kde_individual.shp") |> 
  dplyr::filter(season=="winter") |> 
  st_make_valid() |> 
  st_centroid() |> 
  st_transform(crs=3857) |> 
  mutate(id = paste0(bird, "-", year, "", season)) |> 
  arrange(id, cluster)

kd.xy <- st_coordinates(kd.wint)

traj <- as.ltraj(xy=kd.xy[,c("X", "Y")],
                 id=kd.wint$id,
                 date=kd.wint$cluster,
                 typeII=FALSE,
                 proj4string = CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"))
names(traj) <- unique(kd.wint$id)

dat.traj <- rbindlist(traj, idcol=TRUE) |> 
  rename(id='.id', X=x, Y=y) |> 
  mutate(id=as.numeric(id)) |> 
  arrange(id, date) |> 
  dplyr::select(dist) |> 
  cbind(kd.wint) |> 
  dplyr::filter(!is.na(dist)) |> 
  mutate(dist = dist/1000) |> 
  mutate(region = case_when(kdecluster==2 ~ "West", 
                            kdecluster==1 ~ "Intermountain",
                            kdecluster==3 ~ "Plains"))

ggplot(dat.traj) +
  geom_histogram(aes(x=dist, fill=region)) +
  xlab("Distance between\nsuccessive winter home ranges")

ggsave(filename="Figs/WinteringHomeRangeDistances.jpeg", width = 3, height = 3)

#8f. Included individuals----
locs <- read.csv("Data/LBCUMCLocations.csv")
dat <- read.csv("Data/LBCU_FilteredData_Segmented.csv")

locs.n <-locs.raw  |> 
  dplyr::filter(season %in% c("breed", "winter")) |> 
  group_by(id, season) |> 
  summarize(n=n()) |> 
  group_by(id) |> 
  summarize(n=n()) |> 
  dplyr::filter(n==2) |> 
  ungroup()

individuals <- locs.n |> 
  dplyr::select(id) |> 
  unique() |> 
  left_join(dat |> 
              mutate(date = lubridate::ymd_hms(date)) |> 
              group_by(study, id, sensor, sex) |> 
              summarize(start = min(date),
                        end = max(date)) |> 
              ungroup())

write.csv(individuals, "Data/LBCUIndividualsIncluded.csv", row.names = FALSE)

#APPENDICES####

#1. Appendix 1: Deployment metadata----
birds <- read.csv("Data/LBCUMCLocations.csv") |> 
  group_by(study, id, year) |> 
  summarize(days=sum(days)) |> 
  group_by(study, id) |> 
  summarize(startyear = min(year),
            endyear = max(year),
           days = sum(days)) |> 
  ungroup() |> 
  mutate(use = ifelse(id %in% c(46768277, 33088, 129945787, 46770723, 46769927, 86872), "No", "Yes"), 
         studyid = case_when(study=="BC" ~ 6,
                             study=="IW" ~ 4,
                             study=="Jay" ~ 9,
                             study=="MN" ~ 7,
                             study=="MX" ~ 2,
                             study=="NB" ~ 3,
                             study=="TX" ~ 5,
                             study=="USGS" ~ 1,
                             study=="WY" ~ 8)) |> 
  dplyr::select(studyid, id, startyear, endyear, days, use) |> 
  arrange(studyid, id) |> 
  left_join(read.csv("Data/LBCU_FilteredData_Segmented.csv") |> 
                       dplyr::select(id, sex) |> 
                       unique())

write.csv(birds, "Data/LBCUIndividualMetadata.csv", row.names = FALSE)

lee <- read.csv("Data/lbcu_id_for_Elly.csv")

#2. Appendix 6 - AIC results----

#2a. Movement----
dredged1 <- read.csv("Results/MigrationCharsAIC.csv")

format1 <- dredged1 |> 
  dplyr::select(nclust, mod, model, df, logLik, AICc, delta, weight) |> 
  mutate(mod = factor(mod, levels=c("departure", "arrival", "rate", "stopover", "stopoverduration", "winterhrs", "area"), labels = c("migration departure timing", "migration arrival timing", "migration rate", "number of migration stopovers", "seasonal stopover duration", "number of nonbreeding home ranges", "natural log use of area")),
         nclust = factor(nclust, levels = c("flyway", "expert", "3"))) |> 
  arrange(nclust, mod, delta)

write.csv(format1, "Results/MigrationCharsAIC_formatted.csv", row.names = FALSE)


#2b.RSF----
dredged2 <- read.csv("Results/RSFAIC.csv") |> 
  dplyr::select(nclust, season, model, df, logLik, AICc, delta, weight) |> 
  mutate(season = factor(season, levels=c("breed", "fallmig", "winter", "springmig"), labels = c("breeding", "postbreeding migration", "nonbreeding", "prebreeding migration")),
         nclust = factor(nclust, levels = c("flyway", "expert", "3"))) |> 
  arrange(nclust, season, delta)

format2 <- dredged2 |> 
  mutate(m1 = case_when(!is.na(crop.group) ~ "crop*group",
                        is.na(crop.group) & !is.na(crop) ~ "crop",
                        is.na(crop) ~ NA),
         m2 = case_when(!is.na(grass.group) ~ "grass*group",
                        is.na(grass.group) & !is.na(grass) ~ "grass",
                        is.na(grass) ~ NA),
         m3 = case_when(!is.na(group.wetland) ~ "wetland*group",
                        is.na(group.wetland) & !is.na(wetland) ~ "wetland",
                        is.na(wetland) ~ NA)) |> 
  mutate(model = paste(m1, m2, m3, sep = " + ")) |> 
  rowwise() |> 
  mutate(model = str_replace_all(string = model, pattern = "NA [+] ", replacement = ""),
         model = str_replace_all(string = model, pattern = "[+] NA", replacement = ""),
         model = ifelse(model=="NA", "1", model)) |> 
  dplyr::select(nclust, season, model, df, logLik, AICc, delta, weight) |> 
  mutate(season = factor(season, levels=c("breed", "fallmig", "winter", "springmig"), labels = c("breeding", "postbreeding migration", "nonbreeding", "prebreeding migration")),
         nclust = factor(nclust, levels = c("flyway", "expert", "3"))) |> 
  arrange(nclust, season, delta) |> 
  group_by(nclust, season) |> 
  dplyr::filter(row_number() %in% c(1:5)) |> 
  ungroup()

write.csv(format2, "Results/RSFAIC_formatted.csv", row.names = FALSE)

#GRAPHICAL ABSTRACT###########

#1. Get tracking data----
loc <- read.csv("Data/LBCUMCLocations.csv") |> 
  dplyr::filter(season=="winter") |> 
  mutate(manual = case_when(X > -10960000 & distance < 100000 ~ 4,
                            lon < -105 & lon > -108 & distance < 10000 ~ 2,
                            lon < -108 & lon > -118 ~ 2,
                            lon < -118 ~ 1,
                            !is.na(lon) ~ 3)) |> 
  dplyr::select(id, manual) |> 
  unique() |> 
  full_join(read.csv("Data/LBCUMCLocations.csv")) |> 
  mutate(season = factor(season, levels=c("winter", "springmig", "breed", "fallmig"))) |> 
  dplyr::filter(season %in% c("breed", "winter")) |> 
  arrange(id, year, season, cluster) 

#2. Get trend data & wrangle----
trend.list <- read.csv("Data/LBCU_trend_gamye.csv") |> 
  st_as_sf(coords=c("X", "Y"), crs=3857) |> 
  st_transform(crs=4326) |> 
  st_coordinates() |> 
  data.frame() |> 
  rename(lat = Y, lon = X) |> 
  cbind(read.csv("Data/LBCU_trend_gamye.csv")) |> 
  rename(group = knncluster) |> 
  mutate(Region = case_when(nclust=="3" & group==1 ~ "Intermountain",
                            nclust=="3" & group==2 ~ "West",
                            nclust=="3" & group==3 ~ "Plains",
                            nclust=="flyway" & group==1 ~ "Pacific",
                            nclust=="flyway" & group==2 ~ "Midcontinent",
                            nclust=="expert" & group==1 ~ "West",
                            nclust=="expert" & group==2 ~ "Intermountain",
                            nclust=="expert" & group==3 ~ "Chihuahuan",
                            nclust=="expert" & group==4 ~ "Gulf"),
         Region = factor(Region, levels = c("Pacific", "West", "Intermountain", "Midcontinent", "Chihuahuan", "Gulf", "Plains")))

trend <- trend.list |> 
  dplyr::select(nclust, Region, Trend, 'Trend_Q0.025', 'Trend_Q0.975') |> 
  unique() |> 
  rename(up = 'Trend_Q0.975', down = 'Trend_Q0.025') |> 
  mutate(wd = case_when(nclust=="Clustering" ~ 0.45,
                        nclust=="Expert" ~ 0.6,
                        nclust=="Flyway" ~ 0.3))

clust.raw <- read.csv("Data/LBCUKDEClusters.csv")

clust <-  clust.raw |> 
  dplyr::filter(season=="breed", nclust %in% c("3", "expert", "flyway")) |> 
  group_by(nclust, id, group) |>
  summarize(n=n()) |> 
  group_by(id) |> 
  dplyr::filter(n == max(n)) |> 
  ungroup() |> 
  left_join(clust.raw, multiple="all") |> 
  group_by(nclust, id, group) |>
  summarize(lat=mean(lat),
            lon=mean(lon))  |> 
  ungroup() |> 
  mutate(Region = case_when(nclust=="3" & group==1 ~ "Intermountain",
                            nclust=="3" & group==2 ~ "West",
                            nclust=="3" & group==3 ~ "Plains",
                            nclust=="flyway" & group==1 ~ "Pacific",
                            nclust=="flyway" & group==2 ~ "Midcontinent",
                            nclust=="expert" & group==1 ~ "West",
                            nclust=="expert" & group==2 ~ "Intermountain",
                            nclust=="expert" & group==3 ~ "Chihuahuan",
                            nclust=="expert" & group==4 ~ "Gulf"),
         Region = factor(Region, levels = c("Pacific", "West", "Intermountain", "Midcontinent", "Chihuahuan", "Gulf", "Plains")))  |> 
  left_join(trend |> 
              dplyr::select(nclust, Region, Trend) |> 
              unique()) |> 
  mutate(nclust = factor(nclust, levels = c("flyway", "expert", "3"),
                         labels = c("Based on migratory flyway\n(2 groups)", "From expert knowledge\n(4 groups)", "Unsupervised clustering of\nmovement paths (3 groups)")))

#3. Put together----
all <- inner_join(loc, clust |> 
                    dplyr::select(-lat, -lon),
                  multiple="all")

#4. Plot maps----
ggplot(all) + 
  geom_polygon(data=country, aes(x=long, y=lat, group=group), fill="gray90", colour = "gray70", linewidth=0.3) +
  geom_polygon(data=lake, aes(x=long, y=lat, group=group), fill="white", colour = "grey70", linewidth=0.3) +
  coord_sf(xlim=c(-128, -92), ylim=c(19, 57), expand = FALSE, crs=4326) +
  geom_line(data=all, aes(x=lon, y=lat, group=id), colour="grey30", alpha = 0.2) +
  geom_point(data=all, aes(x=lon, y=lat, fill=Trend), pch=21, colour="grey30", size=3) +
  map.theme +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size=14, face="bold"), 
        strip.text = element_text(size=12),
        legend.position = "bottom",
        panel.spacing = unit(0.5, "lines")) +
#  ggtitle("Approaches for delineating of groups for annual cycle management of long-billed curlew") +
  facet_wrap(~nclust) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, name="Population\ntrend")

ggsave(filename="Figs/GraphicalAbstractFig_Legend.jpeg", width=10, height=4.5)
