library(tidyverse)
library(ggmap)
library(sf)

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
        axis.text = element_blank())

whemi <- map_data("world", region=c("Canada", 
                                    "USA", 
                                    "Mexico",
                                    "Guatemala", 
                                    "Belize", 
                                    "El Salvador",
                                    "Honduras", 
                                    "Nicaragua", 
                                    "Costa Rica",
                                    "Panama", 
                                    "Jamaica", 
                                    "Cuba", 
                                    "The Bahamas",
                                    "Haiti", 
                                    "Dominican Republic", 
                                    "Antigua and Barbuda",
                                    "Dominica", 
                                    "Saint Lucia", 
                                    "Saint Vincent and the Grenadines", 
                                    "Barbados",
                                    "Grenada",
                                    "Trinidad and Tobago",
                                    "Colombia",
                                    "Venezuela",
                                    "Guyana",
                                    "Suriname",
                                    "Ecuador",
                                    "Peru",
                                    "Brazil",
                                    "Bolivia",
                                    "Paraguay",
                                    "Chile",
                                    "Argentina",
                                    "Uruguay")) %>% 
  dplyr::filter(!group%in%c(258:264))

nam <- map_data("world", region=c("Canada", 
                                  "USA", 
                                  "Mexico",
                                  "Guatemala", 
                                  "Belize", 
                                  "El Salvador",
                                  "Honduras", 
                                  "Nicaragua", 
                                  "Costa Rica",
                                  "Panama", 
                                  "Jamaica", 
                                  "Cuba", 
                                  "The Bahamas",
                                  "Haiti", 
                                  "Dominican Republic", 
                                  "Antigua and Barbuda",
                                  "Dominica", 
                                  "Saint Lucia", 
                                  "Saint Vincent and the Grenadines", 
                                  "Barbados",
                                  "Grenada",
                                  "Trinidad and Tobago")) %>% 
  dplyr::filter(!group%in%c(258:264))

#1. Study area figure----

#2. Connectivity figure----
mn <- read.csv("Data/LBCUMCLocations.csv")

plot.bw <- ggplot(mn %>% dplyr::filter(season%in% c("breed", "winter"))) +
  geom_polygon(data=nam, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  geom_path(aes(x=lon, y=lat, group=id), colour="grey20") +
  geom_point(aes(x=lon, y=lat, fill=season), pch=21, colour="grey20", size=2.5, alpha = 0.7) +
  coord_sf(xlim=c(-170, -30), expand = FALSE, crs=4326) +
  scale_fill_manual(values=c("cyan4", "darkgoldenrod1")) +
  xlab("") +
  ylab("") +
  map.theme
plot.bw

ggsave(plot.bw, filename = "Figs/Breed2Winter.jpeg", width = 10, height = 8, device="jpeg")

plot.bf <- ggplot(mn %>% dplyr::filter(season %in% c("breed", "fallmig"))) +
  geom_polygon(data=nam, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  geom_path(aes(x=lon, y=lat, group=id), colour="grey20") +
  geom_point(aes(x=lon, y=lat, fill=season), pch=21, colour="grey20", size=2.5, alpha = 0.7) +
  coord_sf(xlim=c(-170, -30), expand = FALSE, crs=4326) +
  scale_fill_manual(values=c("cyan4", "darkred")) +
  xlab("") +
  ylab("") +
  map.theme
plot.bf

ggsave(plot.bf, filename = "Figs/Breed2Fall.jpeg", width = 10, height = 8, device="jpeg")

plot.fw <- ggplot(mn %>% dplyr::filter(season %in% c("winter", "fallmig"))) +
  geom_polygon(data=nam, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  geom_path(aes(x=lon, y=lat, group=id), colour="grey20") +
  geom_point(aes(x=lon, y=lat, fill=season), pch=21, colour="grey20", size=2.5, alpha = 0.7) +
  coord_sf(xlim=c(-170, -30), expand = FALSE, crs=4326) +
  scale_fill_manual(values=c("darkred", "darkgoldenrod1")) +
  xlab("") +
  ylab("") +
  map.theme
plot.fw

ggsave(plot.fw, filename = "Figs/Fall2Wint.jpeg", width = 10, height = 8, device="jpeg")

plot.ws <- ggplot(mn %>% dplyr::filter(season %in% c("winter", "springmig"))) +
  geom_polygon(data=nam, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  geom_path(aes(x=lon, y=lat, group=id), colour="grey20") +
  geom_point(aes(x=lon, y=lat, fill=season), pch=21, colour="grey20", size=2.5, alpha = 0.7) +
  coord_sf(xlim=c(-170, -30), expand = FALSE, crs=4326) +
  scale_fill_manual(values=c("darkolivegreen4", "darkgoldenrod1")) +
  xlab("") +
  ylab("") +
  map.theme
plot.ws

ggsave(plot.ws, filename = "Figs/Wint2Spring.jpeg", width = 10, height = 8, device="jpeg")

plot.sb <- ggplot(mn %>% dplyr::filter(season %in% c("breed", "springmig"))) +
  geom_polygon(data=nam, aes(x=long, y=lat, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  geom_path(aes(x=lon, y=lat, group=id), colour="grey20") +
  geom_point(aes(x=lon, y=lat, fill=season), pch=21, colour="grey20", size=2.5, alpha = 0.7) +
  coord_sf(xlim=c(-170, -30), expand = FALSE, crs=4326) +
  scale_fill_manual(values=c("cyan4", "darkolivegreen4")) +
  xlab("") +
  ylab("") +
  map.theme
plot.sb

ggsave(plot.sb, filename = "Figs/Spring2Breed.jpeg", width = 10, height = 8, device="jpeg")

#3. Moving window figure----


mantel$width <- factor(mantel$width, labels=c("100 km\nwindow width", "500 km\nwindow width", "1000 km\nwindow width", "2000 km\nwindow width", "3000 km\nwindow width", "4000 km\nwindow width"))

mantel$overlap <- factor(mantel$overlap, labels=c("25% window overlap", "50% window overlap", "75% window overlap"))

plot.r <- ggplot(mn.utm) +
  #  geom_polygon(data=nam.utm, aes(x=X, y=Y, group=group), colour = "gray85", fill = "gray75", size=0.3) +
  geom_rect(data=mantel,
            aes(xmin= start, xmax=end, ymin=2000000, ymax=8000000, fill=r),
            alpha = 0.8) +
  geom_path(aes(x=X, y=Y, group=id), colour="grey50") +
  geom_point(aes(x=X, y=Y, colour=season)) +
  facet_grid(width ~ overlap) +
  scale_fill_viridis_c(name = "Mantel\ncorrelation\ncoefficient") +
  scale_colour_manual(values=c("grey70", "grey30"), name="Season") +
  my.theme +
  theme(legend.position = "right",
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

ggsave(plot.r, filename = "Figs/MigratoryDivides_WithPts.jpeg", width = 12, height = 10, device="jpeg")