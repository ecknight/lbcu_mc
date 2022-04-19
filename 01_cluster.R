library(tidyverse)
library(dtwclust)

#1. Load data----
dat.raw <- read.csv("Data/LBCU_FilteredData_Segmented.csv") %>% 
  dplyr::filter(!id %in% c(46768277, 33088)) %>% 
  mutate(legid = paste(id, year, sep="-"))

#2. Pick tag years with enough data----
ggplot(dat.raw) +
  geom_point(aes(y=legid, x=doy, colour=factor(season)))

ggsave(filename="Figs/Dotplot.jpeg", width=20, height=12)

dat.day <- dat.raw %>% 
  dplyr::filter(doy >= 35, doy <= 220) %>% 
  group_by(legid) %>% 
  summarize(days=n()) %>% 
  ungroup()

#3. Select one year per id----
dat.i <- dat.day %>% 
  dplyr::filter(days==max(days)) %>% 
  left_join(dat.raw) %>% 
  dplyr::select(id, year) %>% 
  unique() %>% 
  group_by(id) %>% 
  sample_n(1) %>% 
  ungroup() %>% 
  left_join(dat.raw) %>% 
  arrange(id, doy) %>% 
  dplyr::select(legid, X, Y)

#4. Split into a list----
dat.s <- group_split(dat.i, legid, .keep=FALSE)
names(dat.s) <- unique(dat.i$legid)

#5. Cluster----
clust.s <- tsclust(dat.s, k = 2L:9L,
                   distance = "dtw_basic", centroid = "pam",
                   seed = 94L)
names(clust.s) <- paste0("k_", 2L:9L)

#6. Apply clusters to data----
dat.clust <- data.frame(legid = unique(dat.i$legid),
                        clust2 = clust.s$k_2@cluster,
                        clust3 = clust.s$k_3@cluster,
                        clust4 = clust.s$k_4@cluster,
                        clust5 = clust.s$k_5@cluster,
                        clust6 = clust.s$k_6@cluster,
                        clust7 = clust.s$k_7@cluster,
                        clust8 = clust.s$k_8@cluster,
                        clust9 = clust.s$k_9@cluster) %>% 
  full_join(dat.i) %>% 
  left_join(dat.raw)

#7. Visualize----
ggplot(dat.clust) +
  geom_line(data=dat.i, aes(x=X, y=Y, group=legid)) +
  geom_point(aes(x=X, y=Y, colour=factor(clust4)))

ggplot(dat.clust) +
#  geom_line(data=dat.i, aes(x=X, y=Y, group=legid)) +
  geom_point(aes(x=X, y=Y, colour=factor(clust4))) + 
  facet_wrap(~season)
