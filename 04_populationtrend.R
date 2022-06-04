#https://bookdown.org/adam_smith2/bbsbayes_intro_workshop/

library(tidyverse)
library(bbsBayes)
library(shinystan)
library(gridExtra)

source("functions.R")

#1. Load clusters for BBS routes with LBCU on them----
clust <- read.csv("Data/LBCUBBSClusters.csv") 

#2. Import bbs monitoring data----
load("~/Library/Application Support/bbsBayes/bbs_raw_data.RData")

#3. Wrangle bbs data into just routes of interest and add cluster attribution----
bbs <- list()

bbs[["route_strat"]] <- bbs_data[["route"]] %>% 
  mutate(id = paste(countrynum, statenum, Route, sep="-")) %>% 
  inner_join(clust) %>% 
  mutate(strat_name=knncluster,
         rt.uni=paste(statenum, Route, sep="-"),
         rt.uni.y=paste(rt.uni, Year, sep="-"))

bbs[["species_strat"]] <- bbs_data[["species"]]

bbs[["bird_strat"]] <- bbs_data[["bird"]] %>% 
  mutate(id = paste(countrynum, statenum, Route, sep="-")) %>% 
  dplyr::filter(id %in% clust$id) %>% 
  mutate(rt.uni=paste(statenum, Route, sep="-"),
         rt.uni.y=paste(rt.uni, Year, sep="-"))

bbs$stratify_by <- "bcr"  

#4. Prepare data for bbsbayes----
dat <- prepare_data(strat_data = bbs,
                      species_to_run = "Long-billed Curlew",
                      model = "gam",
                      heavy_tailed=TRUE)

#6. Run bbsbayes model----
#10:58
mod <- run_model(jags_data = dat,
                       parallel = TRUE,
                       parameters_to_save = c("n","n3"))
mod$stratify_by <- "cluster"
  
write_rds(mod, "bbsBayesModels/LBCU_cluster_gam.rds")
mod <- read_rds("bbsBayesModels/LBCU_cluster_gam.rds") 

#7. Create annual indices----
all_area_weights <- utils::read.csv("Data/area_weight.csv") %>% 
  mutate(area_sq_km = area/1000) %>% 
  rename(region = knncluster) %>% 
  dplyr::select(region, area_sq_km)
write.csv(all_area_weights, "/Library/Frameworks/R.framework/Versions/4.1/Resources/library/bbsBayes/composite-regions/cluster.csv", row.names = FALSE)

indices.j <- generate_indices(jags_mod = mod,
                              jags_data = dat,
                              regions="stratum")

tp <- plot_indices(indices.j)

pdf(file = "Figs/Trajectories.pdf")
print(tp)
dev.off()

#8. Calculate trends----
trends.j <- generate_trends(indices = indices.j,
                            slope=TRUE,
                            Min_year = 1970,
                            Max_year = 2019)

#9. Save output----
trend.list <- data.frame(route = dat$route,
                         count = dat$count) %>%  
  mutate(pres = ifelse(count > 0, 1, 0)) %>% 
  group_by(route) %>% 
  summarize(pres = sum(count),
            count.mn = mean(count),
            count.sd = sd(count),
            count.max = max(count),
            years = n()) %>% 
  ungroup() %>% 
  left_join(bbs$route_strat %>% 
              rename(route = rt.uni) %>% 
              dplyr::select(Country, State, route, id, knncluster, knnprob, nclust, X, Y) %>% 
              unique()) %>% 
  left_join(trends.j %>% 
              mutate(knncluster = as.numeric(Region)))

write.csv(trend.list, "Data/LBCU_cluster_gam.csv")

#10. Plot----
plot.map <- ggplot(trend.list) +
  geom_point(aes(x=X, y=Y, colour=Trend)) +
  scale_colour_gradient2(high="blue", low="red")
plot.map

trend <- trend.list %>% 
  dplyr::select(knncluster, Trend, 'Trend_Q0.025', 'Trend_Q0.975') %>% 
  unique() %>% 
  rename(up = 'Trend_Q0.975', down = 'Trend_Q0.025')

plot.bar <- ggplot(trend) +
  geom_hline(aes(yintercept=0), linetype="dashed", colour="grey30") +
  geom_point(aes(x=knncluster, y=Trend, colour=Trend), size=2) +
  geom_errorbar(aes(x=knncluster, ymin=down, ymax=up, colour=Trend)) +
  scale_colour_gradient2(high="blue", low="red")
plot.bar

ggsave(grid.arrange(plot.map, plot.bar), height=8, width=6, units='in', filename="Figs/Trend.jpeg")
