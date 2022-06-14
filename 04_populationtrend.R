#https://bookdown.org/adam_smith2/bbsbayes_intro_workshop/

library(tidyverse)
library(bbsBayes)
library(shinystan)
library(gridExtra)

source("functions.R")

#TO DO: CHANGE INTERPRETATION TO TRAJECTORY INSTEAD OF TREND####

#1. Load clusters for BBS routes with LBCU on them----
clust.j <- read.csv("Data/LBCUBBSClusters.csv") 

#2. Import bbs monitoring data----
load("~/Library/Application Support/bbsBayes/bbs_raw_data.RData")

#3. Wrangle bbs data into just routes of interest and add cluster attribution----
bbs.j <- list()

bbs.j[["route_strat"]] <- bbs_data[["route"]] %>% 
  mutate(id = paste(countrynum, statenum, Route, sep="-")) %>% 
  inner_join(clust.j) %>% 
  mutate(strat_name=svmcluster,
         rt.uni=paste(statenum, Route, sep="-"),
         rt.uni.y=paste(rt.uni, Year, sep="-"))

bbs.j[["species_strat"]] <- bbs_data[["species"]]

bbs.j[["bird_strat"]] <- bbs_data[["bird"]] %>% 
  mutate(id = paste(countrynum, statenum, Route, sep="-")) %>% 
  dplyr::filter(id %in% clust.j$id) %>% 
  mutate(rt.uni=paste(statenum, Route, sep="-"),
         rt.uni.y=paste(rt.uni, Year, sep="-"))

bbs.j$stratify_by <- "bcr"  

#4. Prepare data for bbsbayes----
dat.j <- prepare_data(strat_data = bbs.j,
                      species_to_run = "Long-billed Curlew",
                      model = "gamye",
                      heavy_tailed=TRUE)
 
#6. Run bbsbayes model----
start.time <- Sys.time()
mod.j <- run_model(jags_data = dat.j,
                 parallel = TRUE,
                 parameters_to_save = c("n","n3"),
                 # n_burnin = 100,
                 # n_save_steps = 10,
                 # n_iter = 100
                 )
mod.j$stratify_by <- "cluster"
end.time <- Sys.time()

write_rds(mod.j, "bbsBayesModels/LBCU_cluster_gamye.rds")
mod.j <- read_rds("bbsBayesModels/LBCU_cluster_gamye.rds") 

#7. Create annual indices----
all_area_weights <- utils::read.csv("Data/area_weight.csv") %>% 
  mutate(area_sq_km = area/1000) %>% 
  rename(region = svmcluster) %>% 
  dplyr::select(region, area_sq_km)

write.csv(all_area_weights, "/Library/Frameworks/R.framework/Versions/4.1/Resources/library/bbsBayes/composite-regions/cluster.csv", row.names = FALSE)

indices <- generate_indices(jags_mod = mod.j,
                            jags_data = dat.j,
                            regions="stratum")

write.csv(indices$data_summary, "Data/LBCU_indices_gamye.csv")

#8. Calculate trends----
trends <- generate_trends(indices = indices,
                            slope=TRUE,
                            Min_year = 1970,
                            Max_year = 2019)

#9. Save output----
trend.list <- data.frame(route = dat.j$route,
                         count = dat.j$count) %>%  
  mutate(pres = ifelse(count > 0, 1, 0)) %>% 
  group_by(route) %>% 
  summarize(pres = sum(count),
            count.mn = mean(count),
            count.sd = sd(count),
            count.max = max(count),
            years = n()) %>% 
  ungroup() %>% 
  left_join(bbs.j$route_strat %>% 
              rename(route = rt.uni) %>% 
              dplyr::select(Country, State, route, id, svmcluster,  X, Y) %>% 
              unique()) %>% 
  left_join(trends %>% 
              mutate(svmcluster = as.numeric(Region)))

write.csv(trend.list, "Data/LBCU_trend_gamye.csv")

#10. Plot----
plot.map <- ggplot(trend.list) +
  geom_point(aes(x=X, y=Y, colour=Trend)) +
  scale_colour_gradient2(high="blue", low="red")
plot.map

trend <- trend.list %>% 
  dplyr::select(svmcluster, Trend, 'Trend_Q0.025', 'Trend_Q0.975') %>% 
  unique() %>% 
  rename(up = 'Trend_Q0.975', down = 'Trend_Q0.025')

plot.bar <- ggplot(trend) +
  geom_hline(aes(yintercept=0), linetype="dashed", colour="grey30") +
  geom_point(aes(x=svmcluster, y=Trend, colour=Trend), size=2) +
  geom_errorbar(aes(x=svmcluster, ymin=down, ymax=up, colour=Trend)) +
  scale_colour_gradient2(high="blue", low="red")
plot.bar

plot.index <- ggplot(indices$data_summary) +
  geom_ribbon(aes(x=Year, ymin=Index_q_0.025, ymax=Index_q_0.975, group=factor(Region)), alpha = 0.5) +
  geom_line(aes(x=Year, y=Index, colour=factor(Region)))
plot.index

ggsave(grid.arrange(plot.map, plot.bar, plot.index, ncol=3), height=6, width=18, units='in', filename="Figs/Trend.jpeg")
