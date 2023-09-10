#https://bookdown.org/adam_smith2/bbsbayes_intro_workshop/

library(tidyverse)
library(bbsBayes)
library(shinystan)
library(gridExtra)
library(data.table)

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
                      model = "gamye",
                      heavy_tailed=TRUE)

#6. Run bbsbayes model----
start.time <- Sys.time()
mod <- run_model(jags_data = dat,
                   parallel = TRUE,
                   parameters_to_save = c("n","n3"))
mod$stratify_by <- "cluster"
end.time <- Sys.time()

write_rds(mod, "Results/bbsBayesModels/LBCU_cluster_gamye_3.rds")
mod <- read_rds("Results/bbsBayesModels/LBCU_cluster_gamye_3.rds")

#7. Create annual indices----
all_area_weights <- utils::read.csv("Data/area_weight.csv") %>% 
  mutate(area_sq_km = area/1000) %>% 
  rename(region = knncluster) %>% 
  dplyr::select(region, area_sq_km)

write.csv(all_area_weights, "/Library/Frameworks/R.framework/Versions/4.1/Resources/library/bbsBayes/composite-regions/cluster.csv", row.names = FALSE)

indices <- generate_indices(jags_mod = mod,
                            jags_data = dat,
                            regions="stratum")

#8. Calculate trends----
trends <- generate_trends(indices = indices,
                          slope=TRUE,
                          Min_year = 2007,
                          Max_year = 2019)

#9. Save output----
trend.out <- data.frame(route = dat$route,
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
              dplyr::select(Country, State, route, id, knncluster,  X, Y) %>% 
              unique() %>% 
              mutate(knncluster = as.character(knncluster))) %>% 
  left_join(trends %>% 
              mutate(knncluster = Region))


indices.out <- data.frame(indices$data_summary)

write.csv(indices.out, "Data/LBCU_indices_gamye.csv")
write.csv(trend.out, "Data/LBCU_trend_gamye.csv")

indices.out <- read.csv("Data/LBCU_indices_gamye.csv")
trend.out <- read.csv("Data/LBCU_trend_gamye.csv")

trend.sum <- trend.out %>% 
  dplyr::select(knncluster, Trend, 'Trend_Q0.025', 'Trend_Q0.975') %>% 
  unique() %>% 
  rename(up = 'Trend_Q0.975', down = 'Trend_Q0.025')

#10. Plot----
bbs.use <- read.csv("Data/BBSRoutesToUse.csv")

plot.map <- ggplot(dplyr::filter(trend.out, id %in% bbs.use$id)) +
  geom_point(aes(x=X, y=Y, colour=Trend)) +
  scale_colour_gradient2(high="blue", low="red")
plot.map

plot.bar <- ggplot(trend.sum) +
  geom_hline(aes(yintercept=0), linetype="dashed", colour="grey30") +
  geom_point(aes(x=knncluster, y=Trend, colour=Trend), size=2) +
  geom_errorbar(aes(x=knncluster, ymin=down, ymax=up, colour=Trend)) +
  scale_colour_gradient2(high="blue", low="red")
plot.bar

plot.index <- ggplot(indices.out) +
  geom_ribbon(aes(x=Year, ymin=Index_q_0.025, ymax=Index_q_0.975, group=factor(Region)), alpha = 0.5) +
  geom_line(aes(x=Year, y=Index, colour=factor(Region))) +
  geom_vline(aes(xintercept=2007), linetype="dashed")
plot.index

ggsave(grid.arrange(plot.map, plot.bar, plot.index, nrow=3), height=18, width=12, units='in', filename="Figs/Trend_eBirdwindow.jpeg")
