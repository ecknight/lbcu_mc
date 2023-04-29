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

#3. Set up loop for approach----
clusts <- unique(clust$nclust)
indices.list <- list()
trend.list <- list()
for(j in 2:length(clusts)){
  
  #4. Filter data to approach----
  clust.j <- clust %>% 
    dplyr::filter(nclust==clusts[j])
  
  #3. Wrangle bbs data into just routes of interest and add cluster attribution----
  bbs.j <- list()
  
  bbs.j[["route_strat"]] <- bbs_data[["route"]] %>% 
    mutate(id = paste(countrynum, statenum, Route, sep="-")) %>% 
    inner_join(clust.j) %>% 
    mutate(strat_name=knncluster,
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
                     parameters_to_save = c("n","n3"))
  mod.j$stratify_by <- "cluster"
  end.time <- Sys.time()
  
  write_rds(mod.j, paste0("Results/bbsBayesModels/LBCU_cluster_gamye_", clusts[j], ".rds"))
  mod.j <- read_rds(paste0("Results/bbsBayesModels/LBCU_cluster_gamye_", clusts[j], ".rds")) 
  
  #7. Create annual indices----
  all_area_weights <- utils::read.csv("Data/area_weight.csv") %>% 
    dplyr::filter(nclust==clusts[j]) %>% 
    mutate(area_sq_km = area/1000) %>% 
    rename(region = knncluster) %>% 
    dplyr::select(region, area_sq_km)
  
  write.csv(all_area_weights, "/Library/Frameworks/R.framework/Versions/4.1/Resources/library/bbsBayes/composite-regions/cluster.csv", row.names = FALSE)
  
  indices.list[[j]] <- generate_indices(jags_mod = mod.j,
                              jags_data = dat.j,
                              regions="stratum")
  
  #8. Calculate trends----
  trends <- generate_trends(indices = indices.list[[j]],
                            slope=TRUE,
                            Min_year = 2007,
                            Max_year = 2019)
  
  #9. Save output----
  trend.list[[j]] <- data.frame(route = dat.j$route,
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
                dplyr::select(Country, State, route, id, knncluster,  X, Y) %>% 
                unique() %>% 
                mutate(knncluster = as.character(knncluster))) %>% 
    left_join(trends %>% 
                mutate(knncluster = Region)) %>% 
    mutate(nclust = clusts[j])
  
}

trend <- rbindlist(trend.list)
indices <- rbind(data.frame(indices.list[[1]]$data_summary) %>% 
                   mutate(nclust=clusts[[1]]),
                 data.frame(indices.list[[2]]$data_summary) %>% 
                   mutate(nclust=clusts[[2]]))

write.csv(indices, "Data/LBCU_indices_gamye.csv")
write.csv(trend, "Data/LBCU_trend_gamye.csv")

indices <- read.csv("Data/LBCU_indices_gamye.csv")
trend <- read.csv("Data/LBCU_trend_gamye.csv")

trend.sum <- trend %>% 
  dplyr::select(nclust, knncluster, Trend, 'Trend_Q0.025', 'Trend_Q0.975') %>% 
  unique() %>% 
  rename(up = 'Trend_Q0.975', down = 'Trend_Q0.025')

#10. Plot----
bbs.use <- read.csv("Data/BBSRoutesToUse.csv")

plot.map <- ggplot(dplyr::filter(trend, id %in% bbs.use$id)) +
  geom_point(aes(x=X, y=Y, colour=Trend)) +
  scale_colour_gradient2(high="blue", low="red") +
  facet_wrap(~nclust)
plot.map

plot.bar <- ggplot(trend.sum) +
  geom_hline(aes(yintercept=0), linetype="dashed", colour="grey30") +
  geom_point(aes(x=knncluster, y=Trend, colour=Trend), size=2) +
  geom_errorbar(aes(x=knncluster, ymin=down, ymax=up, colour=Trend)) +
  scale_colour_gradient2(high="blue", low="red") +
  facet_wrap(~nclust, scales="free_x")
plot.bar

plot.index <- ggplot(indices) +
  geom_ribbon(aes(x=Year, ymin=Index_q_0.025, ymax=Index_q_0.975, group=factor(Region)), alpha = 0.5) +
  geom_line(aes(x=Year, y=Index, colour=factor(Region))) +
  geom_vline(aes(xintercept=2007), linetype="dashed") +
  facet_wrap(~nclust)
plot.index

ggsave(grid.arrange(plot.map, plot.bar, plot.index, nrow=3), height=18, width=12, units='in', filename="Figs/Trend_eBirdwindow.jpeg")
