#https://bookdown.org/adam_smith2/bbsbayes_intro_workshop/

library(tidyverse)
library(bbsBayes)
library(shinystan)
library(MuMIn)

source("functions.R")

#1. Load clusters for BBS routes with LBCU on them----
clust <- read.csv("Data/LBCUBBSClusters.csv")

#2. Import bbs monitoring data----
load("~/Library/Application Support/bbsBayes/bbs_raw_data.RData")

#3. Set up loop through # of clusters---
clusters <- c(2:6,8:9)

for(j in 1:length(clusters)){
  
  if(j==1){
    trend.list <- data.frame()
  }
  else
  {
    trend.list <- read.csv("Data/LBCUClusterTrends.csv")
  }

  #4. Subset bbs route clusters----
  clust.j <- clust %>% 
    dplyr::filter(nclust==clusters[j])
  
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
    dplyr::filter(id %in% clust$id) %>% 
    mutate(rt.uni=paste(statenum, Route, sep="-"),
           rt.uni.y=paste(rt.uni, Year, sep="-"))
  
  bbs.j$stratify_by <- "bcr"
  
  #4. Prepare data for bbsbayes----
  dat.j <- prepare_data(strat_data = bbs.j,
                        species_to_run = "Long-billed Curlew",
                        model = "gam",
                        heavy_tailed=TRUE)
  
  #5. Run bbsbayes model----
#  mod.j <- run_model(jags_data = dat.j,
#                        parallel = TRUE,
#                        parameters_to_save = c("n","n3"))
#  mod.j$stratify_by <- "cluster"
  
#  write_rds(mod.j, paste0("bbsBayesModels/LBCU_cluster", clusters[j],"_gam.rds"))
  mod.j <- read_rds(paste0("bbsBayesModels/LBCU_cluster", clusters[j],"_gamye.rds"))
  
  #6. Create annual indices----
  all_area_weights <- utils::read.csv("Data/area_weight.csv") %>% 
    dplyr::filter(nclust==clusters[j]) %>% 
    mutate(area_sq_km = area/1000) %>% 
    rename(region = knncluster) %>% 
    dplyr::select(region, area_sq_km)
  
  indices.j <- generate_indices(jags_mod = mod.j,
                               jags_data = dat.j,
                               regions="stratum")
  
  #7. Calculate trends----
  trends.j <- generate_trends(indices = indices.j,
                            slope=TRUE,
                            Min_year = 1970,
                            Max_year = 2019)
  
  #8. Save output----
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
                dplyr::select(Country, State, route, id, knncluster, knnprob, nclust, X, Y) %>% 
                unique()) %>% 
    left_join(trends.j %>% 
                mutate(knncluster = as.numeric(Region))) %>% 
    rbind(trend.list)
  
  write.csv(trend.list, "Data/LBCUClusterTrends.csv", row.names = FALSE)
  
  print(paste0("Finished cluster ", clusters[j], " of ", length(clusters)))
  
}

#9. Explore drivers of var in trend between clusters-----
trend.var <- trend.list %>% 
  group_by(route, X, Y) %>% 
  summarize(trend.mn = mean(Trend),
            trend.sd = sd(Trend),
            width.mn = mean(Width_of_95_percent_Credible_Interval),
            width.sd = sd(Width_of_95_percent_Credible_Interval),
            abun.mn = mean(Relative_Abundance),
            abun.sd = sd(Relative_Abundance), 
            pres = mean(pres),
            count.mn = mean(count.mn),
            count.sd = mean(count.sd),
            count.max = mean(count.max),
            years = mean(years)) %>% 
  ungroup() %>% 
  dplyr::filter(!is.na(count.sd))

trend.lm <- trend.list %>% 
  group_by(route) %>% 
  lm()

#Plot spatial
mn <- read.csv("Data/LBCUMCLocations.csv")

ggplot(trend.var) +
  geom_path(data=mn, aes(x=X, y=Y, group=id)) +
  geom_point(data=mn, aes(x=X, y=Y), pch=21, size=2, fill="grey50") +
  geom_point(aes(x=X, y=Y, colour=trend.sd), alpha=0.7) +
  scale_colour_viridis_c()

#Correlations
trend.cols <- trend.var %>% 
  dplyr::select(-route, -X, -Y)
cor(trend.cols)
corrplot::corrplot(cor(trend.cols))

ggplot(trend.var) +
  geom_point(aes(x=abun.mn, y=trend.sd))



lm1 <- lm(trend.sd ~ abun.mn, data=trend.var, na.action = "na.fail")
summary(lm1)
