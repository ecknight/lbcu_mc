#https://bookdown.org/adam_smith2/bbsbayes_intro_workshop/

library(tidyverse)
library(bbsBayes)
library(shinystan)

#1. Load clusters for BBS routes with LBCU on them----
clust <- read.csv("LBCUBBSClusters.csv")

#2. Import bbs monitoring data----
load("~/Library/Application Support/bbsBayes/bbs_raw_data.RData")

#3. Set up loop through # of clusters---
clusters <- c(2:10)

for(j in 1:length(clusters)){
  
  #4. Subset bbs route clusters----
  clust.j <- clust %>% 
    dplyr::filter(nclust==clusters[[j]])
  
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
                        model = "gamye",
                        heavy_tailed=TRUE)
  
  #5. Run bbsbayes model----
  mod.j <- run_model(jags_data = dat.j,
                        n_iter = 24000, #higher than the default 10,000
                        n_burnin = 20000,
                        n_chains = 3,
                        n_thin = 20, #saves memory by retaining only 1/20 posterior samples
                        parallel = TRUE,
                        inits = NULL,
                        parameters_to_save = c("n","n3"))
  
  #6. Create annual indices----
  indices.j <- generate_indices(jags_mod = mod.j,
                               jags_data = dat.j,
                               regions="stratum")
  
  #7. Calculate trends----
  trends <- generate_trends(indices = indices.j,
                            slope=TRUE,
                            Min_year = 1970,
                            Max_year = 2019)
  
  #8. Save output----
  
  
  print(paste0("Finished cluster ", clusters[j], " of ", length(clusters)))
  
  
  
}


  