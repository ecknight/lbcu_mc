#https://bookdown.org/adam_smith2/bbsbayes_intro_workshop/

library(bbsBayes)
library(tidyverse)
library(shinystan)

#get data
fetch_bbs_data(force=TRUE)

#stratify
strat_data <- stratify(by = "latlong")

#Prepare data
jags_data <- prepare_data(strat_data = strat_data,
                          species_to_run = "Long-billed Curlew",
                          model = "gamye",
                          heavy_tailed=TRUE)

#Run model
jags_mod <- run_model(jags_data = jags_data,
                      n_iter = 24000, #higher than the default 10,000
                      n_burnin = 20000,
                      n_chains = 3,
                      n_thin = 20, #saves memory by retaining only 1/20 posterior samples
                      parallel = TRUE,
                      inits = NULL,
                      parameters_to_save = c("n","n3"))

#diagnostics
jags_mod$n.eff
jags_mod$Rhat

#my_sso <- shinystan::launch_shinystan(shinystan::as.shinystan(jags_mod$samples, model_name = "my_toy_example"))

#save
write_rds(jags_mod, "bbsBayesModels/LBCU_latlong_gamye.rds")
jags_mod <- read_rds("bbsBayesModels/LBCU_latlong_gamye.rds")

#create annual indices
indices <- generate_indices(jags_mod = jags_mod,
                            jags_data = jags_data)

write.csv(indices$data_summary, "bbsBayesModels/LBCU_trajectories_latlong_gamye.csv")

#Calculate trends
trends <- generate_trends(indices = indices,
                          slope=TRUE,
                          Min_year = 1970,
                          Max_year = 2019)

write.csv(trends, "bbsBayesModels/Trends_gamye.csv")

#Map
mp <- generate_map(trends, 
                   select=TRUE,
                   stratify_by="latlong",
                   species = "Long-billed Curlew") +
  
print(mp)

#Map with clusters
bw.sf <- st_as_sf(bw, coords=c("X_breed", "Y_breed"), crs=3857) %>% 
  st_transform(crs=st_crs(map))  %>% 
  st_coordinates() %>% 
  data.frame() %>% 
  rename(lat=Y, lon=X) %>% 
  cbind(bw)

mp +
  geom_point(data=bw.sf, aes(x=lon, y=lat, colour=factor(bwcluster8))) +
  scale_colour_viridis_d()