library(sf)
library(raster)
library(tidyverse)
library(lubridate)
library(rgee)
library(data.table)

options(scipen=99999)

#PREAMBLE####

#1. Initialize rgee----
ee_Initialize(gcs=TRUE)
ee_check()

#2. Create GCS bucket----
# Create your own container
# project_id <- ee_get_earthengine_path() %>% 
#   list.files(., "\\.json$", full.names = TRUE) %>% 
#   jsonlite::read_json() %>% 
#   '$'(project_id) # Get the Project ID
# 
# googleCloudStorageR::gcs_create_bucket("lbcu", projectId = project_id)

#A. INDIVIDUAL ATTRIBUTES####

#1. Read in shapefile of HRs and get centroid----
#just use one nclust group, no need to do twice
kde.raw <- read_sf("gis/shp/kde_individual.shp") %>% 
  st_make_valid() %>% 
  st_centroid() %>% 
  mutate(year = as.numeric(year)-1) %>% 
  dplyr::filter(nclust=="3")

kde.springt <- kde.raw %>% 
  dplyr::filter(season=="springmig") %>% 
  mutate(year = year+1,
         season = "springmigt")

kde <- rbind(kde.raw, kde.springt)

#2. Set up to loop through years----
#Spring vars - t, t+1
#Breed vars - t+1
#Fall vars - t+1
#Winter vars - t+1
years <- sort(unique(kde$year))
for(i in 1:length(years)){
  
  kde.i <- kde %>% 
    dplyr::filter(year==years[i])
  
  #3. Get drought data----
  tclim.drght <- ee$ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')$filter(ee$Filter$date(paste0(years[i], "-01-01"), paste0(years[i], "-12-31")))$select('pdsi')$mean()
  
  #4. Get MODIS landcover----
  modis <- ee$ImageCollection("MODIS/006/MCD12Q1")$filter(ee$Filter$date(paste0(years[i], "-01-01"), paste0(years[i], "-12-31")))$select('LC_Type1')$mean()
  
  modis.grass <- modis$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0))
  modis.crop <- modis$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0))
  modis.built <- modis$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0))
  modis.wetland <- ee$ImageCollection("MODIS/006/MCD12Q1")$filter(ee$Filter$date(paste0(years[i], "-01-01"), paste0(years[i], "-12-31")))$select('LC_Prop3')$mean()$remap(c(1,2,3,10,20,27,30,40,50), c(0,0,0,0,0,0,0,1,0))
  
  #5. Get global surface water layer (1984 to 2021)----
  gsw.chg<-ee$Image('JRC/GSW1_3/GlobalSurfaceWater')$select('change_norm')
  gsw.mths<-ee$Image('JRC/GSW1_3/GlobalSurfaceWater')$select('seasonality')
  
  #6. Stack them----
  stack <- gsw.chg$addBands(gsw.mths)$addBands(tclim.drght)$addBands(modis.crop)$addBands(modis.built)$addBands(modis.grass)$addBands(modis.wetland)$rename("change", "seasonality", "drought", "crop", "built", "grass", "wetland")
  
  #7. Set up loop for individuals----
  for(j in 1:nrow(kde.i)){
    
    #8. Get centroid and send to gee----
    kde.j <- kde.i[j,]
    data <- sf_as_ee(kde.j)
    
    #9. Buffer----
    buffer_points <- function(feature){
      properties <- c("Date", "ID", "Radius", "Type", "date_millis", "ptID", "ptIDn", "X", "Y", "timestamp", "uniq")
      new_geometry <- feature$geometry()$buffer(kde.j$rad*1000)
      ee$Feature(new_geometry)$copyProperties(feature, properties)
    }
    
    data.buff <- data$map(buffer_points)
    
    #10. Extract mean values---
    stack.mn <- stack$reduceRegions(reducer=ee$Reducer$mean(),
                                    collection=data.buff,
                                    scale=30)
    
    #11. Export the values----
    task_vector <- ee_table_to_gcs(collection=stack.mn,
                                   bucket="lbcu",
                                   fileFormat = "CSV",
                                   fileNamePrefix = kde.j$id)
    task_vector$start()
    ee_monitoring(task_vector, max_attempts=1000)
    
    #12. Download----
    ee_gcs_to_local(task = task_vector, dsn=paste0("gis/output_ind/", kde.j$id, ".csv"))
    
    print(paste0("Finished individual ", j, " of ", nrow(kde.i), " for year ", years[i]))
  }
  
}

#13. Read in output----
files <- list.files("gis/output_ind")

dat <- data.frame()
for(i in 1:length(files)){
  dat <- read.csv(paste0("gis/output/", files[i])) %>%
    dplyr::select(-.geo, -system.index) %>% 
    mutate(id = str_sub(files[i], -100, -5)) %>% 
    rbind(dat)
}

#14. Join to other data and tidy----
covs <- dat  %>% 
  separate(id, into=c("nclust", "group", "season", "bird", "year", "cluster"), remove=FALSE) %>% 
  mutate(year = as.numeric(year)) %>% 
#  mutate(id=paste(group, season, birdid, year, cluster, sep="-")) %>% 
  left_join(kde.raw %>% 
              st_coordinates() %>% 
              data.frame() %>% 
              cbind(data.frame(kde.raw)) %>% 
              dplyr::select(-id, -geometry)) %>% 
  mutate(change = ifelse(is.na(change), 0, change),
         seasonality = ifelse(is.na(seasonality), 0, seasonality),
         crop = ifelse(is.na(crop), 0, crop),
         drought = ifelse(is.na(drought), 0, drought),
         grass = ifelse(is.na(grass), 0, grass),
         wetland = ifelse(is.na(wetland), 0, wetland))

#15. Add manual cluster back in----
n <- c("3", "manual")

clust <- read.csv("Data/LBCUKDEClusters.csv") %>% 
  dplyr::filter(nclust %in% n,
                !is.na(X)) %>% 
  group_by(id, nclust) %>% 
  summarize(group = round(mean(group))) %>% 
  ungroup() %>% 
  rename(bird = id) %>% 
  mutate(bird = as.character(bird))

covs.clust <- clust %>% 
  dplyr::filter(bird %in% covs$bird) %>% 
  full_join(covs %>% 
              dplyr::select(-nclust, -group))

#16. Save-----
write.csv(covs.clust, "Data/LBCU_environvars.csv", row.names = FALSE)


#B. REGION ATTRIBUTES####

#1. Read in shapefile of HRs and get centroid----
kde <- read_sf("gis/shp/kde_region.shp") %>% 
  st_make_valid() 

#2. Send to gee---
data <- sf_as_ee(kde)

#2. Set up to loop through years----
years <- c(2006:2020)
for(i in 1:length(years)){
  
  #3. Get drought data----
  tclim.drght <- ee$ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')$filter(ee$Filter$date(paste0(years[i], "-01-01"), paste0(years[i], "-12-31")))$select('pdsi')$mean()
  
  #4. Get MODIS landcover----
  modis <- ee$ImageCollection("MODIS/006/MCD12Q1")$filter(ee$Filter$date(paste0(years[i], "-01-01"), paste0(years[i], "-12-31")))$select('LC_Type1')$mean()
  
  modis.grass <- modis$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0))
  modis.crop <- modis$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0))
  modis.built <- modis$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0))
  modis.wetland <- ee$ImageCollection("MODIS/006/MCD12Q1")$filter(ee$Filter$date(paste0(years[i], "-01-01"), paste0(years[i], "-12-31")))$select('LC_Prop3')$mean()$remap(c(1,2,3,10,20,27,30,40,50), c(0,0,0,0,0,0,0,1,0))
  
  #5. Get global surface water layer (1984 to 2021)----
  gsw.chg<-ee$Image('JRC/GSW1_3/GlobalSurfaceWater')$select('change_norm')
  gsw.mths<-ee$Image('JRC/GSW1_3/GlobalSurfaceWater')$select('seasonality')
  
  #6. Stack them----
  stack <- gsw.chg$addBands(gsw.mths)$addBands(tclim.drght)$addBands(modis.crop)$addBands(modis.built)$addBands(modis.grass)$addBands(modis.wetland)$rename("change", "seasonality", "drought", "crop", "built", "grass", "wetland")
  
  #10. Extract mean values---
  stack.mn <- stack$reduceRegions(reducer=ee$Reducer$mean(),
                                  collection=data,
                                  scale=1000)
  
  #11. Export the values----
  task_vector <- ee_table_to_gcs(collection=stack.mn,
                                 bucket="lbcu",
                                 fileFormat = "CSV",
                                 fileNamePrefix = years[i])
  task_vector$start()
  ee_monitoring(task_vector, max_attempts=1000)
  
  #12. Download----
  ee_gcs_to_local(task = task_vector, dsn=paste0("gis/output_reg/", years[i], ".csv"))
  
  print(paste0("Finished year ", i, " of ", length(years)))
  
}

#13. Read in output----
files <- list.files("gis/output_reg")

dat <- data.frame()
for(i in 1:length(files)){
  dat <- read.csv(paste0("gis/output_reg/", files[i])) %>%
    dplyr::select(-.geo, -system.index) %>% 
    mutate(year = str_sub(files[i], -100, -5)) %>% 
    rbind(dat)
}

#14. Save-----
write.csv(dat, "Data/LBCU_environvars_region.csv", row.names = FALSE)
