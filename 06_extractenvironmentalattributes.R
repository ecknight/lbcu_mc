library(sf)
library(raster)
library(tidyverse)
library(lubridate)
library(rgee)
library(data.table)
library(ebirdst)

options(scipen=99999)

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

#4. Read in shapefile of HRs and get centroid----
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

#5. Set up to loop through years----
#Spring vars - t, t+1
#Breed vars - t+1
#Fall vars - t+1
#Winter vars - t+1
years <- sort(unique(kde$year))
for(i in 11:length(years)){
  
  kde.i <- kde %>% 
    dplyr::filter(year==years[i])
  
  #7. Get drought data----
  tclim.drght <- ee$ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')$filter(ee$Filter$date(paste0(years[i], "-01-01"), paste0(years[i], "-12-31")))$select('pdsi')$mean()
  
  #8. Get MODIS landcover----
  modis <- ee$ImageCollection("MODIS/006/MCD12Q1")$filter(ee$Filter$date(paste0(years[i], "-01-01"), paste0(years[i], "-12-31")))$select('LC_Type1')$mean()
  
  modis.grass <- modis$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0))
  modis.crop <- modis$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0))
  modis.built <- modis$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0))
  modis.wetland <- ee$ImageCollection("MODIS/006/MCD12Q1")$filter(ee$Filter$date(paste0(years[i], "-01-01"), paste0(years[i], "-12-31")))$select('LC_Prop3')$mean()$remap(c(1,2,3,10,20,27,30,40,50), c(0,0,0,0,0,0,0,1,0))
  
  #9. Get global surface water layer (1984 to 2021)----
  gsw.chg<-ee$Image('JRC/GSW1_3/GlobalSurfaceWater')$select('change_norm')
  gsw.mths<-ee$Image('JRC/GSW1_3/GlobalSurfaceWater')$select('seasonality')
  
  #10. Stack them----
  stack <- gsw.chg$addBands(gsw.mths)$addBands(tclim.drght)$addBands(modis.crop)$addBands(modis.built)$addBands(modis.grass)$addBands(modis.wetland)$rename("change", "seasonality", "drought", "crop", "built", "grass", "wetland")
  
  
  
  #11. Set up loop for individuals----
  for(j in 1:nrow(kde.i)){
    
    #12. Get centroid and send to gee----
    kde.j <- kde.i[j,]
    data <- sf_as_ee(kde.j)
    
    #13. Buffer----
    buffer_points <- function(feature){
      properties <- c("Date", "ID", "Radius", "Type", "date_millis", "ptID", "ptIDn", "X", "Y", "timestamp", "uniq")
      new_geometry <- feature$geometry()$buffer(kde.j$rad*1000)
      ee$Feature(new_geometry)$copyProperties(feature, properties)
    }
    
    data.buff <- data$map(buffer_points)
    
    #14. Extract mean values---
    stack.mn <- stack$reduceRegions(reducer=ee$Reducer$mean(),
                                    collection=data.buff,
                                    scale=30)
    
    #15. Export the values----
    task_vector <- ee_table_to_gcs(collection=stack.mn,
                                   bucket="lbcu",
                                   fileFormat = "CSV",
                                   fileNamePrefix = kde.j$id)
    task_vector$start()
    ee_monitoring(task_vector, max_attempts=1000)
    
    #16. Download----
    ee_gcs_to_local(task = task_vector, dsn=paste0("gis/output/", kde.j$id, ".csv"))
    
    print(paste0("Finished individual ", j, " of ", nrow(kde.i), " for year ", years[i]))
  }
  
}

#17. Read in output----
files <- list.files("gis/output")

dat <- data.frame()
for(i in 1:length(files)){
  dat <- read.csv(paste0("gis/output/", files[i])) %>%
    dplyr::select(-.geo, -system.index) %>% 
    mutate(id = str_sub(files[i], -100, -5)) %>% 
    rbind(dat)
}

#18. Join to other data and tidy----
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

#19. Save-----
write.csv(covs, "Data/LBCU_environvars.csv", row.names = FALSE)
