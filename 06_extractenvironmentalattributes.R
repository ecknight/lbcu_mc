library(sf)
library(raster)
library(tidyverse)
library(lubridate)
library(rgee)
library(data.table)
library(ebirdst)

options(scipen=99999)

#TO DO: CONSIDER GEOPOTENTIAL HEIGHT####
#TO DO: CONSIDER CENTER PIVOT####
#TO DO: FIX RAVEN####

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
kde <- read_sf("gis/shp/kde_individual.shp") %>% 
  mutate(id = paste(kdecluster, season, id, year, cluster, sep="-")) %>% 
  st_make_valid() %>% 
  st_centroid()

#5. Get Raven density----
#set_ebirdst_access_key("dkao2v8gv378", overwrite=TRUE)
#ebirdst_download("comrav")
ebd <- load_raster("/Users/ellyknight/Library/Application Support/ebirdst/comrav-ERD2019-STATUS-20200930-21f8da4a", product="abundance_seasonal", resolution="hr")

#6. Get global surface water layer (1984 to 2021)----
gsw.occ<-ee$Image('JRC/GSW1_3/GlobalSurfaceWater')$select('occurrence')
gsw.chg<-ee$Image('JRC/GSW1_3/GlobalSurfaceWater')$select('change_norm')
gsw.mths<-ee$Image('JRC/GSW1_3/GlobalSurfaceWater')$select('seasonality')
gsw.yrs<-ee$Image('JRC/GSW1_3/GlobalSurfaceWater')$select('recurrence')

#7. Get drought data----
tclim.drght <- ee$ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')$filter(ee$Filter$date('1984-01-01', '2021-12-31'))$select('pdsi')$mean()

#8. Get landcover data----
dw.grass <- ee$ImageCollection("GOOGLE/DYNAMICWORLD/V1")$filter(ee$Filter$date('2015-01-01', '2021-12-31'))$select('grass')$mean()
dw.wetland <- ee$ImageCollection("GOOGLE/DYNAMICWORLD/V1")$filter(ee$Filter$date('2015-01-01', '2021-12-31'))$select('flooded_vegetation')$mean()
dw.crop <- ee$ImageCollection("GOOGLE/DYNAMICWORLD/V1")$filter(ee$Filter$date('2015-01-01', '2021-12-31'))$select('crops')$mean()
dw.built <- ee$ImageCollection("GOOGLE/DYNAMICWORLD/V1")$filter(ee$Filter$date('2015-01-01', '2021-12-31'))$select('built')$mean()

#9. Get landcover change data----
modis1 <- ee$ImageCollection("MODIS/006/MCD12Q1")$filter(ee$Filter$date('2001-01-01', '2001-12-31'))$select('LC_Type1')$mean()
modis2 <- ee$ImageCollection("MODIS/006/MCD12Q1")$filter(ee$Filter$date('2020-01-01', '2020-12-31'))$select('LC_Type1')$mean()

modis1.grass <- modis1$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0))
modis2.crop <- modis2$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0))
modis2.built <- modis2$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0))

modis.crop <- modis1.grass$subtract(modis2.crop)$remap(c(-1,0,1,2), c(0,0,1,0))
modis.built <- modis1.grass$subtract(modis2.built)$remap(c(-1,0,1,2), c(0,0,1,0))

#13. Stack them----
stack <- gsw.occ$addBands(gsw.chg)$addBands(gsw.mths)$addBands(gsw.yrs)$addBands(tclim.drght)$addBands(dw.grass)$addBands(dw.wetland)$addBands(dw.crop)$addBands(dw.built)$addBands(modis.crop)$addBands(modis.built)

#10. Set up loop----
task.list <- list()
raven <- data.frame()
#i <- 1
for(i in 1:nrow(kde)){
  
  #11. Get centroid and send to gee----
  kde.i <- kde[i,]
  data <- sf_as_ee(kde.i)
  
  #12. Buffer----
  buffer_points <- function(feature){
    properties <- c("Date", "ID", "Radius", "Type", "date_millis", "ptID", "ptIDn", "X", "Y", "timestamp", "uniq")
    new_geometry <- feature$geometry()$buffer(kde.i$rad*1000)
    ee$Feature(new_geometry)$copyProperties(feature, properties)
  }
  
  data.buff <- data$map(buffer_points)
  
  #15. Extract mean values---
  stack.mn <- stack$reduceRegions(reducer=ee$Reducer$mean(), 
                                  collection=data.buff,
                                  scale=30)
  
  #16. Export the values----
  task_vector <- ee_table_to_gcs(collection=stack.mn,
                                 bucket="lbcu", 
                                 fileFormat = "CSV",
                                 fileNamePrefix = kde.i$id[i])
  task_vector$start()
  task.list[[i]] <- task_vector
  ee_monitoring(task.list[[i]], max_attempts=1000)
#  ee_gcs_to_local(task = task_vector, dsn=paste0("gis/output/", kde.i$id[i], ".csv"))
  
  #17. Get Raven density----
  kde.ebd <- kde.i %>% 
    st_transform(crs(ebd)) %>% 
#    st_buffer(kde.i$rad*1000) %>% 
    as_Spatial()

  raven <- data.frame(raven = raster::extract(ebd, kde.ebd, buffer=kde.i$rad*1000, fun=mean),
                       id=kde$id[i]) %>% 
  rbind(raven)
  
  print(paste0("Finished individual ", i, " of ", nrow(kde)))
}

save.image("gis/TaskList.R")

#17. Download results----
for(i in 1:nrow(kde)){
  
  try(ee_gcs_to_local(task = task.list[[i]], dsn=paste0("gis/output/", kde$id[i], ".csv")))

}

#18. Read in results----
files <- list.files("gis/output")

dat <- data.frame()
for(i in 1:length(files)){
  dat <- read.csv(paste0("gis/output/", files[i])) %>%
    dplyr::select(-.geo, -system.index) %>% 
    mutate(id = str_sub(files[i], -100, -5)) %>% 
    rbind(dat)
}

#19. Join to other data and tidy----
covs <- dat %>% 
  left_join(raven) %>% 
  left_join(kde %>% 
              st_coordinates() %>% 
              data.frame() %>% 
              cbind(data.frame(kde)) %>% 
              dplyr::select(id, area, rad)) %>% 
  mutate(change_norm = ifelse(is.na(change_norm), 0, change_norm),
         occurrence = ifelse(is.na(occurrence), 0, occurrence),
         recurrence = ifelse(is.na(recurrence), 0, recurrence),
         seasonality = ifelse(is.na(seasonality), 0, seasonality),
         raven = ifelse(is.na(raven), 0, raven)) %>% 
  rename(covcrop = remapped,
         covbuilt = remapped_1) %>% 
  separate(id, into=c("kdecluster", "season", "id", "year", "cluster"))

#20. Save-----
write.csv(covs, "Data/LBCU_environvars.csv", row.names = FALSE)
