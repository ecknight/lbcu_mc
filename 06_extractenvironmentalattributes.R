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

#3. Get list of raster KDE files----
files <- data.frame(filepath = list.files("gis", pattern="*.tif", recursive=TRUE, full.names = TRUE)) %>% 
  separate(filepath, into=c("gis", "level", "id"), sep="/", remove=FALSE) %>% 
  mutate(id = str_sub(id, -100, -5)) %>% 
  dplyr::filter(id!="re")

#4. Read in shapefile of 95% isopleths----
kde <- read_sf("gis/shp/kde_individual.shp")

#5. Get Raven density----
#set_ebirdst_access_key("dkao2v8gv378", overwrite=TRUE)
#ebirdst_download("comrav")
ebd <- load_raster("/Users/ellyknight/Library/Application Support/ebirdst/comrav-ERD2019-STATUS-20200930-21f8da4a", product="abundance_seasonal", resolution="hr")

#5. Set up loop----
#val.gee <- data.frame()
task.list <- list()
raven <- data.frame()
for(i in 48:nrow(files)){
  
  #6. Read in raster----
  r <- raster(files$filepath[i]) %>% 
    projectRaster(crs=4326)
  
  #7. Send raster to gee---
  assetId <- sprintf("%s/%s",ee_get_assethome(), files$id[i])
  r.g <- raster_as_ee(r, bucket="lbcu", assetId, overwrite=TRUE)
  
  #8. Get global surface water layer (1984 to 2021)----
  gsw.occ<-ee$Image('JRC/GSW1_3/GlobalSurfaceWater')$select('occurrence')
  gsw.chg<-ee$Image('JRC/GSW1_3/GlobalSurfaceWater')$select('change_norm')
  gsw.mths<-ee$Image('JRC/GSW1_3/GlobalSurfaceWater')$select('seasonality')
  gsw.yrs<-ee$Image('JRC/GSW1_3/GlobalSurfaceWater')$select('recurrence')
  
  #9. Get drought data----
  tclim.drght <- ee$ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')$filter(ee$Filter$date('1984-01-01', '2021-12-31'))$select('pdsi')$mean()
  
  #10. Get landcover data----
  dw.grass <- ee$ImageCollection("GOOGLE/DYNAMICWORLD/V1")$filter(ee$Filter$date('2015-01-01', '2021-12-31'))$select('grass')$mean()
  dw.wetland <- ee$ImageCollection("GOOGLE/DYNAMICWORLD/V1")$filter(ee$Filter$date('2015-01-01', '2021-12-31'))$select('flooded_vegetation')$mean()
  dw.crop <- ee$ImageCollection("GOOGLE/DYNAMICWORLD/V1")$filter(ee$Filter$date('2015-01-01', '2021-12-31'))$select('crops')$mean()
  dw.built <- ee$ImageCollection("GOOGLE/DYNAMICWORLD/V1")$filter(ee$Filter$date('2015-01-01', '2021-12-31'))$select('built')$mean()
  
  #11. Get landcover change data----
  modis1 <- ee$ImageCollection("MODIS/006/MCD12Q1")$filter(ee$Filter$date('2001-01-01', '2001-12-31'))$select('LC_Type1')$mean()
  modis2 <- ee$ImageCollection("MODIS/006/MCD12Q1")$filter(ee$Filter$date('2020-01-01', '2020-12-31'))$select('LC_Type1')$mean()
  
  modis1.grass <- modis1$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0))
  modis2.crop <- modis2$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0))
  modis2.built <- modis2$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0))
  
  modis.crop <- modis1.grass$subtract(modis2.crop)$remap(c(-1,0,1,2), c(0,0,1,0))
  modis.built <- modis1.grass$subtract(modis2.built)$remap(c(-1,0,1,2), c(0,0,1,0))
  
  #12. Multiply rasters together----
  gsw.occ.r <- r.g$multiply(gsw.occ)
  gsw.chg.r <- r.g$multiply(gsw.chg)
  gsw.mths.r <- r.g$multiply(gsw.mths)
  gsw.yrs.r <- r.g$multiply(gsw.yrs)
  tclim.drght.r <- r.g$multiply(tclim.drght)
  dw.grass.r <- r.g$multiply(dw.grass)
  dw.wetland.r <- r.g$multiply(dw.wetland)
  dw.crop.r <- r.g$multiply(dw.crop)
  dw.built.r <- r.g$multiply(dw.built)
  modis.crop.r <- r.g$multiply(modis.crop)
  modis.built.r <- r.g$multiply(modis.built)
  
  #13. Stack them----
  stack <- gsw.occ.r$addBands(gsw.chg.r)$addBands(gsw.mths.r)$addBands(gsw.yrs.r)$addBands(tclim.drght.r)$addBands(dw.grass.r)$addBands(dw.wetland.r)$addBands(dw.crop.r)$addBands(dw.built.r)$addBands(modis.crop.r)$addBands(modis.built.r)
  
  #14. Send 95% isopleth to gee----
  kde.i <- kde %>% 
    dplyr::filter(id==str_sub(files$id[i], 5, 100))
  kde.g <- sf_as_ee(kde.i)
  
  #15. Take mean of raster pixels within 95% isopleth---
  stack.mn <- stack$reduceRegions(reducer=ee$Reducer$mean(), 
                                  collection=kde.g,
                                  scale=30)
  
  #16. Export the values----
  task_vector <- ee_table_to_gcs(collection=stack.mn,
                                 bucket="lbcu", 
                                 fileFormat = "CSV",
                                 fileNamePrefix = files$id[i])
  task_vector$start()
  task.list[[i]] <- task_vector
#  ee_monitoring(task_vector, max_attempts=1000)
#  ee_gcs_to_local(task = task_vector, dsn=paste0("gis/output/", files$id[i], ".csv"))
  
  #17. Get Raven density----
  kde.ebd <- kde.i %>% 
    st_transform(crs(ebd))
  
  r.ebd <- r %>% 
    projectRaster(crs=crs(ebd))
  
  ebd.i <- ebd %>% 
    crop(r.ebd) %>% 
    resample(r.ebd) %>% 
    mask(r.ebd)
  
  ebd.m <- ebd.i*r.ebd
  raven = data.frame(raven = mean(ebd.m@data@values, na.rm=TRUE),
                     id=str_sub(files$id[i], 5, 100)) %>% 
    rbind(raven)
  
  print(paste0("Finished individual ", i, " of ", nrow(files)))
}

#17. Put it together----
val.gee <- data.frame()
for(i in 1:length(task.list)){
  
  try(ee_gcs_to_local(task = task.list[[i]], dsn=paste0("gis/output/", files$id[i], ".csv")))
  
  val.i <- try(read.csv(paste0("gis/output/", files$id[i], ".csv")))
  
  if(class(val.i)=="FILL IN"){
    val.gee <- val.i %>% 
      dplyr::select(-.geo) %>% 
      rename(water.occur = b1,
             water.change = b1_1,
             water.months = b1_2,
             water.years = b1_3,
             drought = b1_4,
             grass = b1_5,
             wetland = b1_6,
             crop = b1_7,
             built = b1_8,
             cropchg = b1_9,
             builtchg = b1_10) %>% 
      rbind(val.gee)
  }

}

covs <- val.gee %>% 
  left_join(raven)

write.csv(val.gee, "Data/LBCU_environvars.csv", row.names = FALSE)