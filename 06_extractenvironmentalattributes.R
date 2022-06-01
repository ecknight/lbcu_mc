library(sf)
library(raster)
library(tidyverse)
library(lubridate)
library(rgee)
library(data.table)

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

#3. Get list of raster KDE files----
files <- data.frame(filepath = list.files("gis", pattern="*.tif", recursive=TRUE, full.names = TRUE)) %>% 
  separate(filepath, into=c("gis", "level", "id"), sep="/", remove=FALSE) %>% 
  mutate(id = str_sub(id, -100, -5))

#4. Read in shapefile of 95% isopleths----
kde <- read_sf("gis/kde.shp")

#4. Read in raster----
i <- 578
r <- raster(files$filepath[i]) %>% 
  projectRaster(crs=4326)

#5. Send raster to gee---
assetId <- sprintf("%s/%s",ee_get_assethome(), files$id[i])
r.g <- raster_as_ee(r, bucket="lbcu", assetId, overwrite=TRUE)

#6. Get global surface water layer----
gsw.occ<-ee$Image('JRC/GSW1_3/GlobalSurfaceWater')$select('occurrence')
gsw.chg<-ee$Image('JRC/GSW1_3/GlobalSurfaceWater')$select('change_norm')

#7. Multiply rasters together----
gsw.occ.r <- r.g$multiply(gsw.occ)
gsw.chg.r <- r.g$multiply(gsw.chg)

#8. Stack them----

#9. Send 95% isopleth to gee----
kde.i <- kde %>% 
  dplyr::filter(id==str_sub(files$id[i], 5, 100))
kde.g <- sf_as_ee(kde.i)

#10. Take mean of raster pixels within 95% isopleth---
gsw.occ.mn <- gsw.occ.r$reduceRegions(reducer=ee$Reducer$mean(), 
                                      collection=kde.g,
                                      scale=30)
gsw.chg.mn <- gsw.chg.r$reduceRegions(reducer=ee$Reducer$mean(), 
                                      collection=kde.g,
                                      scale=30)

#11. Export the values----
task_vector <- ee_table_to_gcs(collection=gsw.occ.mn,
                               bucket="lbcu", 
                               fileFormat = "CSV",
                               fileNamePrefix = files$id[i])
task_vector$start()
ee_gcs_to_local(task = task_vector, dsn=paste0("gis/output/", files$id[i], ".csv"))
val <- read.csv(paste0("gis/output/", files$id[i], ".csv")) %>% 
  dplyr::select(-.geo)



#Proportion grassland

#Proportion crop

#Proportion developed

#Grassland to crop conversion

#Grassland to developed (tree/urban) conversion

#Surface water

#Center pivot

#Drought

#Storm frequency

#Max temp

#Raven density


#9. Extract the value ----
pt <- colMeans(xyFromCell(r, which(r[]>0))) %>% 
  t() %>% 
  data.frame() %>% 
  st_as_sf(coords = c('x','y'), crs = 4326)
pt.g <- sf_as_ee(pt)

gsw.occ.val <- ee_extract(
  x = gsw.occ.i,
  y = pt
)
