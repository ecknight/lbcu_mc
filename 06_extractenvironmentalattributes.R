library(sf)
library(raster)
library(tidyverse)
library(lubridate)
library(rgee)
library(data.table)
library(adehabitatHR)

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
kde <- read_sf("gis/shp/kde_individual.shp") %>% 
  st_make_valid() %>% 
  st_centroid() %>% 
  mutate(year = as.numeric(year)-1) %>% 
  dplyr::filter(nclust=="3")

#2. Set up to loop through years----
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
  dat <- read.csv(paste0("gis/output_ind/", files[i])) %>%
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

#16. Save-----
write.csv(covs, "Data/LBCU_environvars.csv", row.names = FALSE)

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

#C. RSF - USED AND AVAILABLE####

#1. Set up loop for combo of season * group----
loop <- kde %>% 
  data.frame() %>% 
  dplyr::select(group, season) %>% 
  unique()

#2. MCP of home ranges----
mcp.list <- list()
for(i in 1:nrow(loop)){
  
  kde.i <- kde %>% 
    dplyr::filter(group==loop$group[i],
                  season==loop$season[i]) %>% 
    mutate(id = paste0(group, "-", season)) %>% 
    dplyr::select(id) %>% 
    as_Spatial()
  
  mcp.list[[i]] <- mcp(kde.i) %>% 
    st_as_sf()
  
}
mcp <- do.call(rbind, mcp.list) %>% 
  separate(id, into=c("group", "season"))
plot(mcp)

#3. Available points----
set.seed(1234)
avail.list.list <- list()
for(i in 1:nrow(loop)){
  
  n.i <- kde %>% 
    dplyr::filter(group==loop$group[i],
                  season==loop$season[i]) %>% 
    group_by(year) %>% 
    summarize(n=n()) %>% 
    ungroup()
  
  avail.list <- list()
  for(j in 1:nrow(n.i)){
    avail.list[[j]] <- st_sample(mcp[i,], n.i$n[j]*10) %>% 
      st_coordinates() %>% 
      data.frame() %>%
      mutate(group=loop$group[i],
             season=loop$season[i],
             year=n.i$year[j])
  }
  
  avail.list.list[[i]] <- do.call(rbind, avail.list)

}
avail <- do.call(rbind, avail.list.list)

ggplot() +
  geom_sf(data=mcp, aes(fill=group)) +
  geom_point(data=avail, aes(x=X, y=Y, colour=group)) +
  facet_wrap(~season)

#4. Put used & available together----
pts <- avail %>% 
  mutate(type = "avail",
         id = NA) %>% 
  rbind(kde %>% 
          st_coordinates() %>% 
          cbind(kde) %>% 
          data.frame() %>% 
          dplyr::select(group, season, year, id, X, Y) %>% 
          mutate(type = "used")) %>% 
  arrange(year, season, group) %>% 
  group_by(year, season, group) %>% 
  mutate(index = row_number()-1) %>% 
  ungroup()

#5. Define radius (mean HR size)----
rad <- kde %>% 
  data.frame() %>% 
  group_by(group, season) %>% 
  summarize(rad = mean(rad)) %>% 
  ungroup()

#7. Set up to loop through years----
years <- sort(unique(pts$year))
for(i in 1:length(years)){
  
  pts.i <- pts %>% 
    dplyr::filter(year==years[i])
  
  #8. Get drought data----
  tclim.drght <- ee$ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')$filter(ee$Filter$date(paste0(years[i], "-01-01"), paste0(years[i], "-12-31")))$select('pdsi')$mean()
  
  #9. Get MODIS landcover----
  modis <- ee$ImageCollection("MODIS/006/MCD12Q1")$filter(ee$Filter$date(paste0(years[i], "-01-01"), paste0(years[i], "-12-31")))$select('LC_Type1')$mean()
  
  modis.grass <- modis$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0))
  modis.crop <- modis$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0))
  modis.built <- modis$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0))
  modis.wetland <- ee$ImageCollection("MODIS/006/MCD12Q1")$filter(ee$Filter$date(paste0(years[i], "-01-01"), paste0(years[i], "-12-31")))$select('LC_Prop3')$mean()$remap(c(1,2,3,10,20,27,30,40,50), c(0,0,0,0,0,0,0,1,0))
  
  #10. Get global surface water layer (1984 to 2021)----
  gsw.chg<-ee$Image('JRC/GSW1_3/GlobalSurfaceWater')$select('change_norm')
  gsw.mths<-ee$Image('JRC/GSW1_3/GlobalSurfaceWater')$select('seasonality')
  
  #11. Stack them----
  stack <- gsw.chg$addBands(gsw.mths)$addBands(tclim.drght)$addBands(modis.crop)$addBands(modis.built)$addBands(modis.grass)$addBands(modis.wetland)$rename("change", "seasonality", "drought", "crop", "built", "grass", "wetland")
  
  #12. Set up loop for season * group----
  for(j in 1:nrow(loop)){
    
    #13. Filter data again----
    pts.j <- pts.i %>% 
      dplyr::filter(season==loop$season[j],
                    group==loop$group[j])
    
    if(nrow(pts.j) > 0){
      
      #14. Send to gee----
      data <- pts.j %>% 
        st_as_sf(coords=c("X", "Y"), crs=4326) %>% 
        sf_as_ee()
      
      #15. Buffer----
      buffer_points <- function(feature){
        properties <- c("Date", "ID", "Radius", "Type", "date_millis", "ptID", "ptIDn", "X", "Y", "timestamp", "uniq")
        new_geometry <- feature$geometry()$buffer(rad$rad[j]*1000)
        ee$Feature(new_geometry)$copyProperties(feature, properties)
      }
      
      data.buff <- data$map(buffer_points)
      
      #16. Extract mean values---
      stack.mn <- stack$reduceRegions(reducer=ee$Reducer$mean(),
                                      collection=data.buff,
                                      scale=30)
      
      #17. Export the values----
      task_vector <- ee_table_to_gcs(collection=stack.mn,
                                     bucket="lbcu",
                                     fileFormat = "CSV",
                                     fileNamePrefix = kde.j$id)
      task_vector$start()
      ee_monitoring(task_vector, max_attempts=1000)
      
      #18. Download----
      ee_gcs_to_local(task = task_vector, dsn=paste0("gis/output_rsf/", years[i], "_", loop$season[j], "_", loop$group[j], ".csv"))
      
    }
    
    print(paste0("Finished loop ", j, " of ", nrow(loop), " for year ", years[i]))
  }
  
}

#19. Read in output----
files <- list.files("gis/output_rsf")

dat <- data.frame()
for(i in 1:length(files)){
  dat <- read.csv(paste0("gis/output_rsf/", files[i])) %>%
    dplyr::select(-.geo) %>% 
    mutate(id = str_sub(files[i], -100, -5)) %>% 
    rename(index = system.index) %>% 
    rbind(dat)
}

#14. Join to other data and tidy----
covs <- dat %>% 
  separate(id, into=c("year", "season", "group"), remove=TRUE) %>% 
  mutate(year = as.numeric(year)) %>% 
  arrange(year, season, group, index) %>% 
  full_join(pts) %>% 
  mutate(change = ifelse(is.na(change), 0, change),
         seasonality = ifelse(is.na(seasonality), 0, seasonality),
         crop = ifelse(is.na(crop), 0, crop),
         drought = ifelse(is.na(drought), 0, drought),
         grass = ifelse(is.na(grass), 0, grass),
         wetland = ifelse(is.na(wetland), 0, wetland),
         built = ifelse(is.na(built), 0, built))

#16. Save-----
write.csv(covs, "Data/LBCU_environvars_RSF.csv", row.names = FALSE)
