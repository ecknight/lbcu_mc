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

#3. Set wd----
setwd("G:/My Drive/SMBC")

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
dat.out <- list()
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
  stack.mn <- list()
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
    stack.mn[[j]] <- ee_extract(stack, data.buff, fun=ee$Reducer$mean(), scale = 30)
    
    print(paste0("Finished individual ", j, " of ", nrow(kde.i), " for year ", years[i]))
  }
  
  dat.out[[i]] <- data.table::rbindlist(stack.mn, fill=TRUE)
  
}

#13. Read in output----

dat <- do.call(rbind, dat.out)

#14. Join to other data and tidy----
covs <- cbind(dat, kde)  %>% 
  st_drop_geometry() |> 
  left_join(kde %>% 
              st_coordinates() %>% 
              data.frame() %>% 
              cbind(data.frame(kde)) %>% 
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

#2. Set up to loop through nclust---
clusts <- unique(kde$nclust)[c(1:3)]

dat.out <- list()
for(h in 1:length(clusts)){
  
  #3. Filter data----
  kde.h <- kde |> 
    dplyr::filter(nclust==clusts[h])
  
  #4. Send to gee---
  data <- sf_as_ee(kde.h)
  
  #5. Set up to loop through years----
  years <- c(2006:2020)
  stack.mn <- list()
  for(i in 1:length(years)){
    
    #6. Get drought data----
    tclim.drght <- ee$ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')$filter(ee$Filter$date(paste0(years[i], "-01-01"), paste0(years[i], "-12-31")))$select('pdsi')$mean()
    
    #7. Get MODIS landcover----
    modis <- ee$ImageCollection("MODIS/006/MCD12Q1")$filter(ee$Filter$date(paste0(years[i], "-01-01"), paste0(years[i], "-12-31")))$select('LC_Type1')$mean()
    
    modis.grass <- modis$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0))
    modis.crop <- modis$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0))
    modis.built <- modis$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0))
    modis.wetland <- ee$ImageCollection("MODIS/006/MCD12Q1")$filter(ee$Filter$date(paste0(years[i], "-01-01"), paste0(years[i], "-12-31")))$select('LC_Prop3')$mean()$remap(c(1,2,3,10,20,27,30,40,50), c(0,0,0,0,0,0,0,1,0))
    
    #8. Get global surface water layer (1984 to 2021)----
    gsw.chg<-ee$Image('JRC/GSW1_3/GlobalSurfaceWater')$select('change_norm')
    gsw.mths<-ee$Image('JRC/GSW1_3/GlobalSurfaceWater')$select('seasonality')
    
    #9. Stack them----
    stack <- gsw.chg$addBands(gsw.mths)$addBands(tclim.drght)$addBands(modis.crop)$addBands(modis.built)$addBands(modis.grass)$addBands(modis.wetland)$rename("change", "seasonality", "drought", "crop", "built", "grass", "wetland")
    
    #10. Extract mean values---
    stack.mn[[i]] <- ee_extract(stack, data,
                           fun = ee$Reducer$mean(),
                           scale = 1000)

    print(paste0("Finished year ", i, " of ", length(years)))
    
  }
  
  dat.out[[h]] <- do.call(rbind, stack.mn)
  
  print(paste0("Finished nclust ", h, " of ", length(clusts)))
  
}

#11. Collapse----
dat <- do.call(rbind, dat.out)

#14. Save-----
write.csv(dat, "Data/LBCU_environvars_region.csv", row.names = FALSE)

#C. RSF - USED AND AVAILABLE####

#1. Read in shapefile of HRs and get centroid----
#just use one nclust group, no need to do twice
kde <- read_sf("gis/shp/kde_individual.shp") %>% 
  st_make_valid() %>% 
  st_centroid() %>% 
  mutate(year = as.numeric(year)-1)

#1. Set up loop for combo of season * group----
loop <- kde %>% 
  data.frame() %>% 
  dplyr::filter(nclust %in% c("3", "expert", "flyway")) |> 
  dplyr::select(nclust, group, season) %>% 
  unique()

#2. MCP of home ranges----
mcp.list <- list()
for(i in 1:nrow(loop)){
  
  kde.i <- kde %>% 
    dplyr::filter(group==loop$group[i],
                  season==loop$season[i],
                  nclust==loop$nclust[i]) %>% 
    mutate(id = paste0(nclust, "-", group, "-", season)) %>% 
    dplyr::select(id) %>% 
    as_Spatial()
  
  mcp.list[[i]] <- mcp(kde.i) %>% 
    st_as_sf()
  
}

mcp <- do.call(rbind, mcp.list) %>% 
  separate(id, into=c("nclust", "group", "season")) |> 
  st_make_valid()
plot(mcp) 

#3. Available points----
set.seed(1234)
avail.list.list <- list()
for(i in 1:nrow(loop)){
  
  n.i <- kde %>% 
    dplyr::filter(group==loop$group[i],
                  season==loop$season[i],
                  nclust==loop$nclust[i]) %>% 
    group_by(year) %>% 
    summarize(n=n()) %>% 
    ungroup()
  
  avail.list <- list()
  for(j in 1:nrow(n.i)){
    
    ids.j <- kde %>% 
      dplyr::filter(group==loop$group[i],
                    season==loop$season[i],
                    nclust==loop$nclust[i],
                    year==n.i$year[j])
    
    avail.list[[j]] <- st_sample(mcp[i,], n.i$n[j]*10) %>% 
      st_coordinates() %>% 
      data.frame() %>%
      mutate(group=loop$group[i],
             season=loop$season[i],
             nclust=loop$nclust[i],
             year=n.i$year[j]) %>% 
      cbind(data.frame(bird = rep(ids.j$bird)))
  }
  
  avail.list.list[[i]] <- do.call(rbind, avail.list)

}
avail <- do.call(rbind, avail.list.list)

ggplot() +
  geom_sf(data=mcp, aes(fill=group)) +
  geom_point(data=avail, aes(x=X, y=Y, colour=group)) +
  facet_grid(nclust~season)

#4. Put used & available together----
pts <- avail %>% 
  mutate(type = "avail") %>% 
  rbind(kde %>% 
          st_coordinates() %>% 
          cbind(kde) %>% 
          data.frame() %>% 
          dplyr::select(nclust, group, season, year, bird, X, Y) %>% 
          mutate(type = "used")) %>% 
  arrange(year, season, group) %>% 
  group_by(year, season, group) %>% 
  mutate(index = row_number()-1) %>% 
  ungroup()

#5. Define radius (mean HR size)----
rad <- kde %>% 
  data.frame() %>% 
  group_by(nclust, group, season) %>% 
  summarize(rad = mean(rad)) %>% 
  ungroup()

#6. Calculate mean months for each season----
raw <- read.csv("Data/LBCU_FilteredData_Segmented.csv") %>% 
  dplyr::filter(!id %in% c(46768277, 33088, 129945787, 46770723, 46769927, 86872)) %>% 
  dplyr::filter(segment %in% c("depart", "arrive")) %>% 
  mutate(season = case_when(season=="breed" ~ "springmig",
                            season=="winter" ~ "fallmig",
                            !is.na(season) ~ season)) %>% 
  group_by(segment, season) %>% 
  summarize(mean = mean(doy)) %>% 
  ungroup()

#7. Set up to loop through year*season*group----
years <- pts %>% 
  dplyr::select(year, season, group, nclust) %>% 
  unique() %>% 
  mutate(month.start = case_when(season=="springmig" ~ "02",
                                 season=="breed" ~ "04",
                                 season=="fallmig" ~ "06",
                                 season=="winter" ~ "07"),
         month.end = case_when(season=="springmig" ~ "03",
                               season=="breed" ~ "06",
                               season=="fallmig" ~ "07",
                               season=="winter" ~ "02"),
         year.start = year,
         year.end = ifelse(season=="winter", year+1, year)) %>% 
  left_join(rad)

stack.list <- list()
for(i in 1:nrow(years)){
  
  #8. Filter points----
  pts.i <- pts %>% 
    dplyr::filter(year==years$year[i],
                  season==years$season[i],
                  nclust==years$nclust[i],
                  group==years$group[i])
  
  #9. Send to gee----
  data <- pts.i %>% 
    st_as_sf(coords=c("X", "Y"), crs=4326) %>% 
    sf_as_ee()
  
  #10. Get drought data----
  tclim.drght <- ee$ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')$filter(ee$Filter$date(paste0(years$year.start[i], "-", years$month.start[i], "-01"), paste0(years$year.end[i], "-", years$month.end[i], "-28")))$select('pdsi')$mean()
  
  #11. Get MODIS landcover----
  modis <- ee$ImageCollection("MODIS/006/MCD12Q1")$filter(ee$Filter$date(paste0(years$year.start[i], "-01-01"), paste0(years$year.end[i], "-12-31")))$select('LC_Type1')$mean()
  
  modis.grass <- modis$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0))
  modis.crop <- modis$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0))
  modis.built <- modis$remap(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0))
  modis.wetland <- ee$ImageCollection("MODIS/006/MCD12Q1")$filter(ee$Filter$date(paste0(years$year.start[i], "-01-01"), paste0(years$year.end[i], "-12-31")))$select('LC_Prop3')$mean()$remap(c(1,2,3,10,20,27,30,40,50), c(0,0,0,0,0,0,0,1,0))
  
  #12. Get global surface water layer (1984 to 2021)----
  gsw.chg<-ee$Image('JRC/GSW1_3/GlobalSurfaceWater')$select('change_norm')
  gsw.mths<-ee$Image('JRC/GSW1_3/GlobalSurfaceWater')$select('seasonality')
  
  #13. Stack them----
  stack <- gsw.chg$addBands(gsw.mths)$addBands(tclim.drght)$addBands(modis.crop)$addBands(modis.built)$addBands(modis.grass)$addBands(modis.wetland)$rename("change", "seasonality", "drought", "crop", "built", "grass", "wetland")
  
  #14. Buffer----
  buffer_points <- function(feature){
    properties <- c("Date", "ID", "Radius", "Type", "date_millis", "ptID", "ptIDn", "X", "Y", "timestamp", "uniq")
    new_geometry <- feature$geometry()$buffer(years$rad[i]*1000)
    ee$Feature(new_geometry)$copyProperties(feature, properties)
  }
  
  data.buff <- data$map(buffer_points)
  
  #15. Extract mean values---
  stack.list[[i]] <- ee_extract(stack, data.buff,
                              fun = ee$Reducer$mean(),
                              scale = 30) |> 
    cbind(pts.i)
  
  print(paste0("Finished loop ", i, " of ", nrow(years)))
  
}

#19. Read in output----
dat <- data.table::rbindlist(stack.list, fill=TRUE)

#14. Join to other data and tidy----
covs <- dat %>% 
  mutate(change = ifelse(is.na(change), 0, change),
         seasonality = ifelse(is.na(seasonality), 0, seasonality),
         crop = ifelse(is.na(crop), 0, crop),
         drought = ifelse(is.na(drought), 0, drought),
         grass = ifelse(is.na(grass), 0, grass),
         wetland = ifelse(is.na(wetland), 0, wetland),
         built = ifelse(is.na(built), 0, built))

#16. Save-----
write.csv(covs, "Data/LBCU_environvars_RSF.csv", row.names = FALSE)
