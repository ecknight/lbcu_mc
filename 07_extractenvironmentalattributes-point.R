library(sf)
library(tidyverse)
library(lubridate)
library(rgee)
library(data.table)
library(terra)

#1. Import & wrangle data----
dat.clust <- read.csv("Data/LBCUKDEClusters.csv") %>% 
  dplyr::filter(nclust==2,
                !is.na(X)) %>% 
  group_by(id, season, nclust) %>% 
  summarize(X = mean(X, na.rm=FALSE), 
            Y = mean(Y, na.rm=FALSE),
            kdecluster = round(mean(kdecluster))) %>% 
  ungroup()

dat <- read.csv("Data/LBCU_FilteredData_Segmented.csv") %>% 
  dplyr::filter(id %in% dat.clust$id,
                !(season %in% c("springmig", "fallmig") & stopover==0),
                winter!="wintermig") %>% 
  mutate(seasoncluster = ifelse(winter!="other", str_sub(winter, -1, -1), stopovercluster),
         seasoncluster = ifelse(is.na(seasoncluster), 0, as.numeric(seasoncluster)),
         kdeid = paste(season, id, seasoncluster, sep="-")) %>% 
  left_join(dat.clust %>% 
              dplyr::select(id, kdecluster) %>% 
              unique()) %>% 
  group_by(kdeid) %>% 
  mutate(n=n()) %>% 
  ungroup() %>% 
  arrange(X, Y) %>% 
  mutate(row = row_number(),
         n = ceiling(row/1000))

#2. Initialize rgee----
#ee_install()
ee_Initialize(gcs=TRUE)
ee_check()

#3. Get global surface water layer (1984 to 2021)----
gsw.occ<-ee$Image('JRC/GSW1_3/GlobalSurfaceWater')$select('occurrence')
gsw.chg<-ee$Image('JRC/GSW1_3/GlobalSurfaceWater')$select('change_norm')
gsw.mths<-ee$Image('JRC/GSW1_3/GlobalSurfaceWater')$select('seasonality')
gsw.yrs<-ee$Image('JRC/GSW1_3/GlobalSurfaceWater')$select('recurrence')

#4. Get drought data----
tclim.drght <- ee$ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')$filter(ee$Filter$date('1984-01-01', '2021-12-31'))$select('pdsi')$mean()

#5. Get landcover data----
dw.grass <- ee$ImageCollection("GOOGLE/DYNAMICWORLD/V1")$filter(ee$Filter$date('2015-01-01', '2021-12-31'))$select('grass')$mean()
dw.wetland <- ee$ImageCollection("GOOGLE/DYNAMICWORLD/V1")$filter(ee$Filter$date('2015-01-01', '2021-12-31'))$select('flooded_vegetation')$mean()
dw.crop <- ee$ImageCollection("GOOGLE/DYNAMICWORLD/V1")$filter(ee$Filter$date('2015-01-01', '2021-12-31'))$select('crops')$mean()
dw.built <- ee$ImageCollection("GOOGLE/DYNAMICWORLD/V1")$filter(ee$Filter$date('2015-01-01', '2021-12-31'))$select('built')$mean()

#6. Get landcover change data----
modis1 <- ee$ImageCollection("MODIS/006/MCD12Q1")$filter(ee$Filter$date('2001-01-01', '2001-12-31'))$select('LC_Type1')$mean()
modis2 <- ee$ImageCollection("MODIS/006/MCD12Q1")$filter(ee$Filter$date('2020-01-01', '2020-12-31'))$select('LC_Type1')$mean()

#7. Stack them----
stack <- gsw.occ$addBands(gsw.chg)$addBands(gsw.mths)$addBands(gsw.yrs)$addBands(tclim.drght)$addBands(dw.grass)$addBands(dw.wetland)$addBands(dw.crop)$addBands(dw.built)$addBands(modis1)$addBands(modis2)

#8. Set up loops----
loops <- max(dat$n)

data.cov <- list()
for(i in 1:loops){
  
  start_time <- Sys.time()
  
  dat.i <- dat %>% 
    dplyr::filter(n==i)
  
  #9. Create sf object----
  datasf <- st_as_sf(dat.i, coords = c('X','Y'), crs = 3857) %>% 
    st_transform(crs=4326)
  
  #10. Send data to GEE----
  data <- sf_as_ee(datasf)
  
  #11. Extract covariate point values----
  data.cov[[i]] <- ee_extract(
    x = stack,
    y = datasf,
    scale = 500,
    sf = FALSE
  )
  
  end_time <- Sys.time()
  
  print(paste0("Finished loop ", i, " of ", loops, " in ", end_time - start_time, " minutes"))

}

#12. Extract raven density----
ebd <- rast("/Users/ellyknight/Library/Application Support/ebirdst/comrav-ERD2019-STATUS-20200930-21f8da4a/abundance_seasonal/comrav-ERD2019-STATUS-20200930-21f8da4a_hr_2019_abundance-seasonal_resident.tif")

raven <- rbindlist(data.cov, fill=TRUE) %>% 
  st_as_sf(coords=c("lon", "lat"), crs=4326) %>% 
  st_transform(crs=crs(ebd)) %>% 
  st_coordinates() %>% 
  terra::extract(x=ebd)
colnames(raven) <- "raven"

#13. Collapse & wrangle----
cov <- rbindlist(data.cov, fill=TRUE) %>% 
  cbind(raven) %>% 
  rename(modis1 = LC_Type1, modis2 = LC_Type1_1) %>% 
  mutate(mgrass1 = ifelse(modis1==10, 1, 0),
         mcrop2 = ifelse(modis2==11, 1, 0),
         mbuilt2 = ifelse(modis2==12, 1, 0),
         convcrop = ifelse((mgrass1 - mcrop2)==1, 1, 0),
         convbuilt = ifelse((mgrass1 - mbuilt2)==1, 1, 0),
         change_norm = ifelse(is.na(change_norm), 0, change_norm),
         occurrence = ifelse(is.na(occurrence), 0, occurrence),
         recurrence = ifelse(is.na(recurrence), 0, recurrence),
         seasonality = ifelse(is.na(seasonality), 0, seasonality),
         raven = ifelse(is.na(raven), 0, raven)) 

#14. Save----
write.csv(cov, "Data/LBCU_environvars_pt.csv", row.names = FALSE)
