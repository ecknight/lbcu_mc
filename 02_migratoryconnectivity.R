library(tidyverse)
library(data.table)
library(vegan)
library(adehabitatHR)
library(sf)
library(sp)
library(raster)
library(MigConnectivity)
library(ebirdst)

#1. Load clustered data----
dat <- read.csv("Data/LBCUKDEClusters.csv")

#2. Extract relative abundance information from eBird----
#set_ebirdst_access_key("e7ld1bagh8n1")
#ebirdst_download("lobcur")
ebd <- load_raster("/Users/ellyknight/Library/Application Support/ebirdst/lobcur-ERD2019-STATUS-20200930-da308f90", product="abundance_seasonal")

dat.sf <- dat %>% 
  st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
  st_transform("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs")

dat.abun <- raster::extract(ebd, dat.sf, na.rm=TRUE) %>% 
  data.frame() %>% 
  cbind(dat)

abun <- dat.abun %>% 
  group_by(nclust, group, season) %>% 
  summarize(breeding = mean(breeding, na.rm=TRUE),
            postbreeding_migration = mean(postbreeding_migration, na.rm=TRUE),
            nonbreeding = mean(nonbreeding, na.rm=TRUE),
            prebreeding_migration = mean(prebreeding_migration, na.rm=TRUE),
            n=n()) %>% 
  ungroup()

#3. Set up loop through # of clusters---
clusters <- unique(dat$nclust)

mc.df <- data.frame()
for(j in 1:length(clusters)){
  
  #4. Wrangle data----
  breed.j <- dat %>% 
    dplyr::filter(season=="breed",
                  nclust==clusters[j])
  
  winter.j <- dat %>% 
    dplyr::filter(season=="winter",
                  nclust==clusters[j])
  
  fall.j <- dat %>% 
    dplyr::filter(season=="fallmig",
                  nclust==clusters[j])
  
  spring.j <- dat %>% 
    dplyr::filter(season=="springmig",
                  nclust==clusters[j])
  
  bw.j <- rbind(breed.j, winter.j) %>% 
    dplyr::select(id, group, X, Y, season) %>% 
    pivot_wider(id_cols=id:group, names_from=season, values_from=X:Y) %>% 
    dplyr::filter(!is.na(X_winter) & !is.na(X_breed)) %>% 
    rename(bird=id)
  
  bf.j <- rbind(breed.j, fall.j) %>% 
    dplyr::select(id, group, X, Y, season) %>% 
    pivot_wider(id_cols=id:group, names_from=season, values_from=X:Y) %>% 
    dplyr::filter(!is.na(X_fallmig) & !is.na(X_breed)) %>% 
    rename(bird=id)
  
  wf.j <- rbind(winter.j, fall.j) %>% 
    dplyr::select(id, group, X, Y, season) %>% 
    pivot_wider(id_cols=id:group, names_from=season, values_from=X:Y) %>% 
    dplyr::filter(!is.na(X_fallmig) & !is.na(X_winter)) %>% 
    rename(bird=id)
  
  bs.j <- rbind(breed.j, spring.j) %>% 
    dplyr::select(id, group, X, Y, season) %>% 
    pivot_wider(id_cols=id:group, names_from=season, values_from=X:Y) %>% 
    dplyr::filter(!is.na(X_springmig) & !is.na(X_breed)) %>% 
    rename(bird=id)
  
  ws.j <- rbind(winter.j, spring.j) %>% 
    dplyr::select(id, group, X, Y, season) %>% 
    pivot_wider(id_cols=id:group, names_from=season, values_from=X:Y) %>% 
    dplyr::filter(!is.na(X_springmig) & !is.na(X_winter)) %>% 
    rename(bird=id)
  
  #5. Set up target grids for MC estimation----
  
  #Breed:Winter
  ptsbw <- st_as_sf(bw.j, coords=c("X_winter", "Y_winter"), crs=3857) %>% 
    as_Spatial()
  
  bb <- bbox(ptsbw)
  cs <- c(500000, 500000)  # cell size 
  cc <- bb[, 1] + (cs/2)  # cell offset
  cd <- ceiling(diff(t(bb))/cs)  # number of cells per direction
  grd <- GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd)
  spgrdbw <- SpatialGridDataFrame(grd,
                                  data=data.frame(id=1:prod(cd)),
                                  proj4string=CRS(proj4string(ptsbw)))
  
  idbw.j <- bw.j %>% 
    rename(breed_id=group) %>% 
    cbind(id=over(ptsbw, spgrdbw)) %>% 
    dplyr::rename(winter_id=id)
  
  #Winter:Breed
  ptswb <- st_as_sf(bw.j, coords=c("X_breed", "Y_breed"), crs=3857) %>% 
    as_Spatial()
  
  bb <- bbox(ptswb)
  cs <- c(500000, 500000)  # cell size 
  cc <- bb[, 1] + (cs/2)  # cell offset
  cd <- ceiling(diff(t(bb))/cs)  # number of cells per direction
  grd <- GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd)
  spgrdwb <- SpatialGridDataFrame(grd,
                                  data=data.frame(id=1:prod(cd)),
                                  proj4string=CRS(proj4string(ptswb)))
  
  idwb.j <- bw.j %>% 
    rename(winter_id=group) %>% 
    cbind(id=over(ptswb, spgrdwb)) %>% 
    dplyr::rename(breed_id=id)
  
  #Breed:Fall
  ptsbf <- st_as_sf(bf.j, coords=c("X_fallmig", "Y_fallmig"), crs=3857) %>% 
    as_Spatial()
  
  bb <- bbox(ptsbf)
  cs <- c(500000, 500000)  # cell size 
  cc <- bb[, 1] + (cs/2)  # cell offset
  cd <- ceiling(diff(t(bb))/cs)  # number of cells per direction
  grd <- GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd)
  spgrdbf <- SpatialGridDataFrame(grd,
                                  data=data.frame(id=1:prod(cd)),
                                  proj4string=CRS(proj4string(ptsbf)))
  
  idbf.j <- bf.j %>% 
    rename(breed_id=group) %>% 
    cbind(id=over(ptsbf, spgrdbf)) %>% 
    dplyr::rename(fallmig_id=id)
  
  #Winter:Fall
  ptswf <- st_as_sf(wf.j, coords=c("X_fallmig", "Y_fallmig"), crs=3857) %>% 
    as_Spatial()
  
  bb <- bbox(ptswf)
  cs <- c(500000, 500000)  # cell size 
  cc <- bb[, 1] + (cs/2)  # cell offset
  cd <- ceiling(diff(t(bb))/cs)  # number of cells per direction
  grd <- GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd)
  spgrdwf <- SpatialGridDataFrame(grd,
                                  data=data.frame(id=1:prod(cd)),
                                  proj4string=CRS(proj4string(ptswf)))
  
  idwf.j <- wf.j %>% 
    rename(winter_id=group) %>% 
    cbind(id=over(ptswf, spgrdwf)) %>% 
    dplyr::rename(fallmig_id=id)
  
  #Breed:Spring
  ptsbs <- st_as_sf(bs.j, coords=c("X_springmig", "Y_springmig"), crs=3857) %>% 
    as_Spatial()
  
  bb <- bbox(ptsbs)
  cs <- c(500000, 500000)  # cell size 
  cc <- bb[, 1] + (cs/2)  # cell offset
  cd <- ceiling(diff(t(bb))/cs)  # number of cells per direction
  grd <- GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd)
  spgrdbs <- SpatialGridDataFrame(grd,
                                  data=data.frame(id=1:prod(cd)),
                                  proj4string=CRS(proj4string(ptsbs)))
  
  idbs.j <- bs.j %>% 
    rename(breed_id=group) %>% 
    cbind(id=over(ptsbs, spgrdbs)) %>% 
    dplyr::rename(springmig_id=id)
  
  #Winter:Spring
  ptsws <- st_as_sf(ws.j, coords=c("X_springmig", "Y_springmig"), crs=3857) %>% 
    as_Spatial()
  
  bb <- bbox(ptsws)
  cs <- c(500000, 500000)  # cell size 
  cc <- bb[, 1] + (cs/2)  # cell offset
  cd <- ceiling(diff(t(bb))/cs)  # number of cells per direction
  grd <- GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd)
  spgrdws <- SpatialGridDataFrame(grd,
                                  data=data.frame(id=1:prod(cd)),
                                  proj4string=CRS(proj4string(ptsws)))
  
  idws.j <- ws.j %>% 
    rename(winter_id=group) %>% 
    cbind(id=over(ptsws, spgrdws)) %>% 
    dplyr::rename(springmig_id=id)
  
  #6. Distance matrices between centroids----
  
  #breed:winter
  distbwb.j <- idbw.j %>% 
    group_by(breed_id) %>% 
    summarize(X=mean(X_breed),
              Y=mean(Y_breed)) %>% 
    ungroup() %>% 
    arrange(breed_id) %>% 
    dplyr::select(X, Y) %>% 
    data.frame() %>% 
    distFromPos("plane")
  dimnames(distbwb.j) <- list(sort(unique(idbw.j$breed_id)), sort(unique(idbw.j$breed_id)))
  
  distbww.j <- idbw.j %>% 
    group_by(winter_id) %>% 
    summarize(X=mean(X_winter),
              Y=mean(Y_winter)) %>% 
    ungroup() %>% 
    arrange(winter_id) %>% 
    dplyr::select(X, Y) %>% 
    data.frame() %>% 
    distFromPos("plane")
  dimnames(distbww.j) <- list(sort(unique(idbw.j$winter_id)), sort(unique(idbw.j$winter_id)))
  
  #winter:breed
  distwbw.j <- idwb.j %>% 
    group_by(winter_id) %>% 
    summarize(X=mean(X_winter),
              Y=mean(Y_winter)) %>% 
    ungroup() %>% 
    arrange(winter_id) %>% 
    dplyr::select(X, Y) %>% 
    data.frame() %>% 
    distFromPos("plane")
  dimnames(distwbw.j) <- list(sort(unique(idwb.j$winter_id)), sort(unique(idwb.j$winter_id)))
  
  distwbb.j <- idwb.j %>% 
    group_by(breed_id) %>% 
    summarize(X=mean(X_breed),
              Y=mean(Y_breed)) %>% 
    ungroup() %>% 
    arrange(breed_id) %>% 
    dplyr::select(X, Y) %>% 
    data.frame() %>% 
    distFromPos("plane")
  dimnames(distwbb.j) <- list(sort(unique(idwb.j$breed_id)), sort(unique(idwb.j$breed_id)))
  
  #breed:fall
  distbfb.j <- idbf.j %>% 
    group_by(breed_id) %>% 
    summarize(X=mean(X_breed),
              Y=mean(Y_breed)) %>% 
    ungroup() %>% 
    arrange(breed_id) %>% 
    dplyr::select(X, Y) %>% 
    data.frame() %>% 
    distFromPos("plane")
  dimnames(distbfb.j) <- list(sort(unique(idbf.j$breed_id)), sort(unique(idbf.j$breed_id)))
  
  distbff.j <- idbf.j %>% 
    group_by(fallmig_id) %>% 
    summarize(X=mean(X_fallmig),
              Y=mean(Y_fallmig)) %>% 
    ungroup() %>% 
    arrange(fallmig_id) %>% 
    dplyr::select(X, Y) %>% 
    data.frame() %>% 
    distFromPos("plane")
  dimnames(distbff.j) <- list(sort(unique(idbf.j$fallmig_id)), sort(unique(idbf.j$fallmig_id)))
  
  #winter:fall
  distwfw.j <- idwf.j %>% 
    group_by(winter_id) %>% 
    summarize(X=mean(X_winter),
              Y=mean(Y_winter)) %>% 
    ungroup() %>% 
    arrange(winter_id) %>% 
    dplyr::select(X, Y) %>% 
    data.frame() %>% 
    distFromPos("plane")
  dimnames(distwfw.j) <- list(sort(unique(idwf.j$breed_id)), sort(unique(idwf.j$breed_id)))
  
  distwff.j <- idwf.j %>% 
    group_by(fallmig_id) %>% 
    summarize(X=mean(X_fallmig),
              Y=mean(Y_fallmig)) %>% 
    ungroup() %>% 
    arrange(fallmig_id) %>% 
    dplyr::select(X, Y) %>% 
    data.frame() %>% 
    distFromPos("plane")
  dimnames(distwff.j) <- list(sort(unique(idwf.j$fallmig_id)), sort(unique(idwf.j$fallmig_id)))
  
  #breed:spring
  distbsb.j <- idbs.j %>% 
    group_by(breed_id) %>% 
    summarize(X=mean(X_breed),
              Y=mean(Y_breed)) %>% 
    ungroup() %>% 
    arrange(breed_id) %>% 
    dplyr::select(X, Y) %>% 
    data.frame() %>% 
    distFromPos("plane")
  dimnames(distbsb.j) <- list(sort(unique(idbs.j$breed_id)), sort(unique(idbs.j$breed_id)))
  
  distbss.j <- idbs.j %>% 
    group_by(springmig_id) %>% 
    summarize(X=mean(X_springmig),
              Y=mean(Y_springmig)) %>% 
    ungroup() %>% 
    arrange(springmig_id) %>% 
    dplyr::select(X, Y) %>% 
    data.frame() %>% 
    distFromPos("plane")
  dimnames(distbss.j) <- list(sort(unique(idbs.j$springmig_id)), sort(unique(idbs.j$springmig_id)))
  
  #winter:spring
  distwsw.j <- idws.j %>% 
    group_by(winter_id) %>% 
    summarize(X=mean(X_winter),
              Y=mean(Y_winter)) %>% 
    ungroup() %>% 
    arrange(winter_id) %>% 
    dplyr::select(X, Y) %>% 
    data.frame() %>% 
    distFromPos("plane")
  dimnames(distwsw.j) <- list(sort(unique(idws.j$breed_id)), sort(unique(idws.j$breed_id)))
  
  distwss.j <- idws.j %>% 
    group_by(springmig_id) %>% 
    summarize(X=mean(X_springmig),
              Y=mean(Y_springmig)) %>% 
    ungroup() %>% 
    arrange(springmig_id) %>% 
    dplyr::select(X, Y) %>% 
    data.frame() %>% 
    distFromPos("plane")
  dimnames(distwss.j) <- list(sort(unique(idws.j$springmig_id)), sort(unique(idws.j$springmig_id)))
  
  #7. Create point objects for individual locations----
  
  #Breed:Winter
  ptbwb.j <- idbw.j %>% 
    st_as_sf(coords=c("X_breed", "Y_breed"), crs=3857) %>% 
    dplyr::select(geometry)
  
  ptbww.j <- idbw.j %>% 
    st_as_sf(coords=c("X_winter", "Y_winter"), crs=3857) %>% 
    dplyr::select(geometry)
  
  #Winter:breed
  ptwbb.j <- idwb.j %>% 
    st_as_sf(coords=c("X_breed", "Y_breed"), crs=3857) %>% 
    dplyr::select(geometry)
  
  ptwbw.j <- idwb.j %>% 
    st_as_sf(coords=c("X_winter", "Y_winter"), crs=3857) %>% 
    dplyr::select(geometry)
  
  #breed:fall
  ptbfb.j <- idbf.j %>% 
    st_as_sf(coords=c("X_breed", "Y_breed"), crs=3857) %>% 
    dplyr::select(geometry)
  
  ptbff.j <- idbf.j %>% 
    st_as_sf(coords=c("X_fallmig", "Y_fallmig"), crs=3857) %>% 
    dplyr::select(geometry)
  
  #winter:fall
  ptwfw.j <- idwf.j %>% 
    st_as_sf(coords=c("X_winter", "Y_winter"), crs=3857) %>% 
    dplyr::select(geometry)
  
  ptwff.j <- idwf.j %>% 
    st_as_sf(coords=c("X_fallmig", "Y_fallmig"), crs=3857) %>% 
    dplyr::select(geometry)
  
  #breed:Spring
  ptbsb.j <- idbs.j %>% 
    st_as_sf(coords=c("X_breed", "Y_breed"), crs=3857) %>% 
    dplyr::select(geometry)
  
  ptbss.j <- idbs.j %>% 
    st_as_sf(coords=c("X_springmig", "Y_springmig"), crs=3857) %>% 
    dplyr::select(geometry)
  
  #winter:Spring
  ptwsw.j <- idws.j %>% 
    st_as_sf(coords=c("X_winter", "Y_winter"), crs=3857) %>% 
    dplyr::select(geometry)
  
  ptwss.j <- idws.j %>% 
    st_as_sf(coords=c("X_springmig", "Y_springmig"), crs=3857) %>% 
    dplyr::select(geometry)
  
  #8. Create origin region polygons----
  #Breed:Winter - breed
  nbw.j <- idbw.j %>% 
    group_by(breed_id) %>% 
    summarize(n=n()) %>% 
    ungroup() %>% 
    summarize(n = min(n))
  
  sitesbwb.j <- idbw.j %>% 
    group_by(breed_id) %>% 
    summarize(X=mean(X_breed),
              Y=mean(Y_breed)) %>% 
    ungroup() %>% 
    st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
    st_buffer(500000) %>% 
    st_cast("MULTIPOLYGON")
  
  #winter:breed - winter
  nwb.j <- idwb.j %>% 
    group_by(winter_id) %>% 
    summarize(n=n()) %>% 
    ungroup() %>% 
    summarize(n = min(n))
  
  siteswbw.j <- idwb.j %>% 
    group_by(winter_id) %>% 
    summarize(X=mean(X_winter),
              Y=mean(Y_winter)) %>% 
    ungroup() %>% 
    st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
    st_buffer(500000) %>% 
    st_cast("MULTIPOLYGON")
  
  #breed:Fall - breed
  nbf.j <- idbf.j %>% 
    group_by(breed_id) %>% 
    summarize(n=n()) %>% 
    ungroup() %>% 
    summarize(n = min(n))
  
  sitesbfb.j <- idbf.j %>% 
    group_by(breed_id) %>% 
    summarize(X=mean(X_breed),
              Y=mean(Y_breed)) %>% 
    ungroup() %>% 
    st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
    st_buffer(500000) %>% 
    st_cast("MULTIPOLYGON")
  
  #winter:Fall - winter
  nwf.j <- idwf.j %>% 
    group_by(winter_id) %>% 
    summarize(n=n()) %>% 
    ungroup() %>% 
    summarize(n = min(n))
  
  siteswfw.j <- idwf.j %>% 
    group_by(winter_id) %>% 
    summarize(X=mean(X_winter),
              Y=mean(Y_winter)) %>% 
    ungroup() %>% 
    st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
    st_buffer(500000) %>% 
    st_cast("MULTIPOLYGON")
  
  #breed:spring - breed
  nbs.j <- idbs.j %>% 
    group_by(breed_id) %>% 
    summarize(n=n()) %>% 
    ungroup() %>% 
    summarize(n = min(n))
  
  sitesbsb.j <- idbs.j %>% 
    group_by(breed_id) %>% 
    summarize(X=mean(X_breed),
              Y=mean(Y_breed)) %>% 
    ungroup() %>% 
    st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
    st_buffer(500000) %>% 
    st_cast("MULTIPOLYGON")
  
  #winter:spring - winter
  nws.j <- idws.j %>% 
    group_by(winter_id) %>% 
    summarize(n=n()) %>% 
    ungroup() %>% 
    summarize(n = min(n))
  
  siteswsw.j <- idws.j %>% 
    group_by(winter_id) %>% 
    summarize(X=mean(X_winter),
              Y=mean(Y_winter)) %>% 
    ungroup() %>% 
    st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
    st_buffer(500000) %>% 
    st_cast("MULTIPOLYGON")
  
  #9. Create target region polygons----
  
  #Breed:Winter - winter
  sitesbww.j <- spgrdbw %>% 
    raster() %>% 
    rasterToPolygons() %>% 
    st_as_sf() %>% 
    dplyr::filter(id %in% unique(idbw.j$winter_id)) %>% 
    dplyr::select(geometry) %>% 
    st_cast("MULTIPOLYGON")
  
  #winter:breed - breed
  siteswbb.j <- spgrdwb %>% 
    raster() %>% 
    rasterToPolygons() %>% 
    st_as_sf() %>% 
    dplyr::filter(id %in% unique(idwb.j$breed_id)) %>% 
    dplyr::select(geometry) %>% 
    st_cast("MULTIPOLYGON")
  
  #breed:Fall - fall
  sitesbff.j <- spgrdbf %>% 
    raster() %>% 
    rasterToPolygons() %>% 
    st_as_sf() %>% 
    dplyr::filter(id %in% unique(idbf.j$fallmig_id)) %>% 
    dplyr::select(geometry)
  
  #winter:Fall - fall
  siteswff.j <- spgrdwf %>% 
    raster() %>% 
    rasterToPolygons() %>% 
    st_as_sf() %>% 
    dplyr::filter(id %in% unique(idwf.j$fallmig_id)) %>% 
    dplyr::select(geometry)
  
  #breed:Spring - spring
  sitesbss.j <- spgrdbs %>% 
    raster() %>% 
    rasterToPolygons() %>% 
    st_as_sf() %>% 
    dplyr::filter(id %in% unique(idbs.j$springmig_id)) %>% 
    dplyr::select(geometry)
  
  #winter:Spring - spring
  siteswss.j <- spgrdws %>% 
    raster() %>% 
    rasterToPolygons() %>% 
    st_as_sf() %>% 
    dplyr::filter(id %in% unique(idws.j$springmig_id)) %>% 
    dplyr::select(geometry)
  
  #10. Add relative abundance----
  
  #Breed
  bbsb.j <- abun %>% 
    dplyr::filter(nclust==clusters[j], season=="breed") %>% 
    mutate(abun = breeding/(sum(breeding)))
  abunb.j <- bbsb.j$abun
  
  #Winter
  bbsw.j <- abun %>% 
    dplyr::filter(nclust==clusters[j], season=="winter") %>% 
    mutate(abun = nonbreeding/(sum(nonbreeding)))
  abunw.j <- bbsw.j$abun
  
  #11. Define tag type----
  telbw.j <- rep(TRUE, nrow(idbw.j))
  telbf.j <- rep(TRUE, nrow(idbf.j))
  telbs.j <- rep(TRUE, nrow(idbs.j))
  telwb.j <- rep(TRUE, nrow(idwb.j))
  telwf.j <- rep(TRUE, nrow(idwf.j))
  telws.j <- rep(TRUE, nrow(idws.j))
  
  #12. Estimate MC----
  
  #Breed:Winter
  if(class(sitesbwb.j)[1]=="sf"){
    set.seed(1234)
    transbw.j <- estTransition(originSites = sitesbwb.j, 
                               targetSites = sitesbww.j, 
                               originPoints = ptbwb.j,
                               targetPoints = ptbww.j,
                               originAssignment = as.numeric(as.factor(idbw.j$breed_id)),
                               targetAssignment = as.numeric(as.factor(idbw.j$winter_id)),
                               originNames = sort(unique(idbw.j$breed_id)),
                               targetNames = sort(unique(idbw.j$winter_id)),
                               nSamples = 1000,
                               isTelemetry = telbw.j)
    
    mcbw.j <- estStrength(originDist = distbwb.j,
                          targetDist = distbww.j,
                          originRelAbund = abunb.j,
                          psi = transbw.j,
                          sampleSize = nrow(idbw.j),
                          originNames = sort(unique(idbw.j$breed_id)),
                          targetNames = sort(unique(idbw.j$winter_id)),
                          nSamples = 1000)
  }
  
  mc.df <-  data.frame(MC = mcbw.j$MC$mean,
                       MClow = mcbw.j$MC$simpleCI[1],
                       MChigh = mcbw.j$MC$simpleCI[2],
                       nclust=clusters[j],
                       originseason="breed",
                       targetseason="winter") %>% 
    rbind(mc.df)
  
  #Winter:breed
  if(class(siteswbb.j)[1]=="sf"){
    set.seed(1234)
    transwb.j <- estTransition(originSites = siteswbw.j, 
                               targetSites = siteswbb.j, 
                               originPoints = ptwbw.j,
                               targetPoints = ptwbb.j,
                               originAssignment = as.numeric(as.factor(idwb.j$winter_id)),
                               targetAssignment = as.numeric(as.factor(idwb.j$breed_id)),
                               originNames = sort(unique(idwb.j$winter_id)),
                               targetNames = sort(unique(idwb.j$breed_id)),
                               nSamples = 1000,
                               isTelemetry = telbw.j)
    
    mcwb.j <- estStrength(originDist = distwbw.j,
                          targetDist = distwbb.j,
                          originRelAbund = abunw.j,
                          psi = transwb.j,
                          sampleSize = nrow(idwb.j),
                          originNames = sort(unique(idwb.j$winter_id)),
                          targetNames = sort(unique(idwb.j$breed_id)),
                          nSamples = 1000)
  }
  
  mc.df <-  data.frame(MC = mcwb.j$MC$mean,
                       MClow = mcwb.j$MC$simpleCI[1],
                       MChigh = mcwb.j$MC$simpleCI[2],
                       nclust=clusters[j],
                       originseason="winter",
                       targetseason="breed") %>% 
    rbind(mc.df)
  
  #Breed:Fall
  if(class(sitesbfb.j)[1]=="sf"){
    set.seed(1234)
    transbf.j <- estTransition(originSites = sitesbfb.j, 
                               targetSites = sitesbff.j, 
                               originPoints = ptbfb.j,
                               targetPoints = ptbff.j,
                               originAssignment = as.numeric(as.factor(idbf.j$breed_id)),
                               targetAssignment = as.numeric(as.factor(idbf.j$fallmig_id)),
                               originNames = sort(unique(idbf.j$breed_id)),
                               targetNames = sort(unique(idbf.j$fallmig_id)),
                               nSamples = 1000,
                               isTelemetry = telbf.j)
    
    mcbf.j <- estStrength(originDist = distbfb.j,
                          targetDist = distbff.j,
                          originRelAbund = abunb.j,
                          psi = transbf.j,
                          sampleSize = nrow(idbf.j),
                          originNames = sort(unique(idbf.j$breed_id)),
                          targetNames = sort(unique(idbf.j$fallmig_id)),
                          nSamples = 1000)
  }
  
  mc.df <-  data.frame(MC = mcbf.j$MC$mean,
                       MClow = mcbf.j$MC$simpleCI[1],
                       MChigh = mcbf.j$MC$simpleCI[2],
                       nclust=clusters[j],
                       originseason="breed",
                       targetseason="fallmig") %>% 
    rbind(mc.df)
  
  #Winter:Fall
  if(class(siteswfw.j)[1]=="sf"){
    set.seed(1234)
    transwf.j <- estTransition(originSites = siteswfw.j, 
                               targetSites = siteswff.j, 
                               originPoints = ptwfw.j,
                               targetPoints = ptwff.j,
                               originAssignment = as.numeric(as.factor(idwf.j$winter_id)),
                               targetAssignment = as.numeric(as.factor(idwf.j$fallmig_id)),
                               originNames = sort(unique(idwf.j$winter_id)),
                               targetNames = sort(unique(idwf.j$fallmig_id)),
                               nSamples = 1000,
                               isTelemetry = telwf.j)
    
    mcwf.j <- estStrength(originDist = distwfw.j,
                          targetDist = distwff.j,
                          originRelAbund = abunw.j,
                          psi = transwf.j,
                          sampleSize = nrow(idwf.j),
                          originNames = sort(unique(idwf.j$winter_id)),
                          targetNames = sort(unique(idwf.j$fallmig_id)),
                          nSamples = 1000)
  }
  
  mc.df <-  data.frame(MC = mcwf.j$MC$mean,
                       MClow = mcwf.j$MC$simpleCI[1],
                       MChigh = mcwf.j$MC$simpleCI[2],
                       nclust=clusters[j],
                       originseason="winter",
                       targetseason="fallmig") %>% 
    rbind(mc.df)
  
  #Breed:Spring
  if(class(sitesbsb.j)[1]=="sf"){
    set.seed(1234)
    transbs.j <- estTransition(originSites = sitesbsb.j, 
                               targetSites = sitesbss.j, 
                               originPoints = ptbsb.j,
                               targetPoints = ptbss.j,
                               originAssignment = as.numeric(as.factor(idbs.j$breed_id)),
                               targetAssignment = as.numeric(as.factor(idbs.j$springmig_id)),
                               originNames = sort(unique(idbs.j$breed_id)),
                               targetNames = sort(unique(idbs.j$springmig_id)),
                               nSamples = 1000,
                               isTelemetry = telbs.j)
    
    mcbs.j <- estStrength(originDist = distbsb.j,
                          targetDist = distbss.j,
                          originRelAbund = abunb.j,
                          psi = transbs.j,
                          sampleSize = nrow(idbs.j),
                          originNames = sort(unique(idbs.j$breed_id)),
                          targetNames = sort(unique(idbs.j$springmig_id)),
                          nSamples = 1000)
  }
  
  mc.df <-  data.frame(MC = mcbs.j$MC$mean,
                       MClow = mcbs.j$MC$simpleCI[1],
                       MChigh = mcbs.j$MC$simpleCI[2],
                       nclust=clusters[j],
                       originseason="breed",
                       targetseason="springmig") %>% 
    rbind(mc.df)
  
  #Winter:Spring
  if(class(siteswsw.j)[1]=="sf"){
    set.seed(1234)
    transws.j <- estTransition(originSites = siteswsw.j, 
                               targetSites = siteswss.j, 
                               originPoints = ptwsw.j,
                               targetPoints = ptwss.j,
                               originAssignment = as.numeric(as.factor(idws.j$winter_id)),
                               targetAssignment = as.numeric(as.factor(idws.j$springmig_id)),
                               originNames = sort(unique(idws.j$winter_id)),
                               targetNames = sort(unique(idws.j$springmig_id)),
                               nSamples = 1000,
                               isTelemetry = telws.j)
    
    mcws.j <- estStrength(originDist = distwsw.j,
                          targetDist = distwss.j,
                          originRelAbund = abunw.j,
                          psi = transws.j,
                          sampleSize = nrow(idws.j),
                          originNames = sort(unique(idws.j$breed_id)),
                          targetNames = sort(unique(idws.j$springmig_id)),
                          nSamples = 1000)
  }
  
  mc.df <-  data.frame(MC = mcws.j$MC$mean,
                       MClow = mcws.j$MC$simpleCI[1],
                       MChigh = mcws.j$MC$simpleCI[2],
                       nclust=clusters[j],
                       originseason="winter",
                       targetseason="springmig") %>% 
    rbind(mc.df)
  
  print(paste0("Finished ", j, " of ", length(clusters)))
  
}

write.csv(mc.df, "Data/LBCUMigConnectivity.csv", row.names = FALSE)
mc.df <- read.csv("Data/LBCUMigConnectivity.csv")

#13. MC per season----
ggplot(mc.df) +
  geom_point(aes(x=nclust, y=MC)) +
  geom_errorbar(aes(x=nclust, ymin = MClow, ymax = MChigh)) +
  facet_grid(originseason ~ targetseason, scales="free")

#14. MC winner----
mc.sum <- mc.df %>% 
  group_by(nclust) %>% 
  summarize(MCmean = mean(MC),
            MClowq = quantile(MC, 0.025),
            MChighq = quantile(MC, 0.975),
            MClow = mean(MClow),
            MChigh = mean(MChigh)) %>% 
  ungroup()
mc.sum

ggplot(mc.sum) +
  geom_point(aes(x=nclust, y=MCmean)) +
  geom_errorbar(aes(x=nclust, ymin = MClow, ymax = MChigh))

ggplot(mc.sum) +
  geom_point(aes(x=nclust, y=MCmean)) +
  geom_errorbar(aes(x=nclust, ymin = MClowq, ymax = MChighq))

mc.sum %>% 
  arrange(desc(MCmean))
