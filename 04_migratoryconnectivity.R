#https://bookdown.org/adam_smith2/bbsbayes_intro_workshop/

library(tidyverse)
library(data.table)
library(vegan)
library(adehabitatHR)
library(sf)
library(sp)
library(raster)
library(MigConnectivity)

#TO DO: REPLACE RELATIVE ABUNDANCE WITH BBS RESULTS####
#TO DO: THINK ABOUT WHETHER TO ALSO DO STOPOVERS####

#1. Load clusters for BBS routes with LBCU on them----
dat <- read.csv("Data/LBCUKDEClusters.csv")

#2. Read in bbsBayes results for relative abundance information----
bbs <- read.csv("Data/LBCUClusterTrends.csv")

#2. Set up loop through # of clusters---
clusters <- c(2:5,6,8:9)

mantel.list <- list()
mc.list <- list()
for(j in 1:length(clusters)){
  
  #3. Wrangle data----
  breed.j <- dat %>% 
    dplyr::filter(season=="breed",
                  nclust==clusters[j])
  
  winter.j <- dat %>% 
    dplyr::filter(season=="winter",
                  nclust==clusters[j])
  
  fall.j <- dat %>% 
    dplyr::filter(season=="fallmig",
                  nclust==clusters[j])
  
  fall.j <- dat %>% 
    dplyr::filter(season=="springmig",
                  nclust==clusters[j])
  
  bw.j <- rbind(breed.j, winter.j) %>% 
    dplyr::select(id, kdecluster, X, Y, season) %>% 
    pivot_wider(id_cols=id:kdecluster, names_from=season, values_from=X:Y) %>% 
    dplyr::filter(!is.na(X_winter) & !is.na(X_breed)) %>% 
    rename(bird=id)
  
  bf.j <- rbind(breed.j, fall.j) %>% 
    dplyr::select(id, kdecluster, X, Y, season) %>% 
    pivot_wider(id_cols=id:kdecluster, names_from=season, values_from=X:Y) %>% 
    dplyr::filter(!is.na(X_fall) & !is.na(X_breed)) %>% 
    rename(bird=id)
  
  bs.j <- rbind(breed.j, spring.j) %>% 
    dplyr::select(id, kdecluster, X, Y, season) %>% 
    pivot_wider(id_cols=id:kdecluster, names_from=season, values_from=X:Y) %>% 
    dplyr::filter(!is.na(X_spring) & !is.na(X_breed)) %>% 
    rename(bird=id)
  
  #4. Calculate mantel within regions----
  k.list <- list()
  for(k in 1:clusters[j]){
    
    bw.k <- bw.j %>% 
      dplyr::filter(kdecluster==k)

    b.k <- bw.k %>% 
      dplyr::select(X_breed, Y_breed) %>% 
      vegdist("euclidean")
    
    w.k <- bw.k %>% 
      dplyr::select(X_winter, Y_winter) %>% 
      vegdist("euclidean")
    
    mantel.k <- mantel(b.k, w.k)
    
    k.list[[k]] <- data.frame(r = mantel.k[["statistic"]],
                              p = mantel.k[["signif"]],
                              n=nrow(bw.k),
                              kdecluster = k,
                              nclust=clusters[j])
    
  }
  mantel.list[[j]] <- rbindlist(k.list)
  
  #5. Set up target grid for MC estimation----
  pts.utm <- st_as_sf(bw.j, coords=c("X_winter", "Y_winter"), crs=3857) %>% 
    as_Spatial()
  
  bb <- bbox(pts.utm)
  cs <- c(100000, 100000)  # cell size 
  cc <- bb[, 1] + (cs/2)  # cell offset
  cd <- ceiling(diff(t(bb))/cs)  # number of cells per direction
  grd <- GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd)
  spgrd <- SpatialGridDataFrame(grd,
                                data=data.frame(id=1:prod(cd)),
                                proj4string=CRS(proj4string(pts.utm)))
  
  plot(spgrd)
  plot(pts.utm, add=TRUE)
  
  id.j <- bw.j %>% 
    rename(breed_id=kdecluster) %>% 
    cbind(id=over(pts.utm, spgrd)) %>% 
    dplyr::rename(winter_id=id)
  
  #6. Number of regions----
  nb <- length(unique(id.j$breed_id))
  nw <- length(unique(id.j$winter_id))
  
  #7. Distance matrices between centroids----
  distb.j <- id.j %>% 
    group_by(breed_id) %>% 
    summarize(X=mean(X_breed),
              Y=mean(Y_breed)) %>% 
    ungroup() %>% 
    arrange(breed_id) %>% 
    dplyr::select(X, Y) %>% 
    data.frame() %>% 
    distFromPos("plane")
  dimnames(distb.j) <- list(sort(unique(id.j$breed_id)), sort(unique(id.j$breed_id)))
  
  distw.j <- id.j %>% 
    group_by(winter_id) %>% 
    summarize(X=mean(X_winter),
              Y=mean(Y_winter)) %>% 
    ungroup() %>% 
    arrange(winter_id) %>% 
    dplyr::select(X, Y) %>% 
    data.frame() %>% 
    distFromPos("plane")
  dimnames(distw.j) <- list(sort(unique(id.j$winter_id)), sort(unique(id.j$winter_id)))
  
  #8. Create sp objects for individual locations----
  ptb.j <- id.j %>% 
    st_as_sf(coords=c("X_breed", "Y_breed"), crs=3857) %>% 
    dplyr::select(geometry)
  
  ptw.j <- id.j %>% 
    st_as_sf(coords=c("X_winter", "Y_winter"), crs=3857) %>% 
    dplyr::select(geometry)
  
  #9. Create breeding region polygons----
  sp.j <- SpatialPointsDataFrame(coords=cbind(id.j$X_breed, id.j$Y_breed), 
                                 data=data.frame(ID=id.j$breed_id),
                                 proj4string = CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"))
  
  sitesb.j <- mcp(sp.j[,1], percent=30) %>% 
    st_as_sf() %>% 
    dplyr::select(geometry)
  plot(sitesb.j)

  #10. Create winter region polygons----
  sitesw.j <- spgrd %>% 
    raster() %>% 
    rasterToPolygons() %>% 
    st_as_sf() %>% 
    dplyr::filter(id %in% unique(id.j$winter_id)) %>% 
    dplyr::select(geometry)
  
  #11. Add relative abundance----
  abun.j <- rep(1/nb, nb)
  
  #12. Define tag type----
  gl.j <- rep(FALSE, nrow(id.j))
  
  #13. Define error----
  LongError <- rnorm(100, 20, 5)
  LatError <- rnorm(100, 20, 5)
  geo.error.model <- lm(cbind(LongError,LatError) ~ 1) 
  geo.bias <- coef(geo.error.model)
  geo.vcov <- vcov(geo.error.model)
  
  #14. Estimate MC----
  set.seed(1234)
  mc.j <-try(estMC(originDist = distb.j, 
               targetDist = distw.j, 
               originSites = sitesb.j, 
               targetSites = sitesw.j, 
               originPoints = ptb.j,
               targetPoints = ptw.j,
               nSamples = 100,
               isGL = gl.j,
               geoBias = geo.bias,
               geoVCov = geo.vcov,
               originRelAbund = abun.j,
               verbose = 1))
  
  if (class(mc.j)[1]=="estMC"){

    mc.list[[j]] <-  data.frame(MC = mc.j$MC$mean,
                                MClow = mc.j$MC$simpleCI[1],
                                MChigh = mc.j$MC$simpleCI[2],
                                nclust=clusters[j])
  }

}

#15. Collapse & summarize results----
mantel <- rbindlist(mantel.list)
mantel.sum <- mantel %>% 
  group_by(nclust) %>% 
  summarize(mean=mean(r),
            sd=sd(r))

mc <- rbindlist(mc.list)

results <- full_join(mc, mantel)
sum <- full_join(mc, mantel.sum)

ggplot(sum) +
  geom_point(aes(x=MC, y=mean, colour=factor(nclust)), size=5)

#16. Save results----
write.csv(results, "LBCUMigConnectivity.csv")
