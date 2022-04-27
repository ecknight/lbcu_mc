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

#3. Set up loop through # of clusters---
clusters <- c(2:6,8:9)

mantel.list <- list()
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
    dplyr::select(id, kdecluster, X, Y, season) %>% 
    pivot_wider(id_cols=id:kdecluster, names_from=season, values_from=X:Y) %>% 
    dplyr::filter(!is.na(X_winter) & !is.na(X_breed)) %>% 
    rename(bird=id)
  
  bf.j <- rbind(breed.j, fall.j) %>% 
    dplyr::select(id, kdecluster, X, Y, season) %>% 
    pivot_wider(id_cols=id:kdecluster, names_from=season, values_from=X:Y) %>% 
    dplyr::filter(!is.na(X_fallmig) & !is.na(X_breed)) %>% 
    rename(bird=id)
  
  bs.j <- rbind(breed.j, spring.j) %>% 
    dplyr::select(id, kdecluster, X, Y, season) %>% 
    pivot_wider(id_cols=id:kdecluster, names_from=season, values_from=X:Y) %>% 
    dplyr::filter(!is.na(X_springmig) & !is.na(X_breed)) %>% 
    rename(bird=id)
  
  #5. Calculate mantel within regions----
  mantel.df <- data.frame()
  for(k in 1:clusters[j]){
    
    #Winter
    bw.k <- bw.j %>% 
      dplyr::filter(kdecluster==k)

    bwb.k <- bw.k %>% 
      dplyr::select(X_breed, Y_breed) %>% 
      vegdist("euclidean")
    
    bww.k <- bw.k %>% 
      dplyr::select(X_winter, Y_winter) %>% 
      vegdist("euclidean")

    mantelbw.k <- try(mantel(bwb.k, bww.k))
    
    if(class(mantelbw.k)=="mantel"){
      mantel.df <- data.frame(r = mantelbw.k[["statistic"]],
                                p = mantelbw.k[["signif"]],
                                n=nrow(bw.k),
                                kdecluster = k,
                                nclust=clusters[j],
                                season="winter") %>% 
        rbind(mantel.df)
    }

    
    #Fallmig
    bf.k <- bf.j %>% 
      dplyr::filter(kdecluster==k)
    
    bfb.k <- bf.k %>% 
      dplyr::select(X_breed, Y_breed) %>% 
      vegdist("euclidean")
    
    bff.k <- bf.k %>% 
      dplyr::select(X_fallmig, Y_fallmig) %>% 
      vegdist("euclidean")
    
    mantelbf.k <- try(mantel(bfb.k, bff.k))

    if(class(mantelbf.k)=="mantel"){
      mantel.df <- data.frame(r = mantelbf.k[["statistic"]],
                              p = mantelbf.k[["signif"]],
                              n=nrow(bf.k),
                              kdecluster = k,
                              nclust=clusters[j],
                              season="fallmig") %>% 
        rbind(mantel.df)
    }
    
    #Springmig
    bs.k <- bs.j %>% 
      dplyr::filter(kdecluster==k)
    
    bsb.k <- bs.k %>% 
      dplyr::select(X_breed, Y_breed) %>% 
      vegdist("euclidean")
    
    bss.k <- bs.k %>% 
      dplyr::select(X_springmig, Y_springmig) %>% 
      vegdist("euclidean")
    
    mantelbs.k <- try(mantel(bsb.k, bss.k))
    
    if(class(mantelbs.k)=="mantel"){
      mantel.df <- data.frame(r = mantelbs.k[["statistic"]],
                              p = mantelbs.k[["signif"]],
                              n=nrow(bs.k),
                              kdecluster = k,
                              nclust=clusters[j],
                              season="springmig") %>% 
        rbind(mantel.df)
    }
  }
  mantel.list[[j]] <- mantel.df
  
  #6. Set up target grids for MC estimation----
  
  #Winter
  ptsw <- st_as_sf(bw.j, coords=c("X_winter", "Y_winter"), crs=3857) %>% 
    as_Spatial()
  
  bb <- bbox(ptsw)
  cs <- c(500000, 500000)  # cell size 
  cc <- bb[, 1] + (cs/2)  # cell offset
  cd <- ceiling(diff(t(bb))/cs)  # number of cells per direction
  grd <- GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd)
  spgrdw <- SpatialGridDataFrame(grd,
                                data=data.frame(id=1:prod(cd)),
                                proj4string=CRS(proj4string(ptsw)))
  
  idw.j <- bw.j %>% 
    rename(breed_id=kdecluster) %>% 
    cbind(id=over(ptsw, spgrdw)) %>% 
    dplyr::rename(winter_id=id)
  
  #Fall
  ptsf <- st_as_sf(bf.j, coords=c("X_fallmig", "Y_fallmig"), crs=3857) %>% 
    as_Spatial()
  
  bb <- bbox(ptsf)
  cs <- c(500000, 500000)  # cell size 
  cc <- bb[, 1] + (cs/2)  # cell offset
  cd <- ceiling(diff(t(bb))/cs)  # number of cells per direction
  grd <- GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd)
  spgrdf <- SpatialGridDataFrame(grd,
                                data=data.frame(id=1:prod(cd)),
                                proj4string=CRS(proj4string(ptsf)))

  idf.j <- bf.j %>% 
    rename(breed_id=kdecluster) %>% 
    cbind(id=over(ptsf, spgrdf)) %>% 
    dplyr::rename(fallmig_id=id)
  
  #Spring
  ptss <- st_as_sf(bs.j, coords=c("X_springmig", "Y_springmig"), crs=3857) %>% 
    as_Spatial()
  
  bb <- bbox(ptss)
  cs <- c(500000, 500000)  # cell size 
  cc <- bb[, 1] + (cs/2)  # cell offset
  cd <- ceiling(diff(t(bb))/cs)  # number of cells per direction
  grd <- GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd)
  spgrds <- SpatialGridDataFrame(grd,
                                data=data.frame(id=1:prod(cd)),
                                proj4string=CRS(proj4string(ptss)))
  
  ids.j <- bs.j %>% 
    rename(breed_id=kdecluster) %>% 
    cbind(id=over(ptss, spgrds)) %>% 
    dplyr::rename(springmig_id=id)
  
  #7. Number of regions----
  nwb <- length(unique(idw.j$breed_id))
  nww <- length(unique(idw.j$winter_id))
  
  nfb <- length(unique(idf.j$breed_id))
  nff <- length(unique(idf.j$fallmig_id))
  
  nsb <- length(unique(ids.j$breed_id))
  nss <- length(unique(ids.j$springmig_id))
  
  #8. Distance matrices between centroids----
  
  #winter
  distwb.j <- idw.j %>% 
    group_by(breed_id) %>% 
    summarize(X=mean(X_breed),
              Y=mean(Y_breed)) %>% 
    ungroup() %>% 
    arrange(breed_id) %>% 
    dplyr::select(X, Y) %>% 
    data.frame() %>% 
    distFromPos("plane")
  dimnames(distwb.j) <- list(sort(unique(idw.j$breed_id)), sort(unique(idw.j$breed_id)))
  
  distww.j <- idw.j %>% 
    group_by(winter_id) %>% 
    summarize(X=mean(X_winter),
              Y=mean(Y_winter)) %>% 
    ungroup() %>% 
    arrange(winter_id) %>% 
    dplyr::select(X, Y) %>% 
    data.frame() %>% 
    distFromPos("plane")
  dimnames(distww.j) <- list(sort(unique(idw.j$winter_id)), sort(unique(idw.j$winter_id)))
  
  #fall
  distfb.j <- idf.j %>% 
    group_by(breed_id) %>% 
    summarize(X=mean(X_breed),
              Y=mean(Y_breed)) %>% 
    ungroup() %>% 
    arrange(breed_id) %>% 
    dplyr::select(X, Y) %>% 
    data.frame() %>% 
    distFromPos("plane")
  dimnames(distfb.j) <- list(sort(unique(idf.j$breed_id)), sort(unique(idf.j$breed_id)))
  
  distff.j <- idf.j %>% 
    group_by(fallmig_id) %>% 
    summarize(X=mean(X_fallmig),
              Y=mean(Y_fallmig)) %>% 
    ungroup() %>% 
    arrange(fallmig_id) %>% 
    dplyr::select(X, Y) %>% 
    data.frame() %>% 
    distFromPos("plane")
  dimnames(distff.j) <- list(sort(unique(idf.j$fallmig_id)), sort(unique(idf.j$fallmig_id)))
  
  #spring
  distsb.j <- ids.j %>% 
    group_by(breed_id) %>% 
    summarize(X=mean(X_breed),
              Y=mean(Y_breed)) %>% 
    ungroup() %>% 
    arrange(breed_id) %>% 
    dplyr::select(X, Y) %>% 
    data.frame() %>% 
    distFromPos("plane")
  dimnames(distsb.j) <- list(sort(unique(ids.j$breed_id)), sort(unique(ids.j$breed_id)))
  
  distss.j <- ids.j %>% 
    group_by(springmig_id) %>% 
    summarize(X=mean(X_springmig),
              Y=mean(Y_springmig)) %>% 
    ungroup() %>% 
    arrange(springmig_id) %>% 
    dplyr::select(X, Y) %>% 
    data.frame() %>% 
    distFromPos("plane")
  dimnames(distss.j) <- list(sort(unique(ids.j$springmig_id)), sort(unique(ids.j$springmig_id)))
  
  #9. Create point objects for individual locations----
  
  #Winter
  ptwb.j <- idw.j %>% 
    st_as_sf(coords=c("X_breed", "Y_breed"), crs=3857) %>% 
    dplyr::select(geometry)
  
  ptww.j <- idw.j %>% 
    st_as_sf(coords=c("X_winter", "Y_winter"), crs=3857) %>% 
    dplyr::select(geometry)
  
  #fall
  ptfb.j <- idf.j %>% 
    st_as_sf(coords=c("X_breed", "Y_breed"), crs=3857) %>% 
    dplyr::select(geometry)
  
  ptff.j <- idf.j %>% 
    st_as_sf(coords=c("X_fallmig", "Y_fallmig"), crs=3857) %>% 
    dplyr::select(geometry)
  
  #Winter
  ptsb.j <- ids.j %>% 
    st_as_sf(coords=c("X_breed", "Y_breed"), crs=3857) %>% 
    dplyr::select(geometry)
  
  ptss.j <- ids.j %>% 
    st_as_sf(coords=c("X_springmig", "Y_springmig"), crs=3857) %>% 
    dplyr::select(geometry)
  
  #10. Create breeding region polygons----
  
  #Winter
  spw.j <- SpatialPointsDataFrame(coords=cbind(idw.j$X_breed, idw.j$Y_breed), 
                                 data=data.frame(ID=idw.j$breed_id),
                                 proj4string = CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"))
  
  siteswb.j <- try(mcp(spw.j[,1], percent=30) %>% 
    st_as_sf() %>% 
    dplyr::select(geometry))
  
  #Fall
  spf.j <- SpatialPointsDataFrame(coords=cbind(idf.j$X_breed, idf.j$Y_breed), 
                                  data=data.frame(ID=idf.j$breed_id),
                                  proj4string = CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"))
  
  sitesfb.j <- try(mcp(spf.j[,1], percent=30) %>% 
    st_as_sf() %>% 
    dplyr::select(geometry))
  
  #spring
  sps.j <- SpatialPointsDataFrame(coords=cbind(ids.j$X_breed, ids.j$Y_breed), 
                                  data=data.frame(ID=ids.j$breed_id),
                                  proj4string = CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"))
  
  sitessb.j <- try(mcp(sps.j[,1], percent=30) %>% 
    st_as_sf() %>% 
    dplyr::select(geometry))

  #11. Create winter region polygons----
  
  #Winter
  sitesww.j <- spgrdw %>% 
    raster() %>% 
    rasterToPolygons() %>% 
    st_as_sf() %>% 
    dplyr::filter(id %in% unique(idw.j$winter_id)) %>% 
    dplyr::select(geometry)
  
  #Fall
  sitesff.j <- spgrdf %>% 
    raster() %>% 
    rasterToPolygons() %>% 
    st_as_sf() %>% 
    dplyr::filter(id %in% unique(idf.j$fallmig_id)) %>% 
    dplyr::select(geometry)
  
  #Spring
  sitesss.j <- spgrds %>% 
    raster() %>% 
    rasterToPolygons() %>% 
    st_as_sf() %>% 
    dplyr::filter(id %in% unique(ids.j$springmig_id)) %>% 
    dplyr::select(geometry)
  
  #12. Add relative abundance----
  bbs.j <- bbs %>% 
    dplyr::filter(nclust==clusters[j]) %>% 
    mutate(abun = Relative_Abundance/(sum(Relative_Abundance)))
  abun.j <- bbs.j$abun
  
  #13. Define tag type----
  glw.j <- rep(FALSE, nrow(idw.j))
  glf.j <- rep(FALSE, nrow(idf.j))
  gls.j <- rep(FALSE, nrow(ids.j))
  
  #14. Define error----
  LongError <- rnorm(100, 20, 5)
  LatError <- rnorm(100, 20, 5)
  geo.error.model <- lm(cbind(LongError,LatError) ~ 1) 
  geo.bias <- coef(geo.error.model)
  geo.vcov <- vcov(geo.error.model)
  
  #15. Estimate MC----
  
  #Winter
  if(class(siteswb.j)[1]=="sf"){
    set.seed(1234)
    mcw.j <-try(estMC(originDist = distwb.j, 
                      targetDist = distww.j, 
                      originSites = siteswb.j, 
                      targetSites = sitesww.j, 
                      originPoints = ptwb.j,
                      targetPoints = ptww.j,
                      nSamples = 10,
                      isGL = glw.j,
                      geoBias = geo.bias,
                      geoVCov = geo.vcov,
                      originRelAbund = abun.j,
                      verbose = 1))
    
    if (class(mcw.j)[1]=="estMC"){
      
      mc.df <-  data.frame(MC = mcw.j$MC$mean,
                           MClow = mcw.j$MC$simpleCI[1],
                           MChigh = mcw.j$MC$simpleCI[2],
                           nclust=clusters[j],
                           season="winter") %>% 
        rbind(mc.df)
    }
    
  }

  #Fall
  if(class(sitesfb.j)[1]=="sf"){
    
    mcf.j <-try(estMC(originDist = distfb.j, 
                      targetDist = distff.j, 
                      originSites = sitesfb.j, 
                      targetSites = sitesff.j, 
                      originPoints = ptfb.j,
                      targetPoints = ptff.j,
                      nSamples = 10,
                      isGL = glf.j,
                      geoBias = geo.bias,
                      geoVCov = geo.vcov,
                      originRelAbund = abun.j,
                      verbose = 1))
    
    if (class(mcf.j)[1]=="estMC"){
      
      mc.df <-  data.frame(MC = mcf.j$MC$mean,
                           MClow = mcf.j$MC$simpleCI[1],
                           MChigh = mcf.j$MC$simpleCI[2],
                           nclust=clusters[j],
                           season="fallmig") %>% 
        rbind(mc.df)
    }
    
  }
  
  #Spring
  if(class(sitessb.j)[1]=="sf"){
    
    mcs.j <-try(estMC(originDist = distsb.j, 
                      targetDist = distss.j, 
                      originSites = sitessb.j, 
                      targetSites = sitesss.j, 
                      originPoints = ptsb.j,
                      targetPoints = ptss.j,
                      nSamples = 10,
                      isGL = gls.j,
                      geoBias = geo.bias,
                      geoVCov = geo.vcov,
                      originRelAbund = abun.j,
                      verbose = 1))
    
    if (class(mcs.j)[1]=="estMC"){
      
      mc.df <-  data.frame(MC = mcs.j$MC$mean,
                           MClow = mcs.j$MC$simpleCI[1],
                           MChigh = mcs.j$MC$simpleCI[2],
                           nclust=clusters[j],
                           season="springmig") %>% 
        rbind(mc.df)
    }
  }
}

#16. Collapse & summarize results----
mantel <- rbindlist(mantel.list)
mantel.sum <- mantel %>% 
  group_by(nclust, season) %>% 
  summarize(mean=mean(r),
            sd=sd(r))

results <- full_join(mc.df, mantel)
sum <- full_join(mc.df, mantel.sum)

ggplot(sum) +
  geom_point(aes(x=MC, y=mean, colour=factor(nclust)), size=5) +
  facet_wrap(~season)

ggplot(sum) +
  geom_point(aes(x=nclust, y=mean, colour=season)) +
  geom_smooth(aes(x=nclust, y=mean, colour=season))

ggplot(sum) +
  geom_point(aes(x=nclust, y=MC, colour=season)) +
  geom_line(aes(x=nclust, y=MC, colour=season))

#17. Save results----
write.csv(results, "Data/LBCUMigConnectivity.csv")
