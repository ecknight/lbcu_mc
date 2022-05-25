library(tidyverse)
library(data.table)
library(vegan)
library(adehabitatHR)
library(sf)
library(sp)
library(raster)
library(MigConnectivity)
library(ebirdst)

#1. Load clusters for BBS routes with LBCU on them----
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
  group_by(nclust, kdecluster, season) %>% 
  summarize(breeding = mean(breeding, na.rm=TRUE),
            postbreeding_migration = mean(postbreeding_migration, na.rm=TRUE),
            nonbreeding = mean(nonbreeding, na.rm=TRUE),
            prebreeding_migration = mean(prebreeding_migration, na.rm=TRUE),
            n=n()) %>% 
  ungroup()

#3. Set up loop through # of clusters---
clusters <- c(2:9)

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
  nw.j <- idw.j %>% 
    group_by(breed_id) %>% 
    summarize(n=n()) %>% 
    ungroup() %>% 
    summarize(n = min(n))
  
  if(nw.j$n <= 4){
    siteswb.j <- idw.j %>% 
      group_by(breed_id) %>% 
      summarize(X=mean(X_breed),
                Y=mean(Y_breed)) %>% 
      ungroup() %>% 
      st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
      st_buffer(40000) %>% 
      st_cast("MULTIPOLYGON")
  }
  
  if(nw.j$n > 4){
    spw.j <- SpatialPointsDataFrame(coords=cbind(idw.j$X_breed, idw.j$Y_breed), 
                                    data=data.frame(ID=idw.j$breed_id),
                                    proj4string = CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"))
    
    siteswb.j <- try(mcp(spw.j[,1], percent=30) %>% 
                       st_as_sf() %>% 
                       dplyr::select(geometry)) %>% 
      st_cast("MULTIPOLYGON")
  }

  #Fall
  nf.j <- idf.j %>% 
    group_by(breed_id) %>% 
    summarize(n=n()) %>% 
    ungroup() %>% 
    summarize(n = min(n))
  
  if(nf.j$n <= 4){
    sitesfb.j <- idf.j %>% 
      group_by(breed_id) %>% 
      summarize(X=mean(X_breed),
                Y=mean(Y_breed)) %>% 
      ungroup() %>% 
      st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
      st_buffer(40000) %>% 
      st_cast("MULTIPOLYGON")
  }
  
  if(nf.j$n > 4){
    spf.j <- SpatialPointsDataFrame(coords=cbind(idf.j$X_breed, idf.j$Y_breed), 
                                    data=data.frame(ID=idf.j$breed_id),
                                    proj4string = CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"))
    
    sitesfb.j <- try(mcp(spf.j[,1], percent=30) %>% 
                       st_as_sf() %>% 
                       dplyr::select(geometry))
  }

  #spring
  ns.j <- ids.j %>% 
    group_by(breed_id) %>% 
    summarize(n=n()) %>% 
    ungroup() %>% 
    summarize(n = min(n))
  
  if(ns.j$n <= 4){
    sitessb.j <- ids.j %>% 
      group_by(breed_id) %>% 
      summarize(X=mean(X_breed),
                Y=mean(Y_breed)) %>% 
      ungroup() %>% 
      st_as_sf(coords=c("X", "Y"), crs=3857) %>% 
      st_buffer(40000) %>% 
      st_cast("MULTIPOLYGON")
  }
  
  if(ns.j$n > 4){
    sps.j <- SpatialPointsDataFrame(coords=cbind(ids.j$X_breed, ids.j$Y_breed), 
                                    data=data.frame(ID=ids.j$breed_id),
                                    proj4string = CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"))
    
    sitessb.j <- try(mcp(sps.j[,1], percent=30) %>% 
                       st_as_sf() %>% 
                       dplyr::select(geometry))
  }


  #11. Create winter region polygons----
  
  #Winter
  sitesww.j <- spgrdw %>% 
    raster() %>% 
    rasterToPolygons() %>% 
    st_as_sf() %>% 
    dplyr::filter(id %in% unique(idw.j$winter_id)) %>% 
    dplyr::select(geometry) %>% 
    st_cast("MULTIPOLYGON")
  
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
  bbs.j <- abun %>% 
    dplyr::filter(nclust==clusters[j], season=="breed") %>% 
    mutate(abun = breeding/(sum(breeding)))
  abun.j <- bbs.j$abun
  
  #13. Define tag type----
  telw.j <- rep(TRUE, nrow(idw.j))
  telf.j <- rep(TRUE, nrow(idf.j))
  tels.j <- rep(TRUE, nrow(ids.j))
  
  #15. Estimate MC----
  
  #Winter
  if(class(siteswb.j)[1]=="sf"){
    set.seed(1234)
  transw.j <- estTransition(originSites = siteswb.j, 
                            targetSites = sitesww.j, 
                            originPoints = ptwb.j,
                            targetPoints = ptww.j,
                            originAssignment = as.numeric(as.factor(idw.j$breed_id)),
                            targetAssignment = as.numeric(as.factor(idw.j$winter_id)),
                            originNames = sort(unique(idw.j$breed_id)),
                            targetNames = sort(unique(idw.j$winter_id)),
                            nSamples = 100,
                            isTelemetry = telw.j)
  
  mcw.j <- estStrength(originDist = distwb.j,
                      targetDist = distww.j,
                      originRelAbund = abun.j,
                      psi = transw.j,
                      sampleSize = nrow(idw.j),
                      originNames = sort(unique(idw.j$breed_id)),
                      targetNames = sort(unique(idw.j$winter_id)),
                      nSamples = 100)
  }
  
  mc.df <-  data.frame(MC = mcw.j$MC$mean,
                       MClow = mcw.j$MC$simpleCI[1],
                       MChigh = mcw.j$MC$simpleCI[2],
                       nclust=clusters[j],
                       season="winter") %>% 
    rbind(mc.df)
  
  #Fall
  if(class(sitesfb.j)[1]=="sf"){
    set.seed(1234)
    transf.j <- estTransition(originSites = sitesfb.j, 
                              targetSites = sitesff.j, 
                              originPoints = ptfb.j,
                              targetPoints = ptff.j,
                              originAssignment = as.numeric(as.factor(idf.j$breed_id)),
                              targetAssignment = as.numeric(as.factor(idf.j$fallmig_id)),
                              originNames = sort(unique(idf.j$breed_id)),
                              targetNames = sort(unique(idf.j$fallmig_id)),
                              nSamples = 100,
                              isTelemetry = telf.j)
    
    mcf.j <- estStrength(originDist = distfb.j,
                         targetDist = distff.j,
                         originRelAbund = abun.j,
                         psi = transf.j,
                         sampleSize = nrow(idf.j),
                         originNames = sort(unique(idf.j$breed_id)),
                         targetNames = sort(unique(idf.j$fallmig_id)),
                         nSamples = 100)
  }
  
  mc.df <-  data.frame(MC = mcf.j$MC$mean,
                       MClow = mcf.j$MC$simpleCI[1],
                       MChigh = mcf.j$MC$simpleCI[2],
                       nclust=clusters[j],
                       season="fallmig") %>% 
    rbind(mc.df)
  
  #Spring
  if(class(sitessb.j)[1]=="sf"){
    set.seed(1234)
    transs.j <- estTransition(originSites = sitessb.j, 
                              targetSites = sitesss.j, 
                              originPoints = ptsb.j,
                              targetPoints = ptss.j,
                              originAssignment = as.numeric(as.factor(ids.j$breed_id)),
                              targetAssignment = as.numeric(as.factor(ids.j$springmig_id)),
                              originNames = sort(unique(ids.j$breed_id)),
                              targetNames = sort(unique(ids.j$springmig_id)),
                              nSamples = 1000,
                              isTelemetry = tels.j)
    
    mcs.j <- estStrength(originDist = distsb.j,
                         targetDist = distss.j,
                         originRelAbund = abun.j,
                         psi = transs.j,
                         sampleSize = nrow(ids.j),
                         originNames = sort(unique(ids.j$breed_id)),
                         targetNames = sort(unique(ids.j$springmig_id)),
                         nSamples = 1000)
  }
  
  mc.df <-  data.frame(MC = mcs.j$MC$mean,
                       MClow = mcs.j$MC$simpleCI[1],
                       MChigh = mcs.j$MC$simpleCI[2],
                       nclust=clusters[j],
                       season="springmig") %>% 
    rbind(mc.df)

}

#16. Collapse & summarize results----
mantel <- rbindlist(mantel.list)
mantel.sum <- mantel %>% 
  group_by(nclust, season) %>% 
  summarize(mean=mean(r, na.rm=TRUE),
            sd=sd(r))

results <- full_join(mc.df, mantel)
sum <- full_join(mc.df, mantel.sum) %>% 
  mutate(mc.s = (MC - min(MC, na.rm=TRUE))/(max(MC, na.rm=TRUE) - min(MC, na.rm=TRUE)),
         r.s = 1-(mean - min(mean, na.rm=TRUE))/(max(mean, na.rm=TRUE) - min(mean, na.rm=TRUE)),
         f = (mc.s*r.s)/(mc.s+r.s))

ggplot(sum) +
  geom_line(aes(x=nclust, y=mc.s), colour="blue") +
  geom_line(aes(x=nclust, y=r.s), colour="red") +
  geom_line(aes(x=nclust, y=f)) +
  facet_wrap(~season)

sum.sum <- sum %>% 
  group_by(nclust) %>% 
  summarize(mc.s= mean(mc.s),
            r.s = mean(r.s),
            f = (mc.s*r.s)/(mc.s+r.s))

ggplot(sum.sum) +
  geom_line(aes(x=nclust, y=r.s), colour="red") +
  geom_line(aes(x=nclust, y=mc.s), colour="blue") +
  geom_line(aes(x=nclust, y=f))

#17. Save results----
write.csv(results, "Data/LBCUMigConnectivity.csv")
