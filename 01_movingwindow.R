library(tidyverse)
library(sf)
library(data.table)
library(vegan)

options(scipen=9999999)

#1. Import data----
mn <- read.csv("Data/LBCUMCLocations.csv")

#5.Try MC moving window----
#Parameters
width <- c(100000, 500000, 1000000, 1500000, 2000000, 2500000)
overlap <- c(0.25, 0.5, 0.75)

loop <- data.frame(expand.grid(width=width, overlap=overlap)) %>% 
  mutate(sliver = width*(1-overlap))

boundaries <- mn %>% 
  dplyr::filter(season=="breed") %>% 
  summarize(minX = round(min(X), -5),
            maxX = round(max(X), -5),
            minY = round(min(Y), -5),
            maxY = round(max(Y), -5))

mantel.list <- list()
for(i in 1:nrow(loop)){
  
  width.i <- loop$width[i]
  overlap.i <- loop$overlap[i]
  sliver.i <- loop$sliver[i]
  
  #Calculate # of windows needed for window width & overlap
  n.wide <- ceiling(((boundaries$maxX - boundaries$minX) - width.i)/sliver.i)
  n.tall <- ceiling(((boundaries$maxY - boundaries$minY) - width.i)/sliver.i)
  
  j.list <- list()
  for(j in 1:n.wide){
    
    #Calculate beginning of latitude window
    start.j <- boundaries$minX + sliver.i*(j-1)
    end.j <- start.j + width.i
    
    k.list <- list()
    for(k in 1:n.tall){
      start.k <- boundaries$minY + sliver.i*(k-1)
      end.k <- start.k + width.i
      
      #filter data
      breed.k <- mn %>% 
        dplyr::filter(season=="breed",
                      X > start.j & X < end.j,
                      Y > start.k & Y < end.k)
      
      if(nrow(breed.k)<=2) {
        next
      }
      
      else{
        
        breed.mat.k <- breed.k %>% 
          dplyr::select(X, Y) %>% 
          vegdist("euclidean")
        
        winter.k <- mn.utm %>% 
          dplyr::filter(season=="winter",
                        id %in% breed.k$id)
        
        winter.mat.k <- winter.k %>% 
          dplyr::select(X, Y) %>% 
          vegdist("euclidean")
        
        #Calculate mantel
        mantel.k <- mantel(breed.mat.k, winter.mat.k)
        
        k.list[[k]] <- data.frame(r = mantel.k[["statistic"]],
                                  p = mantel.k[["signif"]],
                                  n=nrow(breed.k),
                                  start.x = start.j,
                                  end.x = end.j,
                                  start.y = start.k,
                                  end.y = end.k,
                                  width = width.i,
                                  overlap = overlap.i,
                                  i=i,
                                  j=j,
                                  k=k)
      }
    }
    if(length(k.list)==1){
      j.list[[j]] <- data.frame(k.list)
    }
    if(length(k.list) > 1){
      j.list[[j]] <- rbindlist(k.list)
    }
  }
  mantel.list[[i]] <- rbindlist(j.list)
  
  print(paste0("Finished width ", width.i, " and overlap ", overlap.i," - ", i, " of ", nrow(loop), " iterations"))
  
}

mantel <- rbindlist(mantel.list)

write.csv(mantel, "Data/MovingWindow_XY.csv", row.names = FALSE)

ggplot(mantel) +
  geom_path(data=mn.utm, aes(x=X, y=Y, group=id), colour="grey50") +
  geom_point(data=mn.utm, aes(x=X, y=Y, colour=season)) +
  geom_rect(aes(xmin= start.x, xmax=end.x, ymin=start.y, ymax=end.y, fill=r), alpha=.5) +
  scale_colour_manual(values=c("grey10", "grey90")) +
  facet_grid(width ~ overlap) +
  scale_fill_viridis_c()

ggsave("Figs/MigratoryDivides.jpeg", width = 10, height = 8)
