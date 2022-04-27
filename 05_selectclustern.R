library(tidyverse)
library(gtools)

mc <- read.csv("Data/LBCUMigConnectivity.csv")

ggplot(mc) +
  geom_point(aes(x=kdecluster, y=r, colour=season)) +
  facet_wrap(~nclust)

trend <- read.csv("Data/LBCUClusterTrends.csv")

ggplot(trend) +
  geom_point(aes(x=Region_alt, y=Trend)) +
  geom_errorbar(aes(x=Region_alt, ymin=Trend_Q0.025, ymax=Trend_Q0.975)) +
  geom_hline(aes(yintercept=0)) +
  facet_wrap(~nclust)

overlap <- data.frame(nclust=c(2:6, 8:9),
                      diff = c(1,0,3,5,9,3,10),
                      n=c(1,3,6,10,15,10,36)) %>% 
  mutate(ratio = diff/n)

sum <- mc %>% 
  group_by(nclust, season) %>% 
  summarize(MC=mean(MC, na.rm=TRUE),
            r=mean(r, na.rm=TRUE)) %>% 
  left_join(overlap) %>% 
  mutate(rr = 1-r,
         sum=MC+rr+ratio) 

ggplot(sum) +
  geom_point(aes(x=nclust, y=sum, colour=season))
