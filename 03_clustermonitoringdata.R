library(tidyverse)
library(class)
library(sf)
library(data.table)
library(adehabitatHR)
library(rpart)
library(party)
library(klaR)
library(caret)

options(scipen=9999)

#A. PREAMBLE####

#1. Import tracking data with training clusters----
#track.raw <- read.csv("Data/LBCUKDEClusters.csv")
track.raw <- read.csv("Data/LBCUManualClusters.csv") %>% 
  mutate(nclust=4,
         boot=1,
         kdecluster = factor(group))

#2. Import bbs monitoring data----
load("~/Library/Application Support/bbsBayes/bbs_raw_data.RData")

bbs <- bbs_data[["bird"]] %>% 
  dplyr::filter(AOU==2640,
                Year >= 1970)

routes <- bbs %>% 
  dplyr::select(RouteDataID) %>%
  unique() %>% 
  left_join(bbs_data[["route"]]) %>% 
  group_by(countrynum, statenum, Route, Latitude, Longitude) %>%
  summarize(n = n()) %>% 
  dplyr::filter(!is.na(Latitude)) %>% 
  ungroup()

bbs.sf <- routes %>% 
  st_as_sf(coords=c("Longitude", "Latitude"), crs=4326) %>% 
  st_transform(crs=3857) %>% 
  st_coordinates() %>% 
  cbind(routes) %>% 
  mutate(id=paste(countrynum, statenum, Route, sep="-"),
         type="bbs")

#3. Set # of clusters---
#clusters <- 3
clusters <- 4

#B. PICK ALGORITHM####

#1. Set k for cross-fold validation----
k <- 10

#2. Pull training data using 1st bootstrap, randomize, assign fold----
dat <- track.raw %>% 
  dplyr::filter(nclust==clusters,
                boot==1,
                season=="breed") %>% 
  mutate(rand = rnorm(n())) %>% 
  arrange(rand) %>% 
  mutate(fold = ceiling(row_number()/(n()/10)))

#3. Set up k-fold----
set.seed(1234)
out.k <- list()
for(i in 1:k){
  
  dat.train <- dat %>% 
    dplyr::filter(fold != i)
  
  dat.test <- dat %>% 
    dplyr::filter(fold==i)
  
  #4. KNN----
  knn.k <- dat.test %>% 
    cbind(data.frame(knn1 = knn(train=dat.train[,c("X", "Y")], test=dat.test[,c("X", "Y")], cl=dat.train$kdecluster, k=1, prob=TRUE))) %>% 
    cbind(data.frame(knn2 = knn(train=dat.train[,c("X", "Y")], test=dat.test[,c("X", "Y")], cl=dat.train$kdecluster, k=2, prob=TRUE))) %>% 
    cbind(data.frame(knn3 = knn(train=dat.train[,c("X", "Y")], test=dat.test[,c("X", "Y")], cl=dat.train$kdecluster, k=3, prob=TRUE))) %>% 
    cbind(data.frame(knn4 = knn(train=dat.train[,c("X", "Y")], test=dat.test[,c("X", "Y")], cl=dat.train$kdecluster, k=4, prob=TRUE))) %>% 
    cbind(data.frame(knn5 = knn(train=dat.train[,c("X", "Y")], test=dat.test[,c("X", "Y")], cl=dat.train$kdecluster, k=5, prob=TRUE))) %>% 
    cbind(data.frame(knn6 = knn(train=dat.train[,c("X", "Y")], test=dat.test[,c("X", "Y")], cl=dat.train$kdecluster, k=6, prob=TRUE))) %>% 
    cbind(data.frame(knn7 = knn(train=dat.train[,c("X", "Y")], test=dat.test[,c("X", "Y")], cl=dat.train$kdecluster, k=7, prob=TRUE))) %>% 
    cbind(data.frame(knn8 = knn(train=dat.train[,c("X", "Y")], test=dat.test[,c("X", "Y")], cl=dat.train$kdecluster, k=8, prob=TRUE))) %>% 
    cbind(data.frame(knn9 = knn(train=dat.train[,c("X", "Y")], test=dat.test[,c("X", "Y")], cl=dat.train$kdecluster, k=9, prob=TRUE))) %>% 
    cbind(data.frame(knn10 = knn(train=dat.train[,c("X", "Y")], test=dat.test[,c("X", "Y")], cl=dat.train$kdecluster, k=10, prob=TRUE))) %>% 
    cbind(data.frame(knn20 = knn(train=dat.train[,c("X", "Y")], test=dat.test[,c("X", "Y")], cl=dat.train$kdecluster, k=20, prob=TRUE)))
  
  #5. Rpart classification tree----
  set.seed(i)
  rpart.k <- rpart::rpart(as.factor(kdecluster) ~ X + Y, data=dat.train, method="class") %>% 
    predict(dat.test, type="class") %>% 
    data.frame() %>% 
    rename(rpart = '.') %>% 
    cbind(dat.test)
  
  #6. Party conditional inference tree----
  set.seed(i)
  party.k <- ctree(as.factor(kdecluster) ~ X + Y, data=dat.train) %>% 
    predict(dat.test) %>% 
    data.frame() %>% 
    rename(party = '.') %>% 
    cbind(dat.test)
  
  #7. Naive bayes----
  set.seed(i)
  nb.k <- NaiveBayes(as.factor(kdecluster) ~ X + Y, data=dat.train) %>% 
    predict(dat.test) %>% 
    data.frame() %>% 
    dplyr::select(class) %>% 
    rename(nb = class) %>% 
    cbind(dat.test)
  
  #8. Linear discriminant analysis----
  set.seed(i)
  lda.k <- train(as.factor(kdecluster) ~ X + Y, data=dat.train, method="lda") %>% 
    predict(dat.test) %>% 
    data.frame() %>% 
    rename(lda = '.') %>% 
    cbind(dat.test)
  
  #9. Random forest----
  set.seed(i)
  rf.k <- train(as.factor(kdecluster) ~ X + Y, data=dat.train, method="rf") %>% 
    predict(dat.test) %>% 
    data.frame() %>% 
    rename(rf = '.') %>% 
    cbind(dat.test)
  
  #10. Support vector machine----
  set.seed(i)
  svm.k <- train(as.factor(kdecluster) ~ X + Y, data=dat.train, method="svmLinear") %>% 
    predict(dat.test) %>% 
    data.frame() %>% 
    rename(svm = '.') %>% 
    cbind(dat.test)
  
  #11. Put together----
  set.seed(i)
  out.k[[i]] <- full_join(knn.k, rpart.k) %>% 
    full_join(party.k) %>% 
    full_join(nb.k) %>% 
    full_join(lda.k) %>% 
#    full_join(other.k) %>% 
    full_join(rf.k) %>% 
    full_join(svm.k)
  
}

#12. Evaluate match----
out <- rbindlist(out.k) %>% 
  dplyr::select(-nclust, -boot, -season, -X, -Y, -rand) %>% 
  pivot_longer(cols=knn1:svm, names_to="method", values_to="pred") %>% 
  mutate(match = ifelse(kdecluster==pred, 1, 0))

#13. Summarize----
out %>% 
  group_by(method) %>% 
  summarize(correct=sum(match),
            total = n(),
            percent = correct/total) %>% 
  arrange(-percent)
#svm best

#C. PREDICT TO BBS DATA####

#1. Set up bootstrap loop----
boot <- max(track.raw$boot)

set.seed(1234)
knn.out <- list()
svm.out <- list()
for(i in 1:boot){
  
  #2. Wrangle tracking data----
  track.i <- track.raw %>% 
    dplyr::filter(nclust==clusters,
                  boot==i,
                  season=="breed")
  
  #3. Put together with BBS data----
  bbs.i <- bbs.sf %>% 
    dplyr::select(id, X, Y, type)
  
  #4. KNN----
  knn.i <- knn(train=track.i[,c("X", "Y")], test=bbs.i[,c("X", "Y")], cl=track.i$kdecluster, k=1, prob=TRUE)
  
  knn.out[[i]] <- data.frame(knncluster=knn.i,
                             knnprob=attr(knn.i, "prob"),
                             boot=i) %>% 
    cbind(bbs.i)
  
  #5. Naive bayes----
  set.seed(i)
  svm.out[[i]] <- train(as.factor(kdecluster) ~ X + Y, data=track.i, method="svmLinear") %>% 
    predict(bbs.i) %>% 
    data.frame() %>% 
    rename(svmcluster = '.') %>% 
    cbind(bbs.i)
  
  print(paste0("Finished bootstrap ", i, " of ", boot))
  
}

#9. Collapse results----
knn.all <- rbindlist(knn.out)
svm.all <- rbindlist(svm.out)

#10. Look at variation in cluster membership----
knn.sum <- knn.all %>% 
  group_by(id, knncluster) %>% 
  summarize(n=n()) %>% 
  ungroup()
summary(knn.sum$n)

svm.sum <- svm.all %>% 
  group_by(id, svmcluster) %>% 
  summarize(n=n()) %>% 
  ungroup()
summary(svm.sum$n)

knn.99 <- knn.sum %>% 
  dplyr::filter(n < 100)
table(knn.99$id)
#Only 50 points that have 1 instance of cluster variation

svm.99 <- svm.sum %>% 
  dplyr::filter(n < 100)
table(svm.99$id)

#11. Pick mean dominant cluster ID for each route----
knn.final <- knn.all %>% 
  group_by(id, X, Y, type) %>% 
  summarize(knncluster = round(mean(as.numeric(knncluster)))) %>% 
  ungroup()

svm.final <- svm.all %>% 
  group_by(id, X, Y, type) %>% 
  summarize(svmcluster = round(mean(as.numeric(svmcluster)))) %>% 
  ungroup()

#14. Visualize----
track.final <- track.raw %>% 
  dplyr::filter(season=="breed")

ggplot(knn.all) +
  geom_point(aes(x=X, y=Y, colour=knncluster)) +
  geom_point(data=track.final, aes(x=X, y=Y, colour=factor(kdecluster)), pch=21, fill="white", size=3)

ggplot(svm.all) +
  geom_point(aes(x=X, y=Y, colour=factor(svmcluster))) +
  geom_point(data=track.final, aes(x=X, y=Y, colour=factor(kdecluster)), pch=21, fill="white", size=3) +
  facet_wrap(~kdecluster)

#15. Calculate MCP area for BBSbayes-----
sp <- SpatialPointsDataFrame(coords=cbind(knn.all$X, knn.all$Y), 
                               data=data.frame(ID=knn.all$knncluster),
                               proj4string = CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"))

mp <- mcp(sp[,1], percent=100)

knn.area <- data.frame(mp) %>% 
  rename(knncluster=id)

#15. Save out----
write.csv(knn.all, "Data/LBCUBBSClustersAll.csv", row.names = FALSE)
write.csv(knn.area, "Data/area_weight.csv", row.names = FALSE)
write.csv(knn.final, "Data/LBCUBBSClusters.csv", row.names = FALSE)
