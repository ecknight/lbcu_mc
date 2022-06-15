library(tidyverse)
library(lme4)
library(MuMIn)
library(MASS)

#1. Read in data-----
dat <- read.csv("Data/LBCU_environvars_pt.csv") %>% 
  mutate(region = case_when(kdecluster==1 ~ "east", 
                            kdecluster==2 ~ "west"),
         conv = ifelse((convcrop==1 | convbuilt==1), 1, 0),
         convcrop = as.factor(convcrop),
         convbuilt = as.factor(convbuilt),
         conv = as.factor(conv))

#2. Visualize----
ggplot(dat) +
  geom_boxplot(aes(x=season, y=log(built), fill=region))

ggplot(dat) +
  geom_boxplot(aes(x=season, y=crops, fill=region))

ggplot(dat) +
  geom_boxplot(aes(x=season, y=grass, fill=region))

ggplot(dat) +
  geom_boxplot(aes(x=season, y=log(flooded_vegetation), fill=region))

ggplot(dat) +
  geom_boxplot(aes(x=season, y=change_norm, fill=region))

ggplot(dat) +
  geom_boxplot(aes(x=season, y=occurrence, fill=region))

ggplot(dat) +
  geom_boxplot(aes(x=season, y=recurrence, fill=region))

ggplot(dat) +
  geom_boxplot(aes(x=season, y=seasonality, fill=region))

ggplot(dat) +
  geom_boxplot(aes(x=season, y=raven, fill=region))

ggplot(dat) +
  geom_boxplot(aes(x=season, y=as.numeric(convcrop), fill=region))

ggplot(dat) +
  geom_boxplot(aes(x=season, y=as.numeric(convbuilt), fill=region))

ggplot(dat) +
  geom_boxplot(aes(x=season, y=as.numeric(conv), fill=region))

#3. Regression----
m1 <- lmer(crops ~ region*season + (1|id), data=dat, na.action="na.fail")
dredge(m1)

m2 <- lmer(grass ~ region*season + (1|id), data=dat, na.action="na.fail")
dredge(m2)

m3 <- lmer(built ~ region*season + (1|id), data=dat, na.action="na.fail")
dredge(m3)

m4 <- lmer(flooded_vegetation ~ region*season + (1|id), data=dat, na.action="na.fail")
dredge(m4)

m5 <- lmer(change_norm ~ region*season + (1|id), data=dat, na.action="na.fail")
dredge(m5)

m6 <- lmer(occurrence ~ region*season + (1|id), data=dat, na.action="na.fail")
dredge(m6)

m7 <- lmer(recurrence ~ region*season + (1|id), data=dat, na.action="na.fail")
dredge(m7)

m8 <- lmer(seasonality ~ region*season + (1|id), data=dat, na.action="na.fail")
dredge(m8)

m9 <- lmer(raven ~ region*season + (1|id), data=dat, na.action="na.fail")
dredge(m9)

m10 <- glmer(convcrop ~ region*season + (1|id), data=dat, na.action="na.fail", family="binomial")
dredge(m10)

m11 <- glmer(convbuilt ~ region*season + (1|id), data=dat, na.action="na.fail", family="binomial")
dredge(m11)

m12 <- glmer(conv ~ region*season + (1|id), data=dat, na.action="na.fail", family="binomial")
dredge(m12)

#4. LDA----
dat.lda <- dat %>% 
  mutate(group = paste0(region, "-", season)) %>% 
  dplyr::select(group, crops, grass, built, flooded_vegetation, change_norm, occurrence, recurrence, seasonality, raven, convcrop, convbuilt, conv)
l1 <- lda(group~., dat.lda)
