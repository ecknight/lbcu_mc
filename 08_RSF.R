library(tidyverse)
library(usdm)
library(lme4)
library(MuMIn)
library(ggridges)
library(mgcv)

options(scipen=99999)

setwd("G:/My Drive/SMBC")

#1. Load data----
raw <- read.csv("Data/LBCU_environvars_RSF.csv") 

dat <- raw |> 
  mutate(response = ifelse(type=="used", 1, 0),
         change.raw = change,
         drought.raw = drought,
         seasonality.raw = seasonality,
         group = as.factor(group)) |> 
  mutate(change = (change.raw - min(change.raw))/(max(change.raw) - min(change.raw)),
         drought = (drought.raw - min(drought.raw))/(max(drought.raw) - min(drought.raw)),
         seasonality = (seasonality.raw - min(seasonality.raw))/(max(seasonality.raw) - min(seasonality.raw)))

#2. Test VIF----
loop <- dat |> 
  dplyr::select(nclust, season) |> 
  unique()

v.list <- list()
for(i in 1:nrow(loop)){
  
  cov.i <- dat |>
    dplyr::filter(season == loop$season[i],
                  nclust == loop$nclust[i]) |>
    dplyr::select(crop, grass, wetland)
  
  v.list[[i]] <- vif(cov.i) |>
    mutate(season = loop$season[i],
           nclust = loop$nclust[i])
  
}

v <- do.call(rbind, v.list) |> 
  arrange(-VIF)
#all good

#3. Visualize----
#Check for sufficient distribution
#Check for polynomials

#wetland change
ggplot(dat) +
  geom_histogram(aes(x=change, colour=group)) +
  facet_wrap(nclust~season, scales="free")

ggplot(dat) +
#  geom_jitter(aes(x=change, y=response, colour=group)) +
  geom_smooth(aes(x=change, y=response, colour=group)) +
  facet_wrap(nclust~season)
#maybe for winter

#crop
ggplot(dat) +
  geom_histogram(aes(x=crop, colour=group)) +
  facet_wrap(nclust~season, scales="free")

ggplot(dat) +
#  geom_jitter(aes(x=crop, y=response, colour=group)) +
  geom_smooth(aes(x=crop, y=response, colour=group)) +
  facet_wrap(nclust~season)

#grass
ggplot(dat) +
  geom_histogram(aes(x=grass, colour=group)) +
  facet_wrap(nclust~season, scales="free")

ggplot(dat) +
#  geom_jitter(aes(x=grass, y=response, colour=group)) +
  geom_smooth(aes(x=grass, y=response, colour=group)) +
  facet_wrap(nclust~season)

#seasonality
ggplot(dat) +
  geom_histogram(aes(x=seasonality, colour=group)) +
  facet_wrap(nclust~season, scales="free")

ggplot(dat) +
#  geom_jitter(aes(x=seasonality, y=response, colour=group)) +
  geom_smooth(aes(x=seasonality, y=response, colour=group)) +
  facet_grid(nclust~season)

#wetland
ggplot(dat) +
  geom_histogram(aes(x=wetland, colour=group)) +
  facet_wrap(nclust~season, scales="free")

ggplot(dat) +
  #  geom_jitter(aes(x=wetland, y=response, colour=group)) +
  geom_smooth(aes(x=wetland, y=response)) +
  facet_grid(nclust~season)


#4. Set up loop----
loop <- dat |> 
  dplyr::select(nclust, season) |> 
  unique()

m.list <- list()
s.list <- list()
d.list <- list()
for(i in 1:nrow(loop)){
  
  #5. Subset data----
  dat.i <- dat |> 
    dplyr::filter(season==loop$season[i],
                  nclust==loop$nclust[i])
  
  #6. Model----
  m.list[[i]] <- glm(response ~ crop + crop:group + grass + grass:group + wetland + wetland:group, family = "binomial", data=dat.i, na.action = "na.fail")
  
  #7. Save summary----
  s.list[[i]] <- summary(m.list[[i]])$coefficients |> 
    data.frame() |> 
    mutate(season = loop$season[i],
           nclust = loop$nclust[i],
           cov = row.names(summary(m.list[[i]])$coefficients)) 
  
  #8. Dredge----
  d.list[[i]] <- data.frame(dredge(m.list[[i]])) |> 
    mutate(season = loop$season[i],
           nclust = loop$nclust[i])
  
}

#9. Wrangle output----
summary <- do.call(rbind, s.list) |> 
  dplyr::filter(cov!="(Intercept)") |> 
  rename(p = 'Pr...z..', se = 'Std..Error') |> 
  mutate(sig = ifelse(p < 0.05, 1, 0))

dredged <- do.call(rbind, d.list) |> 
  mutate(delta2 = ifelse(delta < 2, 1, 0)) |> 
  group_by(nclust, season, delta2) |> 
  mutate(mindf = min(df),
         pick = ifelse(df==mindf & delta2==1, 1, 0)) |> 
  ungroup() |> 
  mutate(weight = round(weight, 2),
         logLik = round(logLik, 2),
         AICc = round(AICc, 2),
         delta = round(delta, 2))

dredged.pick <- dredged |> 
  dplyr::filter(pick==1) |> 
  arrange(season, nclust, df) |> 
  group_by(nclust, season) |> 
  dplyr::filter(row_number()==1) |> 
  ungroup()

write.csv(dredged, "Results/RSFAIC.csv", row.names = FALSE)

#10. Run final models----

m.final <- list()
vars <- list()

for(i in 1:nrow(loop)){
  
  dat.i <- dat |> 
    dplyr::filter(season==loop$season[i],
                  nclust==loop$nclust[i])
  
  vars1 <- dredged.pick |> 
    dplyr::filter(season==loop$season[i],
                  nclust==loop$nclust[i]) |> 
    dplyr::select(crop:wetland) |> 
    pivot_longer(crop:wetland, names_to="var", values_to="est") |> 
    dplyr::filter(!is.na(est))
  
  vars2 <- dredged.pick |> 
    dplyr::filter(season==loop$season[i],
                  nclust==loop$nclust[i]) |> 
    dplyr::select('crop.group':'group.wetland') |> 
    pivot_longer('crop.group':'group.wetland', names_to="var", values_to="est") |> 
    dplyr::filter(!is.na(est)) |> 
    mutate(var = gsub(pattern="[.]", replacement=":", x=var))
  
  vars[[i]] <- c(vars1$var, vars2$var)
  
  m.base <- glm(response ~ 1, family = "binomial", data = dat.i, na.action = "na.fail")
  
  m.final[[i]] <- update(m.base, formula = as.formula(paste0("~ ", paste(vars[[i]], collapse = " + "))))
  
  
}

#11. Predictions----

pred.list <- list()
for(i in 1:nrow(loop)){
  
  dat.i <- dat |> 
    dplyr::filter(season==loop$season[i],
                  nclust==loop$nclust[i])
  
  newdat.crop <- expand.grid(group=unique(dat.i$group), crop=seq(0, 1, 0.01), grass=0, seasonality=0, change=0, wetland = 0)
  
  newdat.grass <- expand.grid(group=unique(dat.i$group), crop=0, grass=seq(0, 1, 0.01), seasonality=0, change=0, wetland = 0)
  
  newdat.wetland <- expand.grid(group=unique(dat.i$group), crop=0, grass=0, seasonality=0, change=0, wetland = seq(0, 1, 0.01))
  
  newdat.seasonality <- expand.grid(group=unique(dat.i$group), crop=0, grass=0, seasonality=seq(0, 1, 0.01), change=0, wetland = 0)
  
  newdat.change = expand.grid(group=unique(dat.i$group), crop=0, grass=0, seasonality=0, change=seq(0, 1, 0.01), wetland = 0)
  
  pred.crop <- predict(m.final[[i]], newdat.crop, se=TRUE, type="response")  |> 
    data.frame() |> 
    cbind(newdat.crop |> 
            dplyr::select(group, crop)) |> 
    rename(val=crop) |> 
    mutate(up=fit+1.96*se.fit, low=fit-1.96*se.fit,
           cov="crop")
  
  pred.grass <- predict(m.final[[i]], newdat.grass, se=TRUE, type="response")  |> 
    data.frame() |> 
    cbind(newdat.grass |> 
            dplyr::select(group, grass)) |> 
    rename(val=grass) |> 
    mutate(up=fit+1.96*se.fit, low=fit-1.96*se.fit,
           cov="grass")
  
  pred.seasonality <- predict(m.final[[i]], newdat.seasonality, se=TRUE, type="response")  |> 
    data.frame() |> 
    cbind(newdat.seasonality |> 
            dplyr::select(group, seasonality)) |> 
    rename(val=seasonality) |> 
    mutate(up=fit+1.96*se.fit, low=fit-1.96*se.fit,
           cov="seasonality")
  
  pred.change <- predict(m.final[[i]], newdat.change, se=TRUE, type="response")  |> 
    data.frame() |> 
    cbind(newdat.change |> 
            dplyr::select(group, change)) |> 
    rename(val=change) |> 
    mutate(up=fit+1.96*se.fit, low=fit-1.96*se.fit,
           cov="change")
  
  pred.wetland <- predict(m.final[[i]], newdat.change, se=TRUE, type="response")  |> 
    data.frame() |> 
    cbind(newdat.wetland |> 
            dplyr::select(group, wetland)) |> 
    rename(val=wetland) |> 
    mutate(up=fit+1.96*se.fit, low=fit-1.96*se.fit,
           cov="wetland")
  
  vars.i <- vars[[i]]
  
  pred.list[[i]] <- rbind(pred.crop, pred.grass, pred.seasonality, pred.change, pred.wetland) |> 
    mutate(season=loop$season[i],
           nclust = loop$nclust[i]) |> 
    dplyr::filter(cov %in% vars.i)
    
  
}

pred <- do.call(rbind, pred.list) |> 
  # mutate(group = ifelse((season=="breed" & cov %in% c("crop", "grass")) |
  #                          (season=="fallmig" & cov=="grass") |
  #                          (season=="springmig" & cov=="seasonality"),
  #                        "No group\ndifference",
  #                        as.character(group))) |> 
  mutate(val.use = case_when(cov=="change" ~ val*(max(dat$change.raw)-min(dat$change.raw)) + min(dat$change.raw),
                         cov=="seasonality" ~ val*(max(dat$seasonality.raw)-min(dat$seasonality.raw)) + min(dat$seasonality.raw),
                         !is.na(val) ~ val))

write.csv(pred, "Results/RSFPredictions.csv", row.names = FALSE)

#12. Plot----
plot.3 <- ggplot(pred |> 
         dplyr::filter(nclust=="3")) +
  geom_ribbon(aes(x=val.use, ymin=low, ymax=up, group=group), alpha=0.3) +
  geom_line(aes(x=val.use, y=fit, colour=group)) +
  facet_grid(season ~ cov, scales="free") +
  theme(legend.position = "bottom")

plot.expert <- ggplot(pred |> 
         dplyr::filter(nclust=="expert")) +
  geom_ribbon(aes(x=val.use, ymin=low, ymax=up, group=group), alpha=0.3) +
  geom_line(aes(x=val.use, y=fit, colour=group)) +
  facet_grid(season ~ cov, scales="free") +
  theme(legend.position = "bottom")

plot.flyway <- ggplot(pred |> 
         dplyr::filter(nclust=="flyway")) +
  geom_ribbon(aes(x=val.use, ymin=low, ymax=up, group=group), alpha=0.3) +
  geom_line(aes(x=val.use, y=fit, colour=group)) +
  facet_grid(season ~ cov, scales="free") +
  theme(legend.position = "bottom")

gridExtra::grid.arrange(plot.3, plot.expert, plot.flyway, ncol = 3)
