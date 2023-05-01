#1. Extract index data - tweaked for novel strata

extract_index_data <- function (jags_mod = NULL, alt_n = "n", jags_data = NULL) 
{
  stratify_by <- jags_mod$stratify_by
  
  n <- jags_mod$sims.list[[alt_n]]

  bugs_data <- dat.j

  y_min <- bugs_data$ymin
  y_max <- bugs_data$ymax
  strat_list = unique(data.frame(strat_name = jags_mod$strat_name, 
                                 strat = bugs_data$strat, stringsAsFactors = FALSE))
  strat_list = strat_list[order(strat_list$strat), ]
  strata_used <- strat_list$strat_name
  strata_num <- strat_list$strat
  area_weights <- all_area_weights[which(all_area_weights$region %in% 
                                           strata_used), ]
  area_weights <- area_weights[order(match(area_weights$region, 
                                           strata_used)), ]
  area_weights$num <- strata_num

  original_data = get_prepared_data(jags_data = jags_data)
  
  return(list(n = n, area_weights = area_weights, y_min = y_min, 
              y_max = y_max, r_year = jags_mod$r_year, bugs_data = bugs_data, 
              original_data = original_data))
}



#2. Generate indices

generate_indices <- function (jags_mod = NULL, jags_data = NULL, quantiles = c(0.025, 
                                                           0.05, 0.25, 0.75, 0.95, 0.975), regions = c("stratum", "continental"), 
          alternate_n = "n", startyear = NULL, drop_exclude = FALSE, 
          max_backcast = NULL, alt_region_names = NULL) 
{
  data_list <- extract_index_data(jags_mod = jags_mod, 
                                  alt_n = alternate_n,
                                  jags_data = jags_data)

  n <- data_list$n
  stratify_by <- jags_mod$stratify_by
  original_data = data_list$original_data

  area_weights <- data_list$area_weights
  area_weights$region = as.character(area_weights$region)
  
  y_min <- data_list$y_min
  y_max <- data_list$y_max
  mr_year <- min(data_list$r_year)
  max_backcast <- length(y_min:y_max)

  bugs_data = data_list$bugs_data
  rawall = data.frame(year = bugs_data$year, count = bugs_data$count, 
                      strat = bugs_data$strat)
  rawnonz = rawall[which(rawall$count > 0), ]
  fyearbystrat = tapply(rawnonz$year, rawnonz[, c("strat")], 
                        min, na.rm = TRUE)
  raw = rawall[which(rawall$year >= y_min), ]
  original_data <- original_data[which(original_data$Year_Factored >= 
                                         y_min), ]
  nrts_total_by_strat <- tapply(original_data$Route, original_data$Stratum_Factored, 
                                FUN = function(x) {
                                  length(unique(x))
                                })
  non_zero_weight = bugs_data$nonzeroweight
  n_samples <- dim(n)[1]
  
  region_names <- all_area_weights
  region_names$stratum = region_names$region
  region_names$continental = "Continental"
  data_summary <- data.frame(Year = integer(), Region = character(), 
                             Region_alt = character(), Region_type = character(), 
                             Strata_included = character(), Strata_excluded = character(), 
                             Index = double(), stringsAsFactors = FALSE)
  for (qq in quantiles) {
    data_summary[, paste0("Index_q_", qq)] <- double()
  }
  data_summary$obs_mean <- double()
  data_summary$nrts <- integer()
  data_summary$nnzero <- integer()
  data_summary$nrts_total <- integer()
  N_all <- list()
  n_index <- 0
  for (rr in regions) {
    rrall = unique(region_names[, rr])
    col_region_name <- rr
    for (rrs in rrall) {
      region_alt_name <- as.character(unique(region_names[which(region_names[, 
                                                                             rr] == rrs), col_region_name]))
      region_alt_name = paste("region", region_alt_name, 
                                sep = "_")
      st_sela <- as.character(region_names[which(region_names[, 
                                                              rr] == rrs), "region"])
      st_rem <- NULL
      strata_sel <- area_weights[which(area_weights$region %in% 
                                         st_sela), "num"]
      st_sel <- area_weights[which(area_weights$region %in% 
                                     st_sela), "region"]
      pz_area <- area_weights[, "area_sq_km"] * non_zero_weight
      if (length(strata_sel) < 1) {
        next
      }
      obs_df = data.frame(year = integer(), strat = integer(), 
                          obs_mean = double(), nrts = integer(), nnzero = integer(), 
                          nrts_total = integer(), strata_rem_flag = double())
      for (j in strata_sel) {
        rawst <- raw[which(raw$strat == j), c("year", 
                                              "count")]
        yrs <- data.frame(year = c(y_min:y_max))
        rawst <- merge(rawst, yrs, by = "year", all = TRUE)
        rawst <- rawst[order(rawst$year), ]
        o_mns <- as.numeric(by(rawst[, 2], INDICES = rawst[, 
                                                           1], FUN = mean, na.rm = TRUE))
        nrts <- as.numeric(by(rawst[, 2], INDICES = rawst[, 
                                                          1], FUN = function(x) {
                                                            length(which(!is.na(x)))
                                                          }))
        nnzero <- as.numeric(by(rawst[, 2], INDICES = rawst[, 
                                                            1], FUN = function(x) {
                                                              length(which(x > 0))
                                                            }))
        strata_p <- pz_area[j]/sum(pz_area[strata_sel])
        if (sum(nnzero[1:max_backcast]) < 1 & as.integer(fyearbystrat[j]) > 
            y_min) {
          st_rem <- c(st_rem, as.character(area_weights[which(area_weights$num == 
                                                                j), "region"]))
          strem_flag <- c(rep(strata_p, fyearbystrat[j] - 
                                (y_min - 1)), rep(0, y_max - fyearbystrat[j]))
        }
        else {
          strem_flag <- rep(0, length(y_max:y_min))
        }
        obs_df_t <- data.frame(year = c(y_min:y_max), 
                               strat = j, obs_mean = o_mns * ((area_weights$area_sq_km[which(area_weights$num == 
                                                                                               j)])/sum(area_weights[which(area_weights$num %in% 
                                                                                                                             strata_sel), "area_sq_km"])) * (non_zero_weight[j]), 
                               nrts = nrts, nnzero = nnzero, nrts_total = as.integer(nrts_total_by_strat[j]), 
                               strata_rem_flag = strem_flag)
        obs_df <- rbind(obs_df, obs_df_t)
      }
      if (!is.null(st_rem)) {
        if (drop_exclude) {
          strata_sel <- strata_sel[-which(strata_sel %in% 
                                            area_weights[which(area_weights$region %in% 
                                                                 st_rem), "num"])]
          st_sel <- st_sel[-which(st_sel %in% st_rem)]
        }
      }
      if (length(strata_sel) < 1) {
        next
      }
      n_weight <- array(NA, dim = c(dim(n)[c(1, 2)], length(y_min:y_max)))
      n_weight[, , 1:length(y_min:y_max)] <- n[, , y_min:y_max]
      for (i in 1:n_samples) {
        for (j in strata_sel) {
          n_weight[i, j, ] <- (n_weight[i, j, ] * area_weights$area_sq_km[which(area_weights$num == 
                                                                                  j)])/sum(area_weights[which(area_weights$num %in% 
                                                                                                                strata_sel), "area_sq_km"])
        }
      }
      n_weight <- n_weight[, strata_sel, ]
      if (length(strata_sel) > 1) {
        N <- apply(n_weight, c(1, 3), sum)
      }
      else {
        N <- n_weight
      }
      n_index <- n_index + 1
      N_all[[n_index]] <- N
      names(N_all)[n_index] <- paste(rr, rrs, sep = "_")
      n_median <- apply(N, 2, median)
      data_summaryr <- data.frame(Year = seq(y_min:y_max), 
                                  Region = rrs, Region_alt = region_alt_name, Region_type = rr, 
                                  Strata_included = paste(st_sel, collapse = " ; "), 
                                  Strata_excluded = paste(st_rem, collapse = " ; "), 
                                  Index = n_median, stringsAsFactors = FALSE)
      for (qq in quantiles) {
        data_summaryr[, paste0("Index_q_", qq)] <- apply(N, 
                                                         2, stats::quantile, probs = qq)
      }
      data_summaryr$Year <- as.integer((data_summaryr$Year - 
                                          1) + mr_year)
      data_summaryr$obs_mean <- as.numeric(by(obs_df[, 
                                                     3], INDICES = obs_df[, 1], FUN = sum, na.rm = TRUE))
      data_summaryr$nrts <- as.numeric(by(obs_df[, 4], 
                                          INDICES = obs_df[, 1], FUN = sum, na.rm = TRUE))
      data_summaryr$nnzero <- as.numeric(by(obs_df[, 5], 
                                            INDICES = obs_df[, 1], FUN = sum, na.rm = TRUE))
      data_summaryr$nrts_total <- as.numeric(by(obs_df[, 
                                                       6], INDICES = obs_df[, 1], FUN = sum, na.rm = TRUE))
      data_summaryr$backcast_flag <- 1 - as.numeric(by(obs_df[, 
                                                              7], INDICES = obs_df[, 1], FUN = sum, na.rm = TRUE))
      data_summary = rbind(data_summary, data_summaryr)
    }
  }
  return(list(data_summary = data_summary, samples = N_all, 
              area_weights = area_weights, y_min = y_min, y_max = y_max, 
              startyear = mr_year, regions = regions, raw_data = raw))
}
