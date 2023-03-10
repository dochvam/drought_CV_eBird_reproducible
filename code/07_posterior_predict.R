###############################################################################
# 07_posterior_predict.R
# Author: Ben Goldstein
# Description: This script conducts posterior prediction methodology. Using
#      posterior samples generated during N-mixture model and LMM estimation,
#      an abundance is sampled for each species under drought and non-drought
#      conditions. Interaction terms and diagnostics are also monitored.
###############################################################################


library(nimble)
library(Matrix)
library(tidyverse)
library(abind)

source("code/make_data.R")

use_qnbinom <- FALSE

# Make a data frame of the combos of variables representing diff. counterfactual
# scenarios
conditions <- data.frame(
  name = c("All", "EVI only", "NDWI only", "ppt only", "tmax only", "ppt and tmax"),
  evi  = c(1, 1, 0, 0, 0, 0),
  ndwi = c(1, 0, 1, 0, 0, 0),
  ppt  = c(1, 0, 0, 1, 0, 1),
  tmax = c(1, 0, 0, 0, 1, 1)
)


site_hab_props <- all_model_data %>% 
  select(cellnum, all_of(land_covars)) %>% 
  distinct()


# Scale variables to match the scale used in model fitting
site_hab_props$row_field_ag <-
  site_hab_props$row_field_ag * attr(site_hab_props$row_field_ag, 'scaled:scale') +
  attr(site_hab_props$row_field_ag, 'scaled:center')
site_hab_props$perennial_ag <-
  site_hab_props$perennial_ag * attr(site_hab_props$perennial_ag, 'scaled:scale') +
  attr(site_hab_props$perennial_ag, 'scaled:center')
site_hab_props$riparian <-
  site_hab_props$riparian * attr(site_hab_props$riparian, 'scaled:scale') +
  attr(site_hab_props$riparian, 'scaled:center')
site_hab_props$grass_pasture <-
  site_hab_props$grass_pasture * attr(site_hab_props$grass_pasture, 'scaled:scale') +
  attr(site_hab_props$grass_pasture, 'scaled:center')
site_hab_props$other_natural <-
  site_hab_props$other_natural * attr(site_hab_props$other_natural, 'scaled:scale') +
  attr(site_hab_props$other_natural, 'scaled:center')
site_hab_props$other <- 1 - rowSums(site_hab_props[, land_covars])


# For each result file....
all_result_files <-
  list.files("intermediate/SSAMs", pattern = "MCMC_SSAM", full.names = TRUE)
target_files <- c()

# Loop over species
for (i in 1:length(species_vector)) {
  thisfiles <- all_result_files[grepl(species_vector[i], all_result_files)]
  if (length(thisfiles) > 1) {
    target_files <- c(target_files, thisfiles[
      grepl("rerun", thisfiles)
    ])
  } else {
    target_files <- c(target_files, thisfiles)
  }
}

specvec <- substr(target_files, 30, 33)

# Prepare (non-sampled) data
pp_data <- all_model_data %>% 
  mutate(cellID = as.numeric(as.factor(cellnum))) %>% 
  select(cellnum, cellID, intercept_abd, 
         row_field_ag, perennial_ag, riparian, grass_pasture, other_natural,
         all_of(paste0(land_covars, "_unscaled"))) %>% 
  mutate(other_unscaled = 1 - (row_field_ag_unscaled + 
           perennial_ag_unscaled + riparian_unscaled + 
           grass_pasture_unscaled + other_natural_unscaled)) %>% 
  distinct() %>% 
  mutate(
    year1 = 0, year2 = 0, year3 = 0, year4 = 0, year5 = 0,
    year6 = 0, year7 = 0, year8 = 0, year9 = 0,
    intercept_abd = 1
  )
site_hab_props <- site_hab_props[match(pp_data$cellnum, site_hab_props$cellnum),]


#### Read in posterior samples of four drought-varying covars ####
read_ppred <- function(x, scale_col, cellnum_vec, thin = 10) {
  temp <- readRDS(x)
  stopifnot(all(temp$data_drought$cellnum == temp$data_nodrought$cellnum))
  
  # Get the sites we want based on cellnum_vec
  cols <- match(cellnum_vec, temp$data_drought$cellnum)
  
  # Apply thinning
  temp$ppred_drought <- 
    temp$ppred_drought[(1:(nrow(temp$ppred_drought)/thin))*thin, cols]
  temp$ppred_nodrought <- 
    temp$ppred_nodrought[(1:(nrow(temp$ppred_nodrought)/thin))*thin, cols]
  
  # Transform to modeled scale
  temp$ppred_drought[]   <- (temp$ppred_drought[] - attr(scale_col, "scaled:center")) / 
                              attr(scale_col, "scaled:scale")
    
  temp$ppred_nodrought[] <- (temp$ppred_nodrought[] - attr(scale_col, "scaled:center")) / 
                               attr(scale_col, "scaled:scale")
  
  
  temp
}

ndwi_post <- read_ppred("intermediate/NDWI_LMM_posterior.RDS", 
                        scale_col = all_model_data$NDWI,
                        cellnum_vec = pp_data$cellnum)
evi_post  <- read_ppred("intermediate/EVI_LMM_posterior.RDS", 
                        scale_col = all_model_data$EVI,
                        cellnum_vec = pp_data$cellnum)
tmax_post <- read_ppred("intermediate/tmax_ssn_LMM_posterior.RDS", 
                        scale_col = all_model_data$tmax_ssn,
                        cellnum_vec = pp_data$cellnum)
ppt_post  <- read_ppred("intermediate/ppt_LMM_posterior.RDS", 
                        scale_col = all_model_data$ppt,
                        cellnum_vec = pp_data$cellnum)

#### LOOP OVER COMPARISON CONDITIONS ####
for (cond_iter in 1:nrow(conditions)) {
  cat("Running posterior prediction for condition:", conditions$name[cond_iter], "\n")
  pb <- progress::progress_bar$new(total = length(specvec))
  
  diagnostics_df_list <- list()
  spec_est_list <- list()
  spec_habcount_df_list <- list()
  spec_habint_df_list <- list()
  spec_habint_summary_list <- list()
  
#### LOOP OVER SPECIES ####
for (k in 1:length(target_files)) {
  
  thisspec <- specvec[k]
  
  # Read in this species' results
  res <- readRDS(target_files[k])
  samples_mtx <- do.call(rbind, lapply(res$samples, as.matrix))
  
  params <- colnames(res$samples[[1]])
  colnames(samples_mtx) <- params
  
  nsamp <- nrow(samples_mtx)
  nsite <- nrow(pp_data)
  
  site_ranef_inds <- which(grepl("site_ranef\\[", params))
  site_ranef_cellIDs <- parse_number(params[site_ranef_inds])
  stopifnot(all(site_ranef_cellIDs == 1:nsite))
  site_ranef_mtx <- samples_mtx[, site_ranef_inds]
  
  site_ranef_mtx <- site_ranef_mtx[, pp_data$cellID]
  
  # theta
  theta_index <- which(params == "theta")
  
  # Calculate lambda at each site, for a given counterfactual dataset
  #    Outcome is a nsite x nsample data frame
  
  post_samps <- rep(1:nrow(evi_post$ppred_nodrought), 
                    ceiling(nrow(samples_mtx) / nrow(evi_post$ppred_nodrought)))[
    1:nrow(samples_mtx)
  ]
  
  ### Set up data matrices
  # Pull out individual cover covariates
  land_mtx2 <- matrix(unlist(pp_data[, lambda_covars[2]]), nrow = nsamp, ncol = nsite, byrow = TRUE)
  land_mtx3 <- matrix(unlist(pp_data[, lambda_covars[3]]), nrow = nsamp, ncol = nsite, byrow = TRUE)
  land_mtx4 <- matrix(unlist(pp_data[, lambda_covars[4]]), nrow = nsamp, ncol = nsite, byrow = TRUE)
  land_mtx5 <- matrix(unlist(pp_data[, lambda_covars[5]]), nrow = nsamp, ncol = nsite, byrow = TRUE)
  land_mtx6 <- matrix(unlist(pp_data[, lambda_covars[6]]), nrow = nsamp, ncol = nsite, byrow = TRUE)
  # Baseline climate matrices:
  evi_baseline_mtx  <- evi_post$ppred_nodrought[ post_samps, ]
  ndwi_baseline_mtx <- ndwi_post$ppred_nodrought[post_samps, ]
  tmax_baseline_mtx <- tmax_post$ppred_nodrought[post_samps, ]
  ppt_baseline_mtx  <- ppt_post$ppred_nodrought[ post_samps, ]
  # Choose alt climate values based on "condition"
  if (conditions$evi[cond_iter] == 1) {
    evi_alt_mtx  <- evi_post$ppred_drought[post_samps, ]
  } else {
    evi_alt_mtx <- evi_baseline_mtx
  }
  
  if (conditions$ndwi[cond_iter] == 1) {
    ndwi_alt_mtx  <- ndwi_post$ppred_drought[post_samps, ]
  } else {
    ndwi_alt_mtx <- ndwi_baseline_mtx
  }
  
  if (conditions$tmax[cond_iter] == 1) {
    tmax_alt_mtx  <- tmax_post$ppred_drought[post_samps, ]
  } else {
    tmax_alt_mtx <- tmax_baseline_mtx
  }
  
  if (conditions$ppt[cond_iter] == 1) {
    ppt_alt_mtx  <- ppt_post$ppred_drought[post_samps, ]
  } else {
    ppt_alt_mtx <- ppt_baseline_mtx
  }
  
  
  
  # data array: nsamp x nsite x ncovar
  data_arr_nodrought <- abind(
    matrix(1, nrow = nsamp, ncol = nsite), # Intercept
    
    land_mtx2, # Land covars
    land_mtx3, # Land covars
    land_mtx4, # Land covars
    land_mtx5, # Land covars
    land_mtx6, # Land covars
    evi_baseline_mtx, # EVI
    ndwi_baseline_mtx, # NDWI
    tmax_baseline_mtx, # tmax_ssn
    ppt_baseline_mtx, # ppt
    
    # NDWI * tmax interaction
    ndwi_baseline_mtx * tmax_baseline_mtx,
    
    # NDWI interactions:
    ndwi_baseline_mtx * land_mtx2,
    ndwi_baseline_mtx * land_mtx3,
    ndwi_baseline_mtx * land_mtx4,
    ndwi_baseline_mtx * land_mtx5,
    ndwi_baseline_mtx * land_mtx6,
    # EVI interactions:
    evi_baseline_mtx * land_mtx2,
    evi_baseline_mtx * land_mtx3,
    evi_baseline_mtx * land_mtx4,
    evi_baseline_mtx * land_mtx5,
    evi_baseline_mtx * land_mtx6,
    # tmax interactions:
    tmax_baseline_mtx * land_mtx2,
    tmax_baseline_mtx * land_mtx3,
    tmax_baseline_mtx * land_mtx4,
    tmax_baseline_mtx * land_mtx5,
    tmax_baseline_mtx * land_mtx6,    
    # ppt interactions:
    ppt_baseline_mtx * land_mtx2,
    ppt_baseline_mtx * land_mtx3,
    ppt_baseline_mtx * land_mtx4,
    ppt_baseline_mtx * land_mtx5,
    ppt_baseline_mtx * land_mtx6,
    
    along = 3
  )
  
  data_arr_drought <- abind(
    matrix(1, nrow = nsamp, ncol = nsite), # Intercept
    
    land_mtx2, # Land covars
    land_mtx3, # Land covars
    land_mtx4, # Land covars
    land_mtx5, # Land covars
    land_mtx6, # Land covars
    evi_alt_mtx, # EVI
    ndwi_alt_mtx, # NDWI
    tmax_alt_mtx, # tmax_ssn
    ppt_alt_mtx, # ppt
    
    # NDWI * tmax interaction
    ndwi_alt_mtx * tmax_alt_mtx,
    
    # NDWI interactions:
    ndwi_alt_mtx * land_mtx2,
    ndwi_alt_mtx * land_mtx3,
    ndwi_alt_mtx * land_mtx4,
    ndwi_alt_mtx * land_mtx5,
    ndwi_alt_mtx * land_mtx6,
    # EVI interactions:
    evi_alt_mtx * land_mtx2,
    evi_alt_mtx * land_mtx3,
    evi_alt_mtx * land_mtx4,
    evi_alt_mtx * land_mtx5,
    evi_alt_mtx * land_mtx6,
    # tmax interactions:
    tmax_alt_mtx * land_mtx2,
    tmax_alt_mtx * land_mtx3,
    tmax_alt_mtx * land_mtx4,
    tmax_alt_mtx * land_mtx5,
    tmax_alt_mtx * land_mtx6,
    # ppt interactions:
    ppt_alt_mtx * land_mtx2,
    ppt_alt_mtx * land_mtx3,
    ppt_alt_mtx * land_mtx4,
    ppt_alt_mtx * land_mtx5,
    ppt_alt_mtx * land_mtx6,
    
    
    along = 3
  )
  
  coeff_mtx <- samples_mtx[, paste0("lambda_beta[", 1:31, "]")]  
  
  log_lambda_samples_drought <- log_lambda_samples_nodrought <-
    matrix(NA, nrow = nsite, ncol = nsamp)
  
  #### Calculate lambda and get draws of N ####
  for (i in 1:nsite) {
    log_lambda_samples_drought[i, ]   <- 
      rowSums(data_arr_drought[, i, ] * coeff_mtx) + site_ranef_mtx[, i]
    log_lambda_samples_nodrought[i, ] <- 
      rowSums(data_arr_nodrought[, i, ] * coeff_mtx) + site_ranef_mtx[, i]
  }
  
  
  N_drought <- N_nodrought <- log_lambda_samples_drought
  llam_bounds <- quantile(log_lambda_samples_drought, probs = c(0.0025, 0.9975))

  for (i in 1:ncol(N_drought)) {
    if (use_qnbinom) {
      quants <- runif(n = nrow(log_lambda_samples_drought), 0, 1)
      N_drought[, i] <- qnbinom(
        p = quants,
        size = 1 / samples_mtx[i, theta_index],
        mu = exp(log_lambda_samples_drought[, i])
      )
      N_nodrought[, i] <- qnbinom(
        p = quants,
        size = 1 / samples_mtx[i, theta_index],
        mu = exp(log_lambda_samples_nodrought[, i])
      )
    } else {
      N_drought[, i] <- rnbinom(
        n = nrow(log_lambda_samples_drought),
        size = 1 / samples_mtx[i, theta_index],
        mu = exp(log_lambda_samples_drought[, i])
      )
      N_nodrought[, i] <- rnbinom(
        n = nrow(log_lambda_samples_drought),
        size = 1 / samples_mtx[i, theta_index],
        mu = exp(log_lambda_samples_nodrought[, i])
      )
    }
  }
  
  # Overall count in drought vs. non-drought
  N_drought_total <- colSums(N_drought)
  N_nodrought_total <- colSums(N_nodrought)
  
  if (cond_iter == 1) {
    Ndrought_byhab <- matrix(nrow = ncol(N_drought), ncol = 6)
    Nnodrought_byhab <- matrix(nrow = ncol(N_drought), ncol = 6)
    colnames(Ndrought_byhab) <- c(land_covars, "other")  
    colnames(Nnodrought_byhab) <- c(land_covars, "other")  
    
    for (i in 1:ncol(N_drought)) {
      Ndrought_byhab[i,] <- colSums(N_drought[,i] * site_hab_props[, 2:7])
      Nnodrought_byhab[i,] <- colSums(N_nodrought[,i] * site_hab_props[, 2:7])
    }
    
    Nnodrought_byhab <- as.data.frame(Nnodrought_byhab) %>% 
      mutate(species = thisspec, i = row_number(), cond_iter = cond_iter,
             type = "nodrought")
    Ndrought_byhab <- as.data.frame(Ndrought_byhab) %>% 
      mutate(species = thisspec, i = row_number(), cond_iter = cond_iter,
             type = "drought")
    
    spec_habcount_df_list[[k]] <- bind_rows(Nnodrought_byhab, Ndrought_byhab)
  }
  
  # How many sites have 0 counts for each draw?
  prop_zero <- apply(cbind(N_drought, N_nodrought), 2, function(x) {mean(x == 0)})
  
  # In each iter, what fraction of sites is responsible for 50% of the count?
  frac90 <- apply(cbind(N_drought, N_nodrought), 2, function(x) {
    x <- x[x > 0]
    1 - (max(which(cumsum(sort(x)) < 0.1 * sum(x))) / length(x))
  })
  frac50 <- apply(cbind(N_drought, N_nodrought), 2, function(x) {
    x <- x[x > 0]
    1 - (max(which(cumsum(sort(x)) < 0.5 * sum(x))) / length(x))
  })
  
  Ndiff_total <- N_drought_total - N_nodrought_total
  Ndiff_rel <- Ndiff_total / N_nodrought_total
  
  diagnostics_df_list[[k]] <- data.frame(
    i = 1:length(frac50),
    prop_zero = prop_zero,
    frac50 = frac50,
    frac90 = frac90,
    species = thisspec,
    med_count = apply(cbind(N_drought, N_nodrought), 2, median),
    q90_count = apply(cbind(N_drought, N_nodrought), 2, quantile, probs = 0.9),
    q95_count = apply(cbind(N_drought, N_nodrought), 2, quantile, probs = 0.95),
    q99_count = apply(cbind(N_drought, N_nodrought), 2, quantile, probs = 0.99),
    q999_count = apply(cbind(N_drought, N_nodrought), 2, quantile, probs = 0.999),
    max_count = apply(cbind(N_drought, N_nodrought), 2, max),
    type = rep(c("drought", "no drought"), each = ncol(N_drought)),
    cond_type = conditions$name[cond_iter]
  )
  
  spec_est_list[[k]] <- data.frame(
    i = 1:length(N_drought_total),
    drought_N = N_drought_total,
    nodrought_N = N_nodrought_total,
    decline = N_drought_total - N_nodrought_total < 0,
    species = thisspec
  )
  
  #### Calculate interaction terms per habitat type #####
  habint_result_list <- list()
  ct <- 0
  for (h_ind in 1:length(land_covars)) {
    # Empirical quantiles of the covariate we're working on
    habgrid <- quantile(unlist(all_model_data[, land_covars[h_ind]]), probs = 1:9/10)

    # Get indices of important betas from data N-mixture model
    data_hab_main_ind <- 
      which(params == paste0("lambda_beta[", 
                             which(lambda_covars == land_covars[h_ind]), 
                             "]"))
    
    # warning("Precipitation is currently left out of this")
    ppt_main_ind <- which(params == "lambda_beta[10]")
    ppt_hab_int_ind <- 
      which(params == paste0("lambda_beta[",
                             which(lambda_covars == paste0("ppt_", land_covars[h_ind])),
                             "]"))
    
    evi_main_ind <- which(params == "lambda_beta[7]")
    evi_hab_int_ind <-
      which(params == paste0("lambda_beta[",
                             which(lambda_covars == paste0("EVI_", land_covars[h_ind])),
                             "]"))
    ndwi_main_ind <- which(params == "lambda_beta[8]")
    ndwi_hab_int_ind <-
      which(params == paste0("lambda_beta[",
                             which(lambda_covars == paste0("NDWI_", land_covars[h_ind])),
                             "]"))
    tmax_main_ind <- which(params == "lambda_beta[8]")
    tmax_hab_int_ind <-
      which(params == paste0("lambda_beta[",
                             which(lambda_covars == paste0("tmax_", land_covars[h_ind])),
                             "]"))
    
    for (i in 1:length(habgrid)) {
      if (!duplicated(habgrid)[i]) {
    
        target_interaction <- 
          # evi component
          (unlist(samples_mtx[, evi_hab_int_ind]) * (
            2 * unlist(evi_post$maineff_draws[post_samps, paste0("b_spei:", land_covars[h_ind])]) * habgrid[i] + 
              unlist(evi_post$maineff_draws[post_samps, "b_spei"])
          ) + unlist(evi_post$maineff_draws[post_samps, paste0("b_spei:", land_covars[h_ind])]) * 
                unlist(samples_mtx[, evi_main_ind])) +
          # ndwi component
          (unlist(samples_mtx[, ndwi_hab_int_ind]) * (
            2 * unlist(ndwi_post$maineff_draws[post_samps, paste0("b_spei:", land_covars[h_ind])]) * habgrid[i] + 
              unlist(ndwi_post$maineff_draws[post_samps, "b_spei"])
          ) + unlist(ndwi_post$maineff_draws[post_samps, paste0("b_spei:", land_covars[h_ind])]) * 
            unlist(samples_mtx[, ndwi_main_ind])) +
          # tmax component
          (unlist(samples_mtx[, tmax_hab_int_ind]) * (
            2 * unlist(tmax_post$maineff_draws[post_samps, paste0("b_spei:", land_covars[h_ind])]) * habgrid[i] + 
              unlist(tmax_post$maineff_draws[post_samps, "b_spei"])
          ) + unlist(tmax_post$maineff_draws[post_samps, paste0("b_spei:", land_covars[h_ind])]) * 
            unlist(samples_mtx[, tmax_main_ind])) +
          # ppt component
          (unlist(samples_mtx[, ppt_hab_int_ind]) * (
            2 * unlist(ppt_post$maineff_draws[post_samps, paste0("b_spei:", land_covars[h_ind])]) * habgrid[i] + 
              unlist(ppt_post$maineff_draws[post_samps, "b_spei"])
          ) + unlist(ppt_post$maineff_draws[post_samps, paste0("b_spei:", land_covars[h_ind])]) * 
            unlist(samples_mtx[, ppt_main_ind]))
          
        
        ct <- ct + 1
        habint_result_list[[ct]] <- data.frame(
          hab = land_covars[h_ind],
          H = unname(habgrid[i]),
          Hq = (1:9/10)[i],
          value = unname(target_interaction),
          iter = 1:length(target_interaction)
        )
      }
    }
    spec_habint_df_list[[k]] <- bind_rows(habint_result_list) %>% 
      mutate(species = thisspec)
    spec_habint_summary_list[[k]] <- spec_habint_df_list[[k]] %>% 
      group_by(species, hab, H, Hq) %>% 
      summarize(
        med = median(value),
        Q025 = quantile(value, probs = 0.025),
        Q975 = quantile(value, probs = 0.975),
        .groups = "drop"
      )
      
  }
  
  pb$tick()
}
  
  # Handle all species results for this condition type
  bind_rows(spec_est_list) %>% 
    mutate(cond_iter = cond_iter, cond_type = conditions$name[cond_iter]) %>% 
    write_csv(paste0("intermediate/posterior_drought_diffs_", cond_iter, ".csv"))
  
  bind_rows(diagnostics_df_list) %>% 
    mutate(cond_iter = cond_iter, cond_type = conditions$name[cond_iter]) %>% 
    write_csv(paste0("intermediate/posterior_diagnostics_", cond_iter, ".csv"))
  
  bind_rows(spec_habint_df_list) %>% 
    mutate(cond_iter = cond_iter, cond_type = conditions$name[cond_iter]) %>% 
    write_csv(paste0("intermediate/spec_habint_draws_", cond_iter, ".csv"))

  # Only do hab counts for main counterfactual scenario
  if (cond_iter == 1) {
    bind_rows(spec_habcount_df_list) %>% 
      mutate(cond_type = conditions$name[cond_iter]) %>% 
      write_csv(paste0("intermediate/spec_habcount_draws_", cond_iter, ".csv"))
    
    bind_rows(spec_habcount_df_list) %>% 
      mutate(cond_type = conditions$name[cond_iter]) %>% 
      write_csv(paste0("intermediate/spec_habcount_draws_", cond_iter, ".csv"))
  }
  
}







