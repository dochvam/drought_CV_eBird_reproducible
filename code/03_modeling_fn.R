###############################################################################
# 03_modeling_fn.R
# Description: This script defines the functions used in modeling. The 
#      two major functions are write_nmix_model, which uses inputs to create a
#      script defining nimbleModel code for a given model specification; and
#      fit_one_nmix, which conducts MCMC model estimation for a given species
#      and set of input parameters 
#      
###############################################################################

##### Packages #####
library(nimbleEcology)
library(coda)
library(MCMCvis)
library(tidyverse)


##### Objects #####

# This data frame is a holdover from an earlier study design where we were going
# to do model selection to decide whether to add various aspects of model 
# structure. We ultimately found that model selection was computationally
# infeasible and decided to just go with the most conservative model.
nmix_to_test <- data.frame(
  expand.grid(
    type = "Nmixture",
    betabin = c(FALSE, TRUE),
    negbin = c(FALSE, TRUE),
    site_ranef_occ = c(FALSE, TRUE),
    obs_ranef_det = c(FALSE, TRUE))
) %>% 
  mutate(mi = row_number())


##### Functions #####

# Function to write N-mixture model NIMBLE code based on whether it contains
# various pieces of model structure (see note above)
write_nmix_model <- function(outfile,
                             negbin,
                             betabin,
                             obs_ranef_det,
                             site_ranef_occ
) {
  
  if (negbin) {
    theta_prior <- "theta ~ dunif(0.0001, 25)"
    if (betabin) {
      s_prior <- "s ~ dunif(0.0001, 25)"
      N_distribution <- "
        N[i] ~ dnegbin(size = 1 / theta, prob = 1 / (1 + theta * lambda[i]))
      "
      data_distribution <- "
        y[start[i]:end[i]] ~ dBetaBinom(
          N = N[i],
          shape1 = p[start[i]:end[i]] * s,
          shape2 = s - p[start[i]:end[i]] * s
        )
      "
      oneobs_distribution <- "
        y[start[i]] ~ dBetaBinom_One(
          N = N[i],
          shape1 = p[start[i]] * s,
          shape2 = s - p[start[i]] * s
        )
      "
    } else {
      s_prior <- NULL
      N_distribution <- "
        N[i] ~ dnegbin(size = 1 / theta, prob = 1 / (1 + theta * lambda[i]))
      "
      data_distribution <- "
        for (t in start[i]:end[i]) {
          y[t] ~ dbinom(size = N[i], prob = p[t])
        }
      "
      oneobs_distribution <- "
        y[start[i]] ~ dbinom(size = N[i], prob = p[start[i]])
      "
    }
  } else {
    theta_prior <- NULL
    if (betabin) {
      s_prior <- "s ~ dunif(0.0001, 25)"
      N_distribution <- "
        N[i] ~ dpois(lambda[i])
      "
      data_distribution <- "
        y[start[i]:end[i]] ~ dBetaBinom(
          N = N[i],
          shape1 = p[start[i]:end[i]] * s,
          shape2 = s - p[start[i]:end[i]] * s
        )
      "
      oneobs_distribution <- "
        y[start[i]] ~ dBetaBinom_One(
          N = N[i],
          shape1 = p[start[i]] * s,
          shape2 = s - p[start[i]] * s
        )
      "   
    } else {
      s_prior <- NULL
      N_distribution <- "
        N[i] ~ dpois(lambda[i])
      "
      data_distribution <- "
        for (t in start[i]:end[i]) {
          y[t] ~ dbinom(size = N[i], prob = p[t])
        }
      "
      oneobs_distribution <- "
        y[start[i]] ~ dbinom(size = N[i], prob = p[start[i]])
      "
    }
  }
  
  if (site_ranef_occ) {
    site_ranef_prior <- "
    site_ranef_sd ~ dunif(0.001, 10)
    for (i in 1:ncell) {
      site_ranef[i] ~ dnorm(0, sd = site_ranef_sd)
    }
    "
    lambda_formula <- "    
      log(lambda[i]) <- inprod(lambda_covars[start[i], 1:nLamBeta], lambda_beta[1:nLamBeta]) +
                          site_ranef[cellnum[i]]
    "
  } else {
    site_ranef_prior <- NULL
    lambda_formula <- "    
      log(lambda[i]) <- inprod(lambda_covars[start[i], 1:nLamBeta], lambda_beta[1:nLamBeta])
    "
    
  }
  
  if (obs_ranef_det) {
    obs_ranef_prior <- "
    obs_ranef_sd ~ dunif(0.001, 10)
    for (i in 1:nObserver) {
      obs_ranef[i] ~ dnorm(0, sd = obs_ranef_sd)
    }
    "
    p_formula <- "    
      logit(p[i]) <- min(logit(1-1e-6), 
                         max(logit(1e-6),
                         inprod(p_covars[i, 1:nPBeta], p_beta[1:nPBeta]) +
                         obs_ranef[obsID[i]]))
    "
  } else {
    obs_ranef_prior <- NULL
    p_formula <- "    
      logit(p[i]) <- min(logit(1-1e-6), 
                         max(logit(1e-6),
                         inprod(p_covars[i, 1:nPBeta], p_beta[1:nPBeta])))"
  }
  
  
  writeLines(paste0("model_code <- nimbleCode({
  for (i in 1:nLamBeta) {
    lambda_beta[i] ~ dnorm(0, sd = 2.25)
  }
  p_beta[1] ~ dlogis(0, 1)
  for (i in 2:nPBeta) {
    p_beta[i] ~ dnorm(0, sd = 2.25)
  }
  ",
                    s_prior, "\n\t",
                    theta_prior, "\n\t",
                    obs_ranef_prior, "\n\t",
                    site_ranef_prior,
                    "
  for (i in 1:nobs) {
    ",
                    p_formula, 
                    "  }
  
  # For each spatial unit (n > 1)...
  for (i in 1:(ncy - numOneObs)) {
  ",
                    lambda_formula, "\n",
                    N_distribution,
                    data_distribution,
                    "} 
  
    for (i in (ncy - numOneObs + 1):(ncy)) {
      ",
                    lambda_formula,
                    N_distribution,
                    oneobs_distribution,
                    
                    "}
                    
    
  })"),
             con = outfile
  )
  
}


# Main function used for modeling. Runs MCMC estimation for an N-mixture model
# of a given species. Options to use WAIC and to overwrite. Can provide a vector
# of species, in which case the model will iteratively estimate each. This saves
# time because the NIMBLE model only needs to be built and compiled once. 
# The option "rerun" indicates that the function should look for previous MCMC
# samples and run new chains to append on; this lets us continue sampling for a
# model that hasn't converged yet.
fit_one_nmix <- function(mod_row, spec_vec,
                         use_WAIC = TRUE, 
                         overwrite = FALSE,
                         rerun = FALSE) {
  
  in_modeling <- TRUE
  
  if (!is.data.frame(mod_row) || nrow(mod_row) != 1) {
    stop("mod_row must be a one-row data.frame")
  }
  
  # Write the N-mixture model
  tempfile <- paste0("temp/modfit_", mod_row$mi, "_", spec_vec[1], ".R")
  write_nmix_model(outfile = tempfile, 
                   negbin = mod_row$negbin, 
                   betabin = mod_row$betabin, 
                   obs_ranef_det = mod_row$obs_ranef_det, 
                   site_ranef_occ = mod_row$site_ranef_occ
  )
  source(tempfile)
  
  # Load the data
  source("code/make_data.R")
  
  # Only rerun one chain at a time
  if (rerun) {
    nc <- 1
  }
  
  # Nodes to monitor
  monitors <- c("lambda_beta", "p_beta",
                "site_ranef_sd", "site_ranef",
                "obs_ranef_sd", "obs_ranef"
                )
  if (mod_row$betabin) monitors <- c(monitors, "s")
  if (mod_row$negbin) monitors <- c(monitors, "theta")
  if (use_WAIC) {
    monitors <- c(monitors, "N")
  }
  
  #### Modeling ####
  # Get start and end indices
  start_vec <- end_vec <- max_obs <- numeric(max(all_model_data$cyID))
  
  for (i in 1:max(all_model_data$cyID)) {
    this_id <- i
    start_vec[i] <- min(which(all_model_data$cyID == i))
    end_vec[i] <- max(which(all_model_data$cyID == i))
  }
  
  
  for (i in 1:max(all_model_data$cyID)) {
    max_obs[i] <- max(all_model_data[start_vec[i]:end_vec[i], spec_vec[1]])
  }
  
  inits <- function(max_obs) {
    
    return(list(
      lambda_beta = c(1, rnorm(length(lambda_covars) - 1, sd = 0.05)),
      p_beta = rnorm(length(p_covars), sd = 0.2),
      obs_ranef = rnorm(max(all_model_data$obsID), sd = 0.5),
      obs_ranef_sd = 0.5,
      site_ranef = rnorm(length(unique(all_model_data$cellnum)), sd = 0.5),
      site_ranef_sd = 0.5,
      theta = 1,
      s = 1,
      N = max_obs + 1
    ))
  }
  
  thisspec <- spec_vec[1]
  joint_model <- nimbleModel(
    code = model_code,
    constants = list(
      nObserver = max(all_model_data$obsID),
      obsID = all_model_data$obsID,
      ncy = max(all_model_data$cyID),
      cy = all_model_data$cyID,
      ncell = length(unique(all_model_data$cellnum)),
      cellnum = as.numeric(as.factor(all_model_data$cellnum))[start_vec],
      nLamBeta = length(lambda_covars),
      nPBeta = length(p_covars),
      numOneObs = numOneObs,
      start = start_vec,
      end = end_vec,
      nobs = nrow(all_model_data)
    ), data = list(
      y = unlist(all_model_data[, thisspec]),
      lambda_covars = all_model_data[, c(lambda_covars)],
      p_covars = all_model_data[, c(p_covars)]
    ), 
    inits = inits(max_obs)
  )
  
  joint_mcmc_conf <- configureMCMC(joint_model, print = TRUE)
  joint_mcmc_conf$addMonitors(monitors)

  # Intercepts and random effects
  joint_mcmc_conf$addSampler(target = c("s", "obs_ranef_sd", "p_beta[1]"), 
                             type = "AF_slice")
  joint_mcmc_conf$addSampler(target = c("lambda_beta[1]", "site_ranef_sd", "theta"), 
                      type = "AF_slice")

  # Detection int and date
  joint_mcmc_conf$addSampler(target = c("p_beta[1]", "p_beta[5]", "p_beta[6]"), 
                             type = "AF_slice")

  # Major land covars
  joint_mcmc_conf$addSampler(target = joint_model$expandNodeNames("lambda_beta[2:6]"),
                             type = "AF_slice")
  
  # NDWI, EVI, tmax
  joint_mcmc_conf$addSampler(target = joint_model$expandNodeNames("lambda_beta[7:11]"),
                             type = "AF_slice")
  
  # NDWI and interactions
  joint_mcmc_conf$addSampler(target = c(
                              "lambda_beta[8]",
                              joint_model$expandNodeNames("lambda_beta[12:16]")),
                      type = "AF_slice")
  
  # EVI and interactions
  joint_mcmc_conf$addSampler(target = c(
        "lambda_beta[7]",
        joint_model$expandNodeNames("lambda_beta[17:21]")),
        type = "AF_slice"
      )
  
  # Tmax and interactions
  joint_mcmc_conf$addSampler(target = c(
    "lambda_beta[9]",
    joint_model$expandNodeNames("lambda_beta[22:26]")),
    type = "AF_slice"
  )
  
  # ppt and interactions
  joint_mcmc_conf$addSampler(target = c(
    "lambda_beta[10]",
    joint_model$expandNodeNames("lambda_beta[27:31]")),
    type = "AF_slice"
  )
  
  # Year effects
  joint_mcmc_conf$addSampler(target = c("lambda_beta[1]", 
                                        joint_model$expandNodeNames("lambda_beta[32:40]")), 
                      type = "AF_slice")
  
  # Drop the samplers I want to overwrite with slice samplers
  sampler_list <- unlist(lapply(joint_mcmc_conf$getSamplers(), 
                                function(x) paste(x$target, collapse = ", ")))
  exclude <- which(
    sampler_list == "s" | sampler_list == "obs_ranef_sd" | sampler_list == "p_beta[1]" |
      sampler_list == "p_beta[5]" | sampler_list == "p_beta[6]" | 
      sampler_list == "site_ranef_sd" | sampler_list == "lambda_beta[6]" | sampler_list == "lambda_beta[1]" |
      sampler_list == "theta" | 
      sampler_list %in% paste0("lambda_beta[", 1:40, "]")
  )
  duplicate <- which(grepl("obs_ranef_sd,", sampler_list) |
                     grepl("site_ranef_sd,", sampler_list))
  order <- (1:length(sampler_list))[!1:length(sampler_list) %in% exclude]
  
  
  joint_mcmc_conf$setSamplerExecutionOrder(
    order
  )
  joint_mcmc <- buildMCMC(joint_mcmc_conf)
  
  complist <- compileNimble(joint_model, joint_mcmc)
  
  this_lambda_covars_df <-
    data.frame(
      name = lambda_covars,
      param = paste0("lambda_beta[", 1:length(lambda_covars), "]"),
      parnum = 1:length(lambda_covars),
      layer = "lambda"
    )
  this_p_covars_df <-
    data.frame(
      name = p_covars,
      param = paste0("p_beta[", 1:length(p_covars), "]"),
      parnum = 1:length(p_covars),
      layer = "p"
    )
  this_covars_df <- bind_rows(this_p_covars_df, this_lambda_covars_df)
  
  joint_summary <- list()
  
  for (k in 1:length(spec_vec)) {
    thisspec_start_time <- Sys.time()
    output_file <- paste0("intermediate/SSAMs/MCMC_SSAM_", spec_vec[k], "_", str_pad(mod_row$mi, 2, pad = "0"), ".RDS")
    if (rerun) {
      output_file_final <- paste0(
        substr(output_file, 1, nchar(output_file) - 4),
        "_rerun",
        substr(output_file, nchar(output_file) - 3, nchar(output_file))
      )
      
      if (file.exists(output_file_final)) {
        last_result_samples <- readRDS(output_file_final)$samples
      } else {
        last_result_samples <- readRDS(output_file)$samples
      }
      
    } else {
      output_file_final <- output_file
    }
    
    # Do we need to run samples?
    if (rerun || overwrite || !file.exists(output_file)) {

      # Print a message to console
      writeLines(paste0(spec_vec[k], " - ", mod_row$mi, " - ", format(Sys.time()- 7*3600, "%b %d %X")))
      
      # Update the data
      complist$joint_model$setData(list(
        y = as.numeric(unlist(all_model_data[, spec_vec[[k]]]))
      ))
      
      # Update the "max_obs" value to ensure sensible inits for this species
      for (i in 1:max(all_model_data$cyID)) {
        max_obs[i] <- max(all_model_data[start_vec[i]:end_vec[i], spec_vec[k]])
      }
      
      # Start running MCMC!
      MCMC_output <- runMCMC(complist$joint_mcmc, niter = ni, nburnin = nb,
                             thin = nt, nchains = nc,
                             inits = inits(max_obs)) #started at 11:53
      samples <- MCMC_output
      
      
      if (use_WAIC) {
        WAIC <- calculateWAIC(do.call(rbind, MCMC_output), complist$joint_model)
      } else {
        WAIC <- NA
      }
      
      # Process MCMC samples and produce a usable summary file
      samples_list <- list()
      if (nc > 1) {
        for (i in 1:nc) {
          uniques <- unlist(apply(samples[[i]], 2, function(x) length(unique(x))))
          
          cols_to_keep <- which(uniques > 2 | !grepl("logProb", colnames(samples[[i]])))
          
          samples_list[[i]] <- mcmc(data = samples[[i]][, cols_to_keep], start = nb + 1, end = ni, thin = nt)
          # drop NA log prob columns?
          
        }
      } else {
        samples_list[[1]] <- mcmc(data = samples, start = nb + 1, end = ni, thin = nt)
      }
      
      if (rerun) {
        joint_samples <- do.call(mcmc.list, c(mcmc.list(samples_list), last_result_samples))
      } else {
        joint_samples <- mcmc.list(samples_list)
      }        
      
      
      joint_summary <- MCMCsummary(joint_samples, probs = c(0.05, 0.5, 0.95),
                                        HPD = TRUE, hpd_prob = 0.95) %>%
        mutate(nonzero = sign(`95%_HPDL`) == sign(`95%_HPDU`) & sign(`95%_HPDL`) != 0)
      
      joint_summary$param <- rownames(joint_summary)
      joint_summary$index <- 1:nrow(joint_summary)
      joint_summary$species <- spec_vec[k]
      joint_summary$layer <- lapply(joint_summary$param, function(x) {
        if (grepl("lambda_", x)) {
          "lambda"
        } else if (grepl("p_", x)) {
          "p"
        } else {
          NA
        }
      }) %>% unlist()
      
      joint_summary$parname <- get_parname(joint_summary, this_covars_df)
      
      # Save output
      thisspec_end_time <- Sys.time()
      saveRDS(list(
        samples = joint_samples,
        summary = joint_summary,
        WAIC = WAIC,
        time_taken = thisspec_start_time - thisspec_end_time
      ), file = output_file_final)
    } else {
      writeLines(paste0("Skipping ", spec_vec[k]))
    }
  }
  gc()
}




