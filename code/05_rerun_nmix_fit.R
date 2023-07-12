###############################################################################
# 05_rerun_nmix_fit.R
# Description: This script is a modified version of 04_main_nmix_fit.R designed
#      to handle species for which the original model run gave insufficient
#      samples. It identifies minimum effective sample size for all species and
#      then runs additional chains as needed.
###############################################################################


library(nimbleEcology)
library(Matrix)
library(tidyverse)
library(parallel)
library(MGMM)
library(mvtnorm)


source("code/abundance/03_modeling_fn.R")
ncores <- 4


source("code/abundance/make_data.R")

all_result_files <-
  list.files("intermediate/SSAMs/", pattern = "MCMC_SSAM", full.names = TRUE)
target_files <- c()
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


processed_list <- lapply(
  target_files, function(x) {
    readRDS(x)
  }
)

summary_df <- bind_rows(lapply(processed_list, function(x) x$summary))
species_res_list <- unlist(lapply(processed_list, function(x) x$summary$species[1]))

# Minimum ESS per species, stochastic nodes only
min_ess <- summary_df %>% 
  filter(str_detect(param, "lambda_beta") |
           str_detect(param, "p_beta") |
           param %in% c("s", "theta", "site_ranef_sd", "sy_ranef_sd")) %>% 
  group_by(species) %>% 
  filter(n.eff == min(n.eff)) %>% 
  select(species, n.eff, param)

# Maximum r-hat per species, stochastic nodes only
max_rhat <- summary_df %>% 
  filter(str_detect(param, "lambda_beta") |
           str_detect(param, "p_beta") |
           param %in% c("s", "theta", "site_ranef_sd", "sy_ranef_sd")) %>% 
  group_by(species) %>% 
  filter(Rhat == max(Rhat)) %>% 
  select(species, Rhat, param) %>% 
  arrange(-Rhat)

specs_to_rerun <- data.frame(
  species = unique(c(
    min_ess$species[min_ess$n.eff < 100]
  )),
  date = Sys.Date()
)

nmix_to_test_rows <- 16

source("code/abundance/make_data.R")

target_specs <- specs_to_rerun$species

if (!dir.exists("intermediate/SSAMs")) {
  dir.create("intermediate/SSAMs")
}

if (!dir.exists("temp")) dir.create("temp")
ts_list <- split(target_specs, cut(seq_along(target_specs), ncores, labels = FALSE))

tempfiles <- character(nrow(nmix_to_test) * length(ts_list))

completed <- list.files("intermediate/SSAMs")

needed <- c()

ct <- 0
# for (row in nmix_to_test_rows) {
for (i in 1:length(ts_list)) {
    ct <- ct + 1
    tempfiles[ct] <- paste0("temp/runscript_", ct, ".R")
    
    specs_as_char <- paste0(ts_list[[i]], collapse = '", "')
    
    writeLines(paste0(
      '
      library(tidyverse)
      library(coda)
      library(MCMCvis)
      library(lubridate)
      library(nimbleEcology)
  
      source("code/abundance/03_modeling_fn.R")
  
      spec_vec <- c("', specs_as_char,'")
      fit_one_nmix(mod_row = nmix_to_test[16,], 
                   spec_vec = spec_vec, 
                   overwrite = TRUE, 
                   use_WAIC = FALSE,
                   rerun = TRUE)
      '), con = tempfiles[ct])
}
# }

tempfiles <- tempfiles[1:ct]

ncores <- min(ncores, length(tempfiles))

cl <- parallel::makeCluster(ncores)
parLapply(cl = cl, tempfiles, source)
stopCluster(cl)




