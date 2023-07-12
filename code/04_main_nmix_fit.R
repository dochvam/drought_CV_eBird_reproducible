###############################################################################
# 04_main_nmix_fit.R
# Description: This script calls helper script 03_modeling_fn.R and executes
#      MCMC estimation of the model for all species. It contains some logic
#      to help the user avoid re-running models for species that have already
#      run. It also uses parLapply to fit models for multiple species in 
#      parallel.
###############################################################################

library(nimbleEcology)
library(Matrix)
library(tidyverse)
library(parallel)

source("code/03_modeling_fn.R")
ncores <- 4

# nmix_to_test_rows <- 1:16
nmix_to_test_rows <- 16

source("code/abundance/make_data.R")
fit_species <- list.files("intermediate/SSAMs/", pattern = "MCMC_SSAM_.*16") %>% 
  substr(11, 14)

overwrite <- FALSE

if (!overwrite) {
  target_specs <- sample(species_vector[!species_vector %in% fit_species])
} else {
  target_specs <- sample(species_vector)
}

if (!dir.exists("intermediate/SSAMs")) {
  dir.create("intermediate/SSAMs")
}

if (!dir.exists("temp")) dir.create("temp")
ts_list <- split(target_specs, cut(seq_along(target_specs), ncores, labels = FALSE))

tempfiles <- character(nrow(nmix_to_test) * length(ts_list))

completed <- list.files("intermediate/SSAMs")

needed <- c()

ct <- 0
# Write out a set of files that can be executed in parallel. Each file will run
# MCMC estimation on a vector of species in sequence
for (i in 1:length(ts_list)) {
  any_targets <- TRUE
  if (!overwrite) {
    targets <- 
      paste0("MCMC_SSAM_", unlist(ts_list[i]), "_", str_pad(nmix_to_test$mi[nmix_to_test_rows], 2, pad = "0"), ".RDS")
    if (all(targets %in% completed)) {
      any_targets <- FALSE
    } else {
      needed <- c(needed, targets[!targets %in% completed])
    }
  }
  
  if (any_targets) {
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
                 overwrite = ', overwrite, ', 
                 use_WAIC = FALSE)
    '), con = tempfiles[ct])
  }
}

tempfiles <- tempfiles[1:ct]

ncores <- min(ncores, length(tempfiles))

# Execute the files in parallel
cl <- parallel::makeCluster(ncores)
parLapply(cl = cl, tempfiles, source)
stopCluster(cl)




