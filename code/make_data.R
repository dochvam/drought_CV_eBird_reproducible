###############################################################################
# make_data.R
# Author: Ben Goldstein
# Description: This helper script is used in a number of places in the workflow.
#      Its purpose is to produce the single data frame "all_model_data" that
#      represents the data as modeled for all species, where each row is an
#      eBird checklist with columns for all covariates and for counts of all
#      species.
###############################################################################


# Create helper fn to standardize creating all_data, cdfw_data
library(nimble)
library(Matrix)
library(tidyverse)

# MCMC parameters
ni <- 15000
nb <- 5000
nc <- 3
nt <- 10
use_WAIC_global <- FALSE


ni_brms <- 2000
nb_brms <- 500
nc_brms <- 3
nt_brms <- 2

land_covars  <- c("row_field_ag", "perennial_ag", "riparian", "grass_pasture", 
                  "other_natural")

land_covars_unscaled <- c(paste0(land_covars, "_unscaled"), "urban_unscaled",
                          "other_unscaled")

lambda_covars <- c("intercept_abd", land_covars, 
                   "EVI", "NDWI", "tmax_ssn", "ppt",
                   "NDWI_tmax", 
                   paste("NDWI", land_covars, sep = "_"),
                   paste("EVI", land_covars, sep = "_"),
                   paste("tmax", land_covars, sep = "_"),
                   paste("ppt", land_covars, sep = "_"),
                   paste0("year", 1:9))

p_covars <- c("intercept_det", "duration", "ebd_tod", "ebd_tod_squared",
              "julian_ebd", "julian_ebd_squared", "group_size", "tmax_daily")

cell_vals <- read_csv("intermediate/landscape_covs.csv") %>% 
  filter(year == 2018)

cell_lats <- cell_vals %>% 
  distinct(cellnum, Latitude)

EVI_df <- read_csv("intermediate/EVI_df.csv") %>% 
  mutate(year = year - 2010)
NDWI_df <- read_csv("intermediate/NDWI_df.csv") %>% 
  mutate(year = year - 2010)

birdcodes <- read_csv("intermediate/birdcodes_used.csv")
birdcodes$SCIENTIFIC.NAME[birdcodes$name_clean == "Double-crested_Cormorant"] <-
  "Nannopterum auritum"

species_vector <- birdcodes$alpha.code

ebird_data_raw <- read_csv("intermediate/ebird_info_for_model.csv")



ebird_data_intermediate <- ebird_data_raw %>%
  left_join(cell_lats) %>% 
  select(-LONGITUDE, -LATITUDE) %>% 
  mutate(intercept_abd = 1, intercept_det = 1) %>%
  left_join(EVI_df, by = c("cellnum", "year")) %>% 
  left_join(NDWI_df, by = c("cellnum", "year")) %>% 
  mutate(
    obsID = as.numeric(as.factor(OBSERVER.ID)),
    
    julian_ebd = scale(julian_ebd),
    ebd_tod = scale(ebd_tod),
    duration = scale(duration),
    group_size = scale(group_size),
    intercept_abd = 1, intercept_det_ebd = 1,
    julian_ebd_squared = julian_ebd^2,
    ebd_tod_squared = ebd_tod^2,
    Latitude = scale(Latitude),
    tmax_daily = scale(tmax_daily),
    
    year1 = as.numeric(year == 1),
    year2 = as.numeric(year == 2),
    year3 = as.numeric(year == 3),
    year4 = as.numeric(year == 4),
    year5 = as.numeric(year == 5),
    year6 = as.numeric(year == 6),
    year7 = as.numeric(year == 7),
    year8 = as.numeric(year == 8),
    year9 = as.numeric(year == 9),

    row_field_ag_unscaled = `Row and field crops` + Rice,
    perennial_ag_unscaled = Fruits + Nuts + Vineyards,
    riparian_unscaled = `Riparian vegetation`,
    grass_pasture_unscaled = Grassland,
    other_natural_unscaled = Other_semi_natural,
    
    row_field_ag = scale(row_field_ag_unscaled),
    perennial_ag = scale(perennial_ag_unscaled),
    riparian = scale(riparian_unscaled),
    grass_pasture = scale(grass_pasture_unscaled),
    other_natural = scale(other_natural_unscaled),
    
    NDWI = scale(NDWI),
    EVI = scale(EVI),
    tmax_ssn = scale(tmax_ssn),
    ppt = scale(ppt)
  ) %>%
  mutate(
    EVI_tmax = EVI * tmax_ssn,
    NDWI_tmax = NDWI * tmax_ssn,
    
    NDWI_row_field_ag = NDWI * row_field_ag,
    NDWI_perennial_ag = NDWI * perennial_ag,
    NDWI_riparian = NDWI * riparian,
    NDWI_grass_pasture = NDWI * grass_pasture,
    NDWI_other_natural = NDWI * other_natural,
    
    EVI_row_field_ag = EVI * row_field_ag,
    EVI_perennial_ag = EVI * perennial_ag,
    EVI_riparian = EVI * riparian,
    EVI_grass_pasture = EVI * grass_pasture,
    EVI_other_natural = EVI * other_natural,
    
    tmax_row_field_ag = tmax_ssn * row_field_ag,
    tmax_perennial_ag = tmax_ssn * perennial_ag,
    tmax_riparian = tmax_ssn * riparian,
    tmax_grass_pasture = tmax_ssn * grass_pasture,
    tmax_other_natural = tmax_ssn * other_natural,

    ppt_row_field_ag = ppt * row_field_ag,
    ppt_perennial_ag = ppt * perennial_ag,
    ppt_riparian = ppt * riparian,
    ppt_grass_pasture = ppt * grass_pasture,
    ppt_other_natural = ppt * other_natural
  ) %>% 
  filter(!is.na(cellnum), !is.na(tmax_other_natural), !is.na(grass_pasture_unscaled)) %>% 
  filter(!is.na(EVI), !is.na(NDWI))





land_scaling_factors <- data.frame(
  mean = c(
    attr(ebird_data_intermediate$row_field_ag, "scaled:center"),
    attr(ebird_data_intermediate$perennial_ag, "scaled:center"),
    attr(ebird_data_intermediate$riparian, "scaled:center"),
    attr(ebird_data_intermediate$grass_pasture, "scaled:center"),
    attr(ebird_data_intermediate$other_natural, "scaled:center")
  ),
  sd = c(
    attr(ebird_data_intermediate$row_field_ag, "scaled:scale"),
    attr(ebird_data_intermediate$perennial_ag, "scaled:scale"),
    attr(ebird_data_intermediate$riparian, "scaled:scale"),
    attr(ebird_data_intermediate$grass_pasture, "scaled:scale"),
    attr(ebird_data_intermediate$other_natural, "scaled:scale")
  ),
  param = c(
    "row_field_ag",
    "perennial_ag",
    "riparian",
    "grass_pasture",
    "other_natural"
  )
) %>%
  mutate(zero_pct = (0 - mean) / sd,
         ninety_pct = (0.9 - mean) / sd,
         onehundred_pct = (1 - mean) / sd)

sites <- count(ebird_data_intermediate, cellnum, year) %>%
  arrange(-n) %>%
  mutate(cyID = row_number())

all_model_data <- left_join(ebird_data_intermediate, sites, by = c("cellnum", "year")) %>%
  arrange(cyID)
numOneObs <- sum(sites$n == 1)
all_model_data[is.na(all_model_data)] <- 0

get_specnum <- function(x) {
  if (!is.character(x)) stop("Provide a character string or vector.")
  rtn <- strsplit(x, split = ",") %>%
    lapply(function(x) x[length(x)]) %>%
    unlist() %>%
    gsub(pattern = "[^0-9.-]", replacement = "") %>%
    as.numeric()
  
  rtn[grepl("beta_mu", samples$summary$param) |
        grepl("beta_sd", samples$summary$param) |
        grepl("det_sy_sd", samples$summary$param)] <- NA
  return(rtn)
}

get_parname <- function(df, covars_df) {
  parname <- character(nrow(df))
  
  if ("parnum" %in% colnames(df)) {
    parnum_vec <- df$parnum
  } else {
    parnum_vec <- numeric(nrow(df))
    for (i in 1:nrow(df)) {
      parnum_vec[i] <- as.numeric(gsub(strsplit(df$param[i], split = ",")[[1]][1],
                                       pattern = "[^0-9.-]", replacement = ""))
    }
  }
  
  if (!"layer" %in% colnames(df)) {
    df$layer <- ifelse(grepl("lambda_", df$param), "lambda",
                       ifelse(grepl("p_", df$param), "p", NA))
  }
  
  for (i in 1:nrow(df)) {
    if (!is.na(parnum_vec[i]) && !is.na(df$layer[i])) {
      parname[i] <- covars_df$name[covars_df$layer == df$layer[i] & covars_df$parnum == parnum_vec[i]]
    }
  }
  
  parname
}



rescale <- function(main_df, target_cols, target_df) {
  
  for (i in 1:length(target_cols)) {
    this_ctr <- attr(target_df[, target_cols[i]][[1]], "scaled:center")
    this_scl <- attr(target_df[, target_cols[i]][[1]], "scaled:scale")
    
    if (is.null(this_ctr) || is.null(this_scl)) {
      stop("Column ", target_cols[i], " is not a scaled vector")
    }
    
    main_df[, target_cols[i]] <- (main_df[, target_cols[i]] - this_ctr) / this_scl
  }
  
  main_df
}



calc_psi_from_samples_onespec <- function(samples, species, modnum) {
  source("code/fit_onespec_models/onespec_fn.R")
  
  if (modnum %in% c(5, 6, 9, 10)) {
    stop("Haven't built this for occ ranef models yet")
  }
  
  mod_data_list <- get_mod_data(modnum, all_model_data)
  mod_data <- mod_data_list$mod_data
  numOneObs <- mod_data_list$numOneObs
  
  covars_list <- get_covars(modnum)
  p_covars <- covars_list$p_covars
  psi_covars <- covars_list$psi_covars
  
  site_data <- mod_data %>%
    select(cyID, cellnum, year, all_of(psi_covars)) %>%
    distinct()
  
  samples <- do.call(rbind, samples)
  params <- colnames(samples)
  
  psi_cols <- which(grepl("psi_beta\\[", params))
  
  lpsi <- samples[, psi_cols] %*% t(as.matrix(site_data[, psi_covars]))
  
  psi <- expit(lpsi)
  
  # Why 83% CI? https://rstudio-pubs-static.s3.amazonaws.com/132971_a902bb2b962b407e9e9436559c6f5d36.html
  data.frame(
    year = site_data$year,
    cellnum = site_data$cellnum,
    species = species,
    model = modnum,
    psi_mu = colMeans(psi),
    psi_085 = apply(psi, 2, quantile, probs = 0.085),
    psi_915 = apply(psi, 2, quantile, probs = 0.915)
  )
}

if (exists("in_modeling") && in_modeling) {
  rm(cell_lats)
  rm(cell_vals)
  rm(ebird_data_raw)
  rm(ebird_data_intermediate)
  rm(EVI_df)
  rm(NDWI_df)
}




