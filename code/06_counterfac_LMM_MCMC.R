###############################################################################
# 06_counterfac_LMM_MCMC.R
# Author: Ben Goldstein
# Description: This script prepares data for and runs four Bayesian linear
#      mixed models of environmental variables. It uses prediction methodology
#      to generate posterior predictive distributions of each environmental
#      variable, in each grid cell, under drought and non-drought conditions.
###############################################################################


library(raster)
library(sf)
library(Matrix)
library(tidyverse)
library(brms)

# Should we redo the creation of data or just retrieve cached df?
redo_data <- FALSE
ni <- 15000; nb <- 5000; nc <- 4

# Load eBird spatial grid objects
ebird_data_grid <- raster("intermediate/GV_grid_raster.tif")
sites_w_data <- read_csv("intermediate/ebird_info_for_model.csv") %>% 
  distinct(cellnum) %>% 
  .$cellnum

#### DATA PREP ####
### read in SPEI for all cells / years (intermediate/SPEI_df.csv)
SPEI_df <- read_csv("intermediate/SPEI_df.csv") %>% 
  dplyr::select(-cover_year) %>% 
  filter(cellnum %in% sites_w_data)

# Create the NDWI dataset for all cells
if (!file.exists("intermediate/NDWI_df.csv") | redo_data) {
  NDWI_files <- list.files("data/GV_NDWI", pattern = "Valley_ndwi.*.tif",
                          full.names = TRUE)
  
  NDWI_stack <- stack(NDWI_files)
  names(NDWI_stack) <- paste0("Y", substr(NDWI_files, 26, 29))
  grid_as_pts <- rasterToPoints(ebird_data_grid)
  
  grid_as_poly <- st_read("intermediate/grid_as_poly.shp")
  
  grid_as_pts <- st_as_sf(as.data.frame(grid_as_pts[, 1:2]),
                          crs = crs(grid_as_poly), 
                          coords = c("x", "y"))
  
  ndwi_cellnum <- st_join(grid_as_pts, grid_as_poly, join = st_within)
  
  NDWI_df <- raster::values(NDWI_stack) %>% 
    as.data.frame() %>% 
    mutate(cellnum = ndwi_cellnum$cellnum) %>% 
    filter(!is.na(cellnum), !is.na(Y2010),
           cellnum %in% sites_w_data) %>% 
    pivot_longer(cols = paste0("Y", 2010:2019)) %>%
    rename(year = name, NDWI = value) %>% 
    filter(!is.na(NDWI)) %>% 
    mutate(year = as.numeric(substr(year, 2, 5))) %>% 
    group_by(cellnum) %>% 
    mutate(NDWI_dev = NDWI - mean(NDWI)) %>% 
    ungroup()
  
  write_csv(NDWI_df, "intermediate/NDWI_df.csv")
} else {
  NDWI_df <- read_csv("intermediate/NDWI_df.csv")
}

# Create the EVI dataset for all cells
if (!file.exists("intermediate/EVI_df.csv") | redo_data) {
  EVI_files <- list.files("data/GV_EVI", pattern = "Valley_EVI.*.tif",
                          full.names = TRUE)
  
  EVI_stack <- stack(EVI_files)
  names(EVI_stack) <- paste0("Y", substr(EVI_files, 24, 27))
  grid_as_pts <- rasterToPoints(ebird_data_grid)
  
  grid_as_poly <- st_read("intermediate/grid_as_poly.shp")
  
  grid_as_pts <- st_as_sf(as.data.frame(grid_as_pts[, 1:2]),
                          crs = crs(grid_as_poly), 
                          coords = c("x", "y"))
  
  EVI_cellnum <- st_join(grid_as_pts, grid_as_poly, join = st_within)
  
  EVI_df <- raster::values(EVI_stack) %>% 
    as.data.frame() %>% 
    mutate(cellnum = EVI_cellnum$cellnum) %>% 
    filter(!is.na(cellnum), !is.na(Y2010),
           cellnum %in% sites_w_data) %>% 
    pivot_longer(cols = paste0("Y", 2010:2019)) %>%
    rename(year = name, EVI = value) %>% 
    filter(!is.na(EVI)) %>% 
    mutate(year = as.numeric(substr(year, 2, 5))) %>% 
    group_by(cellnum) %>% 
    mutate(EVI_dev = EVI - mean(EVI)) %>% 
    ungroup()
  
  write_csv(EVI_df, "intermediate/EVI_df.csv")
} else {
  EVI_df <- read_csv("intermediate/EVI_df.csv")
}

if (redo_data) warning(
  "Just a heads up, redo_data will NOT redo tmax or precip; these are in script 02"
)
if (!file.exists("intermediate/precip_bycellyr.csv") |
    !file.exists("intermediate/tmax_yearly.csv")
    ) {
  stop("You need to run script 02 to get the precip and/or tmax files")
} else {
  # Load in precipitation and tmax data by cell
  precip_df <- read_csv("intermediate/precip_bycellyr.csv") %>% 
    filter(cellnum %in% sites_w_data)
  tmax_df   <- read_csv("intermediate/tmax_yearly.csv") %>% 
    filter(cellnum %in% sites_w_data)
}

covar_data <- reduce(
  list(NDWI_df, EVI_df, tmax_df, precip_df, SPEI_df), 
  full_join,
  by = c("cellnum", "year")
)


#### MODELING ####
# Names of habitat types
covgrps <- c("row_field_ag", 
             "perennial_ag", 
             "riparian", 
             "grass_pasture", 
             "other_natural")

# Get habitat types df and reformat
landscape_covs <- read_csv("intermediate/landscape_covs.csv") %>% 
  filter(cellnum %in% sites_w_data) %>% 
  group_by(cellnum) %>% 
  filter(year == max(year)) %>% 
  ungroup() %>% 
  # mutate(Latitude = as.numeric(scale(Latitude))) %>% 
  mutate(
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
    Latitude = scale(Latitude)
  ) %>% 
  dplyr::select(cellnum, all_of(covgrps), Latitude) %>% 
  distinct()
    

all_data <- left_join(covar_data, landscape_covs, by = "cellnum") %>% 
  mutate(cellID_factor = as.factor(cellnum)) %>% 
  na.omit()

covars <- c("Latitude", paste0("spei*", covgrps))


cell_info <- all_data %>% 
  dplyr::select(cellnum, cellID_factor, Latitude,
                all_of(covgrps)) %>% 
  distinct()
spei_archetype_vals <- data.frame(
  drought_level = c("high", "med", "low"),
  spei = c(-1.41, 0, 1.66)
)

# Make data frames from which to predict posterior distributions
newdata_drought <- cell_info %>% 
  mutate(spei = spei_archetype_vals$spei[
    spei_archetype_vals$drought_level == "high"]
  )
newdata_nodrought <- cell_info %>% 
  mutate(spei = spei_archetype_vals$spei[
    spei_archetype_vals$drought_level == "low"]
  )

# Start the actual model-fitting workflow
target_responses <- c(
  "NDWI", "EVI", "ppt", "tmax_ssn"
)

# Loop over the 4 variables, fit LMMs and get posterior predictions
for (i in 1:length(target_responses)) {
  resp <- target_responses[[i]]
  
  f <- as.formula(paste(
    resp, " ~ ", paste(
      covars, collapse = "+"
    ), "+ (1 | cellID_factor)", collapse = ""
  ))
  
  this_brms <- 
    brms::brm(data = all_data, formula = f, 
              iter = ni, warmup = nb, chains = nc, cores = nc,
              prior = set_prior("cauchy(0, 2)", class = "sd") +
                set_prior("cauchy(0, 2)", class = "sd", group = "cellID_factor") +
                set_prior("normal(0, 5)", class = "b")
    )
  
  # What do we want to do with brms output?
  # - check empirical coverage probability
  # - predict NDWI values at new data points:
  
  # Retrieve posterior draws of model coefficients
  draws <- as_draws_df(this_brms)
  draws <- draws[, 1:(min(which(grepl("sd_", colnames(draws)))) - 1)]
  
  
  # Generate posterior predictions, our original goal
  posterior_predicted_drought <-
    posterior_predict(this_brms, newdata = newdata_drought)
  posterior_predicted_nodrought <- 
    posterior_predict(this_brms, newdata = newdata_nodrought)
  
  # Evaluate ECP
  incl <- function(range, val) {
    return(val >= range[1] && val <= range[2])
  }
  evaluate_ECP <- function(input_dat, response,
                           posterior_samples, 
                           progress = TRUE) {
    # Prep output vector
    in_95CI_all <- logical(nrow(input_dat))
    
    if (progress) {
      pb <- progress::progress_bar$new(total = nrow(input_dat))
    }
    
    # Do quantiles row-by-row over all data
    for (f in 1:nrow(input_dat)) {
      if (progress) pb$tick()
      ci95 <- quantile(posterior_samples[, f], probs = c(0.025, 0.975))
      in_95CI_all[f] <- incl(ci95, input_dat[f, response])
    }
    
    return(mean(in_95CI_all))
  }
  
  posterior_predicted_forECP <- posterior_predict(this_brms)
  
  ECP <- evaluate_ECP(input_dat = all_data, response = resp,
                      posterior_samples = posterior_predicted_forECP)
  
  saveRDS(
    list(
      mod_summary = summary(this_brms), 
      maineff_draws = draws,
      ppred_drought = posterior_predicted_drought,
      ppred_nodrought = posterior_predicted_nodrought,
      data_drought = newdata_drought,
      data_nodrought = newdata_nodrought,
      ECP = ECP
    ),
    file = paste0("intermediate/", resp, "_LMM_posterior.RDS")
  )
  
}
