###############################################################################
# 02_covariate_data.R
# Description: This script prepares all covariate data, including: aggregating
#     land cover types into categories used for modeling, retrieving and
#     preparing temperature and precipitation data, processing EVI and NDWI.
#     Then covariate data is joined to eBird data to make the file used
#     in modeling.
###############################################################################


library(rgdal)
library(rgeos)
library(sp)
library(sf)
library(raster)
library(tidyverse)
library(prism)
library(lubridate)
library(taxadb)


#### 1. Read in raw data ####

sepFruitVineyard <- TRUE

min_cls <- 500

grid <- readOGR("intermediate/grid_as_poly.shp")
grid_ras <- raster("intermediate/GV_grid_raster.tif")

grid_poly_centroids <- bind_rows(lapply(grid@polygons, function(x) {
  as.data.frame(t(colMeans(x@Polygons[[1]]@coords[1:4,])))
})) %>% 
  mutate(cellnum = as.numeric(grid$cellnum)) %>% 
  rename(x = V1, y = V2)
grid_poly_cent_pts <- SpatialPoints(coords = grid_poly_centroids[,c("x", "y")], 
                           proj4string = crs(grid))

# eBird data:
ebird_cls <- read_csv("intermediate/area_checklists.csv")
ebird_obs <- read_csv("intermediate/area_obs.csv")

# Birdcodes
breeding_birds <- read_csv("data/all_CA_breeding_birds.csv") %>% 
  left_join(count(ebird_obs, name_clean, SCIENTIFIC.NAME), by = "name_clean") %>% 
  filter(n >= min_cls)
write_csv(breeding_birds, "intermediate/birdcodes_used.csv")

ebird_obs <- ebird_obs %>% filter(name_clean %in% breeding_birds$name_clean)

# Read in CA Great Valley ecoregion shapefile
valley <- readOGR(dsn = "data/greatvalley_outline.shp") %>% 
  spTransform(proj4string(grid))
valley_ll <- spTransform(valley, CRS("+proj=longlat"))
valley_pts <- valley_ll@polygons[[1]]@Polygons[[1]]@coords



# Read in PLS grid for 2018
grid_cover <- read_csv("intermediate/grid_cover_2018.csv")

cov_categories <- readxl::read_xlsx("data/LandUseClassification_AcrossSites_FVEG.xlsx") %>% 
  select(RemoteSensingProduct, Class, SubClass, `Classification Scheme 2`, Class3)

if (sepFruitVineyard) {
  cov_categories$`Classification Scheme 2`[
    cov_categories$`Classification Scheme 2` == "Fruits_Vineyards"
  ] <- 
    cov_categories$Class3[
      cov_categories$`Classification Scheme 2` == "Fruits_Vineyards"
    ]
}

grid_cover_legend_2018 <- read_csv("intermediate/grid_cover_legend_2018.csv")

cov_categories$SubClass[grepl("New", cov_categories$SubClass)] <-
  grid_cover_legend_2018$cover[grepl("New", grid_cover_legend_2018$cover)]
grid_cover_legend_2018$got_it <- grid_cover_legend_2018$cover %in% cov_categories$SubClass
grid_cover_legend_2018$in_dat <- grid_cover_legend_2018$cover %in% grid_cover$cover
main_legend <- left_join(grid_cover_legend_2018, cov_categories,
                         by = c("cover" = "SubClass")) %>% 
  select(class, cover, `Classification Scheme 2`) %>% 
  rename(covgrp = `Classification Scheme 2`) %>% 
  distinct(cover, covgrp)

grid_covariate_long <- left_join(grid_cover, main_legend, by = c("cover")) %>% 
  # group_by(CO_MTRS, covgrp) %>% 
  group_by(cellnum, covgrp) %>% 
  filter(!is.na(covgrp)) %>%
  summarize(value = sum(value, na.rm = T), .groups = "drop")

grid_covariate_wide_2018 <- pivot_wider(grid_covariate_long, 
                                        names_from = covgrp,
                                        values_from = value, 
                                        values_fill = 0) %>% 
  mutate(cover_year = 2018)



grid_covariate_wide <- grid_covariate_wide_2018

#### Process PRISM data ####
# Get temperature data
if (!file.exists("intermediate/tmax_bycelldate.csv")) {
  prism_sampling_dates <- expand.grid(2010:2019, 4:6, 1:31) %>% 
    filter(!(Var2 %in% c(4, 6) & Var3 == 31)) %>% 
    apply(1, function(x) paste(x, collapse = "-")) %>% 
    as.Date()
  
  prism_sampling_dates <- c(prism_sampling_dates, 
                            as.Date(paste0(2010:2019, "-03-31")))
  
  prism_set_dl_dir("data/prism/")
  get_prism_dailys(
    type = "tmax", keepZip = TRUE,
    dates = prism_sampling_dates)
  
  prism_files <- prism_archive_ls()
  tmax_ras <- pd_stack(prism_files[grepl("_tmax_", prism_files)])
  
  cell_tmax_list <- list()
  
  pb <- progress::progress_bar$new(total = nlayers(tmax_ras))
  for (i in 1:nlayers(tmax_ras)) {
    pb$tick()
    this_tmax_resampled <- raster::projectRaster(from = tmax_ras[[i]], 
                                                 to = grid_ras, 
                                                 method = "ngb")
    
    cell_tmax_list[[i]] <- data.frame(
      cellnum = grid_poly_centroids$cellnum,
      value = raster::extract(this_tmax_resampled,
                              grid_poly_cent_pts),
      date = lubridate::as_date(substr(names(this_tmax_resampled), 25, 32), 
                                format = "%Y%m%d")
    )
  }
  
  write_csv(bind_rows(cell_tmax_list), "intermediate/tmax_bycelldate.csv")
}

# Get precip data
if (!file.exists("intermediate/precip_bycellyr.csv")) {
  prism_set_dl_dir("data/prism/")
  
  prism_files <- prism_archive_ls()
  tmax_ras <- pd_stack(prism_files[grepl("_tmax_", prism_files)])
  
  get_prism_monthlys(
    type = "ppt", keepZip = TRUE, years = 2009:2019,
    mon = 1:12
  )
  
  prism_files <- prism_archive_ls()
  ppt_ras <- pd_stack(prism_files[grepl("_ppt_", prism_files)])
  
  cell_ppt_list <- list()
  
  pb <- progress::progress_bar$new(total = nlayers(ppt_ras))
  for (i in 1:nlayers(ppt_ras)) {
    pb$tick()
    this_ppt_resampled <- raster::projectRaster(from = ppt_ras[[i]], 
                                                 to = grid_ras, 
                                                 method = "ngb")
    
    cell_ppt_list[[i]] <- data.frame(
      cellnum = grid_poly_centroids$cellnum,
      value = raster::extract(this_ppt_resampled,
                              grid_poly_cent_pts),
      year  = as.numeric(substr(names(this_ppt_resampled), 24, 27)),
      month = as.numeric(substr(names(this_ppt_resampled), 28, 29))
    )
  }
  
  cell_ppt_df <- bind_rows(cell_ppt_list) %>% 
    mutate(ag_year = ifelse(month <= 6, year, year + 1)) %>% 
    group_by(cellnum, ag_year) %>% 
    summarize(ppt = sum(value)) %>% 
    rename(year = ag_year) %>% 
    filter(year %in% 2010:2019)
  
  
  write_csv(cell_ppt_df, "intermediate/precip_bycellyr.csv")
}

cell_ppt_df <- read_csv("intermediate/precip_bycellyr.csv")


# Get elevation data
r <- raster::getData("worldclim", var = "alt", res = 0.5, path = "data", 
                     lon = -125, lat = 45)

r2 <- raster::getData("worldclim", var = "alt", res = 0.5, path = "data",
                      lon = -115, lat = 45)

r_elev <- raster::mosaic(r, r2, fun = mean)


grid_proj <- spTransform(grid, crs(r_elev))



landcov_cols <- c(colnames(grid_covariate_wide)[2:ncol(grid_covariate_wide)], "spei")

sitelocs <- as.data.frame(geosphere::centroid(grid))
colnames(sitelocs) <- c("Longitude", "Latitude")
sitelocs$cellnum <- grid$cellnum




# Get drought data
spei <- raster::stack("data/spei07.nc")
spei <- raster::subset(spei, c(paste0("X", 2010:2019, ".04.01")))

spei_vals <- data.frame(
  cellnum = rep(unique(grid$cellnum), length(2010:2019)),
  year = rep(2010:2019, each = length(unique(grid$cellnum))),
  spei = NA
)

grid_centers <- sp::coordinates(grid_proj)

for (i in 1:length(2010:2019)) {
  cat(i, "\t")
  spei_vals$spei[spei_vals$year == (2010:2019)[i]] <-
    unlist(raster::extract(spei[[i]], grid_centers, fun = mean))
}

spei_vals <- spei_vals %>% 
  mutate(cover_year = ifelse(
    year >= 2018, 2018, 2018
  ))
write_csv(spei_vals, "intermediate/SPEI_df.csv")


allsite_covs <- left_join(full_join(
  grid_covariate_wide, spei_vals,
  by = c("cellnum", "cover_year")
), sitelocs, by = "cellnum"
)

# Covariates across the whole landscape
write_csv(allsite_covs, "intermediate/landscape_covs.csv")

#### 2. Associate checklists with their grid cells ####

ebird_locs <- ebird_cls %>% distinct(LONGITUDE, LATITUDE, year)
ebird_pts <- SpatialPoints(coords = ebird_locs[,c("LONGITUDE", "LATITUDE")], 
                           proj4string = CRS("+proj=longlat")) %>% 
  spTransform(crs(grid))

# ebird_locs$CO_MTRS <- sp::over(ebird_pts, plsnet_valley)$CO_MTRS
ebird_locs$cellnum <- sp::over(ebird_pts, grid)$cellnum
ebird_cls_wcovs <- left_join(left_join(ebird_cls, ebird_locs), allsite_covs)
# ebird_cls_wcovs <- left_join(ebird_cls, ebird_locs) %>% 
#   left_join(pesticide_dat, by = c("CO_MTRS" = "comtrs", "year"))
# ebird_cls_wcovs <- left_join(ebird_cls, ebird_locs) %>% 
#   left_join(spei_vals, by = c("cellnum", "year"))

grid_ebd_summary <- ebird_cls_wcovs %>% 
  count(cellnum) %>% 
  rename(n_ebd = n)

#### Deal w tmax #### 

if (!(file.exists("intermediate/tmax_yearly.csv") &&
      file.exists("intermediate/tmax_daily_obsDates.csv"))) {
  
  tmax_cd_df <- read_csv("intermediate/tmax_bycelldate.csv")
  
  # Tmax average per year per cell
  tmax_peryear <- tmax_cd_df %>% 
    mutate(year =  lubridate::year(date)) %>% 
    group_by(cellnum, year) %>% 
    summarize(tmax_ssn = mean(value))
  write_csv(ungroup(tmax_peryear), "intermediate/tmax_yearly.csv")
  
  
  # Tmax only in cell/days w/ data
  ebird_sitedates <- distinct(ebird_cls_wcovs, cellnum, OBSERVATION.DATE)
  ebird_sitedate_str <- 
    paste0(ebird_sitedates$cellnum, "_", ebird_sitedates$OBSERVATION.DATE)
  
  tmax_daily_obsDates <- tmax_cd_df %>% 
    filter(cellnum %in% ebird_sitedates$cellnum) %>% 
    filter(date %in% ebird_sitedates$OBSERVATION.DATE) %>% 
    filter(paste0(cellnum, "_", date) %in% ebird_sitedate_str) %>% 
    rename(tmax_daily = value)
  
  nrow(distinct(ebird_cls_wcovs, cellnum, OBSERVATION.DATE)) ==
    nrow(tmax_daily_obsDates)
  
  write_csv(tmax_daily_obsDates, "intermediate/tmax_daily_obsDates.csv")
  
}

tmax_peryear <- read_csv("intermediate/tmax_yearly.csv")
tmax_daily <- read_csv("intermediate/tmax_daily_obsDates.csv")

#### 4. Final processing ####
ebird_info <- ebird_cls_wcovs %>% 
  left_join(sitelocs, by = "cellnum") %>% 
  left_join(tmax_peryear, by = c("cellnum", "year")) %>% 
  left_join(tmax_daily, by = c("cellnum", "OBSERVATION.DATE" = "date")) %>% 
  left_join(cell_ppt_df, by = c("cellnum", "year")) %>% 
  mutate(duration = DURATION.MINUTES, ebd_tod = as.numeric(TIME.OBSERVATIONS.STARTED) / 60,
         julian_ebd = lubridate::yday(OBSERVATION.DATE), group_size = NUMBER.OBSERVERS,
         distance = EFFORT.DISTANCE.KM, protocol = as.numeric(PROTOCOL.CODE == "P22")) %>% 
  select(-DURATION.MINUTES, -TIME.OBSERVATIONS.STARTED, -OBSERVATION.DATE,
         -NUMBER.OBSERVERS, -EFFORT.AREA.HA, -EFFORT.DISTANCE.KM,
         -ALL.SPECIES.REPORTED) %>% 
  mutate(year = year - 2010)

# Add observations to eBird checklist data, one column per modeled species
for (spec in breeding_birds$alpha.code) {
  ebird_name <- breeding_birds$name_clean[breeding_birds$alpha.code == spec]
  
  this_spec_obs <- ebird_obs %>% 
    filter(name_clean == ebird_name) %>% 
    select(SAMPLING.EVENT.IDENTIFIER, total_count)
  colnames(this_spec_obs)[2] <- spec
  
  ebird_info <- left_join(ebird_info, this_spec_obs, 
                          by = "SAMPLING.EVENT.IDENTIFIER")
  
  ebird_info[is.na(ebird_info[, spec]), spec] <- 0
}


#### 5. Write output ####

# eBird checklists w/ covariates
write_csv(ebird_info, "intermediate/ebird_info_for_model.csv")

