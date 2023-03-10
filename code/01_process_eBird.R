###############################################################################
# 01_process_eBird.R
# Author: Ben Goldstein
# Description: This script processes eBird data. Checklists are filtered for
#     quality and aggregated to the Central Valley 1km grid.
###############################################################################

##### eBird filtering parameters #####

complete_only <- TRUE
stationary_only <- TRUE
year_range <- 2010:2019
date_range <- 91:181
max_duration_mins <- 180
max_distance_km <- 1 # ignored if stationary_only == TRUE
max_num_observers <- 8

response_type <- "count" # "detection"
drop_cls_w_X <- response_type == "count"


##### Libraries #####

library(rgdal)
library(sp)
library(rgeos)
library(raster)
library(tidyverse)
library(prism)
library(lubridate)
library(taxalight)
library(nimbleEcology)


if (!response_type %in% c("count", "detection")) {
  stop("Variable response_type must be either 'count' or 'detection'.")
}

#### 1. Load raw data ####

# eBird data:
ebird_CA_obs <- 
  read_csv("data/CA_Aug2020_species_counts.csv")
ebird_CA_checklists <- 
  read_csv("data/CA_Aug2020_checklist_info.csv")
ebird_protocol_codes <- read_csv("data/ebird_protocol_codes.csv")

# Drop checklists that recorded any "X" if we're doing count data
if (drop_cls_w_X) {
  cls_w_X <- ebird_CA_obs %>% 
    filter(total_count == "X") %>% 
    .$SAMPLING.EVENT.IDENTIFIER %>% 
    unique()
  
  ebird_CA_checklists <- ebird_CA_checklists %>% 
    filter(!SAMPLING.EVENT.IDENTIFIER %in% cls_w_X)
  ebird_CA_obs <- ebird_CA_obs %>% 
    filter(!SAMPLING.EVENT.IDENTIFIER %in% cls_w_X)
}


# Read in CA Great Valley ecoregion
valley <- readOGR(dsn = "data/greatvalley_outline.shp")
valley_ll <- spTransform(valley, CRS("+proj=longlat"))
valley_pts <- valley_ll@polygons[[1]]@Polygons[[1]]@coords

#### 2. Process eBird ####
# Filter eBird data to:
#  - GV ecoregion
#  - Complete checklist
#  - 2016-2017
#  - Match date range of ARU surveys
#  - Stationary point count or traveling point count protocols
# & track how much data is lost at each step

# Range of autonomous recorders: Roughly March 25 - July 06
# Filter by time, completeness
checklists_filtered <- ebird_CA_checklists %>% 
  filter(ALL.SPECIES.REPORTED == 1 || !complete_only) %>% 
  filter(year(OBSERVATION.DATE) %in% year_range) %>% 
  filter(yday(OBSERVATION.DATE) %in% date_range)


# Filter by space (inside GV)
unique_locations <- checklists_filtered %>% 
  dplyr::select(LONGITUDE, LATITUDE) %>% 
  distinct()

unique_locations$inside_area <- point.in.polygon(point.x = unique_locations$LONGITUDE,
                                                 point.y = unique_locations$LATITUDE,
                                                 pol.x = valley_pts[,1], 
                                                 pol.y = valley_pts[,2])

checklists_in_area <- left_join(checklists_filtered, unique_locations) %>% 
  filter(inside_area == 1)

# Filter by checklist protocol
acceptable_checklists <- checklists_in_area %>% 
  filter(PROTOCOL.CODE == "P21" | (PROTOCOL.CODE == "P22" && !stationary_only)) %>% 
  mutate(year = year(OBSERVATION.DATE))

# Filter maximum distance traveled to 1km
acceptable_checklists <- acceptable_checklists %>% 
  mutate(EFFORT.DISTANCE.KM = ifelse(is.na(EFFORT.DISTANCE.KM), 0, EFFORT.DISTANCE.KM)) %>% 
  filter(EFFORT.DISTANCE.KM <= max_distance_km,
         DURATION.MINUTES <= max_duration_mins,
         NUMBER.OBSERVERS <= max_num_observers)


#### 4. Write end product: a list of eBird checklists 
# Get obs
area_obs <- ebird_CA_obs %>% 
  filter(SAMPLING.EVENT.IDENTIFIER %in% 
           acceptable_checklists$SAMPLING.EVENT.IDENTIFIER)
rm(ebird_CA_checklists)
rm(ebird_CA_obs)
# rwb_obs <- gv_obs %>% filter(name_clean == "Red-winged_Blackbird")

# Get unique locations
area_locations <- acceptable_checklists %>% 
  count(LONGITUDE, LATITUDE, year) %>% 
  filter(n > 1) %>% 
  dplyr::select(-n) %>% 
  mutate(site_id = row_number())

ebird_species <- area_obs %>% 
  distinct(name_clean, SCIENTIFIC.NAME) %>% 
  filter(!grepl("sp\\.", SCIENTIFIC.NAME)) %>% 
  filter(!grepl("\\/", SCIENTIFIC.NAME)) %>% 
  filter(!grepl("\\(", SCIENTIFIC.NAME)) %>% 
  filter(!grepl("\\)", SCIENTIFIC.NAME)) %>% 
  filter(!grepl(" x ", SCIENTIFIC.NAME)) 


# Write out vaerious products
# List of species:
write_csv(ebird_species, "intermediate/ebird_species.csv")
# List of cell IDs
write_csv(area_locations, "intermediate/area_locations.csv")
# The eBird species observations (1 row per species per checklist)
write_csv(area_obs, "intermediate/area_obs.csv")
# The checklists w associated metadata
write_csv(acceptable_checklists, "intermediate/area_checklists.csv")


