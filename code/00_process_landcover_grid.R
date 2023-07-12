###############################################################################
# 00_process_covars.R
# Description: This script calculates land cover percentages using
#    raw input land cover classes for every grid cell in the Central
#    Valley.
###############################################################################

# Libraries
library(raster)
library(rgdal)
library(rgeos)
library(sf)
library(lubridate)
library(fasterize)
library(SpaDES)
library(spdep)
library(landscapemetrics)
library(tidyverse)

verbose <- TRUE


# Load a shapefile defining the study region
valley <- readOGR(dsn = "data/greatvalley_outline.shp")
valley_ll <- spTransform(valley, CRS("+proj=longlat"))
valley_pts <- valley_ll@polygons[[1]]@Polygons[[1]]@coords

# Put on a the 30x30 grid
valley_extent <- extent(valley) + 1000
valley_extent[1] <- valley_extent[1] - valley_extent[1] %% 30
valley_extent[2] <- valley_extent[2] - valley_extent[2] %% 30
valley_extent[3] <- valley_extent[3] - valley_extent[3] %% 30
valley_extent[4] <- valley_extent[4] - valley_extent[4] %% 30



#### Process 2018 land cover data ####
# Ultimate goal is to produce the file data/LandIQ_fveg_2018.tif that contains
# land cover from LandIQ masking fveg
# If-statement logic provides that some intermediate calculations may not need
# to be re-run. Processing these rasters is slow so it's nice to do it piecewise
if (!file.exists("data/LandIQ_fveg_2018.tif")) {
  # First, create the fveg raster
  if (!file.exists("data/fveg_clipped.grd")) {
    fveg_raw <- raster("data/fveg_WHRNHUM")
    fveg <- raster::crop(fveg_raw, valley_extent)
    values(fveg)[values(fveg) == 255] <- NA
    raster::writeRaster(fveg, "data/fveg_clipped.grd")
    
    fveg_attributes <- read_csv("data/fveg_attributes.txt")
    
    fveg_attributes %>%
      distinct(WHRNUM, WHRNAME) %>%
      arrange(WHRNUM) %>%
      write_csv("data/fveg_legend.csv")
    fveg_legend <- read_csv("data/fveg_legend.csv")
    
  } else {
    fveg <- raster("../CA_GV_birds/data/fveg_clipped.grd")
    fveg_legend <- read_csv("../CA_GV_birds/data/fveg_legend.csv")
  }
  
  # Put LandIQ on the same grid
  if (!file.exists("data/landIQ_clipped_2018.grd")) {
    landIQ_raw <- read_sf("../CA_GV_birds/data/i15_crop_mapping_2018.shp")
    landIQ <- st_transform(landIQ_raw, crs(fveg))
    
    # Determine unique number ID for each unique cropID combination,
    # write csv of unique categories
    uniqland <- landIQ %>%
      as.data.frame() %>%
      distinct(CROPTYP2) %>%
      arrange(CROPTYP2) %>%
      dplyr::mutate(CropID=100+(1:nrow(.)))
    
    write_csv(uniqland, "intermediate/LandIQ_2018_UniqCrop.csv")
    
    landIQ_cropID <- landIQ %>%
      left_join(uniqland, by=c("CROPTYP2"))
    
    # Initialize empty rasters
    ras <- raster()
    # Set the raster extent
    extent(ras) <- extent(spTransform(valley, crs(fveg)))
    # Set raster resolution (meters)
    res(ras) <- 30
    # Define the CRS
    crs(ras) <- crs(landIQ)
    
    landIQ_ras <- fasterize(landIQ_cropID, ras, field="CropID")
    landIQ <- resample(landIQ_ras, fveg, method = "ngb")
    
    writeRaster(landIQ, "data/landIQ_clipped_2018.grd", overwrite = T)
  } else {
    landIQ <- raster("data/landIQ_clipped_2018.grd")
    uniqland <- read_csv("intermediate/LandIQ_2018_UniqCrop.csv")
  }
  
  landuse <- cover(landIQ, fveg)
  writeRaster(landuse, "data/LandIQ_fveg_2018.tif", overwrite=TRUE)
  
} else {
  landuse <- raster("data/LandIQ_fveg_2018.tif")
}

CropIQ_legend <- read_csv("intermediate/LandIQ_2018_UniqCrop.csv") %>%
  rename(CODE = CROPTYP2) %>%
  left_join(read_csv("data/CropIQ_legend.csv"), by = "CODE") %>%
  rename(class = CropID, cover = NAME) %>%
  select(class, cover) %>%
  mutate(source = "CropIQ")
fveg_legend <- read_csv("data/fveg_legend.csv") %>%
  rename(class = WHRNUM, cover = WHRNAME) %>%
  mutate(source = "FVEG")
cover_legend <- bind_rows(CropIQ_legend, fveg_legend)

if (sum(duplicated(cover_legend$class)) > 0) {
  stop("Duplicate cover codes across sources.")
}


# Make the spatial grid
gridres <- 1000 # grid resolution in m
grid_r <- raster(raster::extent(landuse) + gridres, resolution = gridres)
# grid_r_hi <- raster(extent(landuse) + gridres, resolution = gridres / 20)

projection(grid_r) <- projection(landuse)
values(grid_r) <- 0
writeRaster(grid_r, "intermediate/GV_grid_raster.tif", overwrite = TRUE)

grid_as_poly <- rasterToPolygons(grid_r) %>%
  spTransform(crs(landuse))
grid_as_poly@data$in_valley <- (grid_as_poly %over% spTransform(valley, crs(grid_as_poly)))$EPA_REGION
grid_as_poly <- grid_as_poly[!is.na(grid_as_poly@data$in_valley),]
grid_as_poly@data$cellnum <- 1:nrow(grid_as_poly@data)

polygrid_as_raster <- fasterize(st_as_sf(grid_as_poly), landuse, field = "cellnum")

landuse_cts_by_cell <- data.frame(
  cellnum = values(polygrid_as_raster),
  class  = values(landuse)
) %>% 
  filter(!is.na(cellnum)) %>% 
  count(cellnum, class)

cell_totals <- landuse_cts_by_cell %>% 
  group_by(cellnum) %>% 
  summarize(total_cells = sum(n), .groups = "drop")

landuse_cts_by_cell <- left_join(landuse_cts_by_cell, cell_totals, by = "cellnum")

landuse_cts_by_cell$value <- landuse_cts_by_cell$n / landuse_cts_by_cell$total_cells

all_vals_2018 <- landuse_cts_by_cell %>%
  left_join(cover_legend, by = "class")

# Write out the legend and df of cover class proportions
write_csv(cover_legend, "intermediate/grid_cover_legend_2018.csv")
write_csv(all_vals_2018, "intermediate/grid_cover_2018.csv")
