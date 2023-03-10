###############################################################################
# 08_vis_data.R
# Author: Ben Goldstein
# Description: This script creates some summary plots based on input data.
#      
###############################################################################


library(Matrix)
library(tidyverse)
library(ggmap)
library(sf)
library(gridExtra)
library(viridis)
source("code/make_data.R")
source("code/03_modeling_fn.R")

cover_colors <- c("Perennial ag"     = "#1f78b4",
                  "Grassland"        = "#33a02c",
                  "Riparian"         = "#e31a1c",
                  "Row and field ag" = "#fdbf6f",
                  "Other natural"    = "#6a3d9a",
                  "Developed/Other"  = "#cab2d6")

devtools::install_github("UrbanInstitute/urbnmapr")
library(urbnmapr)

#### Spatial plot of the study area and distribution of eBird data ####
valley <- st_read("data/greatvalley_outline.shp")
valley_ll <- st_transform(valley, st_crs("+proj=longlat"))

cali <- get_urbn_map(map = "states", sf = TRUE) %>% 
  filter(state_name == "California")
cali_ll <- st_transform(cali, st_crs("+proj=longlat"))


map_bounds <- unlist(raster::extent(cali_ll) + 3)
coords.map <- get_stamenmap(bbox = map_bounds[c(1,3,2,4)], zoom = 6, 
                            maptype = "terrain-background")
cali_3857 <- st_transform(cali, 3857)
valley_3857 <- st_transform(valley, 3857)

ebird_data_raw <- read_csv("intermediate/ebird_info_for_model.csv")

points <- st_as_sf(count(ebird_data_raw, LONGITUDE, LATITUDE), 
                   coords = c("LONGITUDE", "LATITUDE"),
                   crs = st_crs("+proj=longlat"))


test_map <- ggmap::get_map(location = unname(st_bbox(cali_ll)),
                           source = "stamen")
attr(test_map, "bb")$ll.lat <- st_bbox(cali_3857)["ymin"]
attr(test_map, "bb")$ll.lon <- st_bbox(cali_3857)["xmin"]
attr(test_map, "bb")$ur.lat <- st_bbox(cali_3857)["ymax"]
attr(test_map, "bb")$ur.lon <- st_bbox(cali_3857)["xmax"]

(spaceplot <- 
  ggmap(test_map, extent = "device", legend = "none") +
  geom_sf(data = cali_3857, color = "black", alpha = 0, 
          inherit.aes = FALSE) +
  geom_sf(data = valley_3857, fill = "darkred", 
          color = "darkred", alpha = 0.4, inherit.aes = FALSE) +
  geom_sf(data = points, inherit.aes = FALSE, size = 0.2,
          color = "black"
          ) +
  ggtitle("(A)") +
  theme(panel.border = element_rect(colour = "white", fill=NA, size=5))
)

# Make plot of time-varying components (# checklists, SPEI)
year_counts <- ebird_data_raw %>% 
  group_by(year) %>% 
  summarize(n = n(),
            spei = mean(spei)) %>% 
  mutate(year = year + 2010)

(timeplot <- ggplot(year_counts) +
  geom_col(aes(year, n, fill = -spei)) +
  scale_x_continuous(breaks = seq(2010, 2020, by = 2)) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) +
  xlab("Year") +
  ylab("Number of checklists") +
  scale_fill_gradient2("Drought", mid = "lightgray", high = "darkred", low = "blue") +
  ggtitle("(B)")
)


# Make plot of covariate space
# 

plot_groups <- c("Row and field ag", "Perennial ag", "Grassland", 
                 "Riparian", "Other natural", "Other")
grid_cover <- read_csv("intermediate/landscape_covs.csv") %>% 
  filter(year == 2010, !is.na(Fruits)) %>% 
  mutate(`Row and field ag` = `Row and field crops` + Rice,
         `Perennial ag` = Fruits + Nuts + Vineyards,
         Other = Urban+Other) %>% 
  rename(`Other natural` = "Other_semi_natural",
         Riparian = `Riparian vegetation`)
  
grid_cover_sampled <- grid_cover %>% 
  filter(cellnum %in% all_model_data$cellnum)


cover_sums <- data.frame(
  cover_type = plot_groups,
  prop = colSums(grid_cover[, plot_groups]) / sum(grid_cover[, plot_groups]),
  ncell_nonzero = colSums(grid_cover[, plot_groups] > 0),
  type = "Whole area"
) %>% 
  arrange(-prop)

cover_sums_sampled <- data.frame(
  cover_type = plot_groups,
  prop = colSums(grid_cover_sampled[, plot_groups]) / 
    sum(grid_cover_sampled[, plot_groups]),
  ncell_nonzero = colSums(grid_cover_sampled[, plot_groups] > 0),
  type = "Sampled"
) %>% 
  arrange(-prop) %>% 
  bind_rows(cover_sums) %>% 
  mutate(
    type = factor(type, levels = c("Whole area", "Sampled")),
    cover_type = factor(
      recode(cover_type, "Other_semi_natural" = "Other semi-natural",
             "Other" = "Developed/Other"),
      levels = names(cover_colors))
  )

(coverplot <- cover_sums_sampled %>% 
  ggplot(aes(fill = cover_type, x = type, y = prop)) +
  geom_col(position = "stack", size = 0.2) +
  theme_minimal() +
  ylab("% cover") +
  xlab("") +
  scale_fill_manual("", values  = cover_colors) +
  ggtitle("(C)") +
  theme(legend.key.size = unit(0.4, 'cm'), 
        axis.text.y = element_text(size = 7)
        )
  )



layout_mtx <- matrix(c(
  1,1,1,2,2,2,2,
  1,1,1,2,2,2,2,
  1,1,1,3,3,3,3,
  1,1,1,3,3,3,3
), nrow = 4, byrow = TRUE)

grid.arrange(grobs = list(spaceplot, timeplot, coverplot), 
             layout_matrix  = layout_mtx,
             padding = unit(1, "line")) %>% 
  ggsave(filename = "plots/fig1_data_summary.jpg",
         width = 10, height = 5)



# Columns: species code, common name, number of obs, median nonzero obs

birds_tbl <- read_csv("intermediate/birdcodes_used.csv")[
    , c("name", "alpha.code")
  ] %>% 
  mutate(prop_cls = NA, median_ct = NA)

for (i in 1:nrow(birds_tbl)) {
  birds_tbl$prop_cls[i] <- sum(all_model_data[, birds_tbl$alpha.code[i]] > 0) / nrow(all_model_data)
  birds_tbl$median_ct[i] <- median(
    all_model_data[, birds_tbl$alpha.code[i]][all_model_data[, birds_tbl$alpha.code[i]] > 0]
  )
}


birds_tbl$prop_cls <- sprintf("%.2f", birds_tbl$prop_cls)
birds_tbl$median_ct <- sprintf("%.0f", birds_tbl$median_ct)

write_csv(arrange(birds_tbl, alpha.code), "tables/birds_tbl.csv")


