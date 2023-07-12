###############################################################################
# 09_visualize_postpred.R
# Description: This script creates plots based on posterior prediction output.
###############################################################################

library(Matrix)
library(tidyverse)
library(MGMM)
library(gridExtra)

# Extract legend from a ggplot object
get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

source("code/make_data.R")

cred_colors <- c(
  "Credible decline" = "#FD8700",
  "No credible diff." = "#858585",
  "Credible increase" = "#444E7E"
)
cred_alphas <- c(
  "Credible decline"  = 1,
  "No credible diff." = 0.6,
  "Credible increase" = 1
)

posterior_pred <- do.call(bind_rows, 
                          lapply(list.files("intermediate/", pattern = "posterior_drought_diffs", 
                                            full.names = TRUE), read_csv))

#### Are species declining? ####
diffs <- posterior_pred %>% 
  filter(cond_type == "All") %>% 
  mutate(drought_diff = drought_N - nodrought_N,
         drought_diff_pct = drought_diff / nodrought_N
  ) %>% 
  group_by(species, cond_type) %>%
  summarize(med_diff = median(drought_diff), 
            med_diff_pct = median(drought_diff_pct),
            diff_pct_025 = quantile(drought_diff_pct, probs = 0.025),
            diff_pct_975 = quantile(drought_diff_pct, probs = 0.975),
            decl_rate = mean(decline)) %>% 
  mutate(cred_type = ifelse(decl_rate >= 0.975, "Credible decline",
                            ifelse(decl_rate <= 0.025, "Credible increase", 
                                   "No credible diff."))) %>% 
  arrange(med_diff_pct)


(cred_hist <- diffs %>% 
    mutate(cred_type = factor(cred_type, levels = names(cred_colors))) %>% 
    ggplot() +
    geom_bar(aes(group = cred_type, fill = cred_type, cred_type),
             show.legend = FALSE) +
    scale_fill_manual("", values = cred_colors)+
    ylab("Number of species") +
    xlab("") +
    theme_minimal() +
    ggtitle("(a)")
)

(count_mag_plot_A <- diffs %>% 
    filter(diff_pct_975 < 10) %>% 
    mutate(#diff_pct_975 = ifelse(diff_pct_975 > 5, 5, diff_pct_975),
      species = factor(species, levels = diffs$species)) %>% 
    ungroup() %>% 
    mutate(rn = row_number(),
           diff_pct_025 = diff_pct_025 * 100,
           diff_pct_975 = diff_pct_975 * 100,
           med_diff_pct = med_diff_pct * 100
    ) %>% 
    filter(rn < 29) %>% 
    ggplot(aes(species, med_diff_pct,
               ymin = diff_pct_025, ymax = diff_pct_975,
               col = cred_type, alpha = cred_type)) +
    geom_hline(color = "black", size = 1, yintercept = 0) +
    geom_pointrange() +
    theme_minimal() +
    scale_color_manual("", values = cred_colors) +
    scale_alpha_manual("", values = cred_alphas) +
    theme(panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank()) +
    xlab("") + ylab("") +
    coord_flip() +
    theme(legend.position = "none")
)

(count_mag_plot_B <- diffs %>% 
    filter(diff_pct_975 < 10) %>% 
    mutate(#diff_pct_975 = ifelse(diff_pct_975 > 5, 5, diff_pct_975),
      species = factor(species, levels = diffs$species)) %>% 
    ungroup() %>% 
    filter(species != "NSHO", species != "PUMA", species != "RUHU") %>% 
    mutate(rn = row_number(),
           diff_pct_025 = diff_pct_025 * 100,
           diff_pct_975 = diff_pct_975 * 100,
           med_diff_pct = med_diff_pct * 100
    ) %>% 
    filter(rn >= 29) %>% 
    ggplot(aes(species, med_diff_pct,
               ymin = diff_pct_025, ymax = diff_pct_975,
               col = cred_type, alpha = cred_type)) +
    geom_hline(color = "black", size = 1, yintercept = 0) +
    geom_pointrange() +
    theme_minimal() +
    scale_color_manual("", values = cred_colors) +
    scale_alpha_manual("", values = cred_alphas) +
    theme(panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank()) +
    xlab("Species") + ylab("") +
    coord_flip()
)

cmp_legend <- get_legend(count_mag_plot_B)
count_mag_plot_B <- count_mag_plot_B + theme(legend.position = "none")

ggsave("plots/decline_hists.jpg", cred_hist, width = 4, height = 2)
ggsave("plots/fig_s1_count_mags_all.jpg", 
       grid.arrange(count_mag_plot_B, count_mag_plot_A, right = cmp_legend,
                    nrow = 1, bottom = "Difference in abundance during drought (%)"), 
       width = 8, height = 6)


cred_specs <- diffs %>% 
  filter(cred_type != "No credible diff.") %>% 
  left_join(birdcodes, by = c("species" = "alpha.code")) %>% 
  arrange(med_diff_pct)

(cred_spec_mag_plot <- cred_specs %>% 
    mutate(name = factor(name, levels = cred_specs$name)) %>% 
    mutate(diff_pct_025 = diff_pct_025 * 100,
           diff_pct_975 = diff_pct_975 * 100,
           med_diff_pct = med_diff_pct * 100) %>% 
    ggplot(aes(name, med_diff_pct,
               ymin = diff_pct_025, ymax = diff_pct_975,
               col = cred_type, alpha = cred_type)) +
    geom_hline(color = "black", size = 1, yintercept = 0) +
    geom_hline(color = "gray", size = 1, yintercept = -1) +
    geom_hline(color = "gray", size = 1, yintercept = 1) +
    geom_pointrange() +
    theme_minimal() +
    scale_color_manual("", values = cred_colors) +
    scale_alpha_manual("", values = cred_alphas) +
    theme(panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank()) +
    coord_flip() + 
    xlab("") + ylab("Change in abundance during drought (%)") +
    ggtitle("(b)")
  # annotate(geom="text", x=22, y=-0.5, label="(a)",size=8,family="serif")
)

ggsave(cred_spec_mag_plot, filename = "plots/cred_specs.jpg",
       width = 6, height = 5)

ggsave(
  grid.arrange(cred_hist, cred_spec_mag_plot, nrow = 1, widths = c(0.6, 1.2)), 
  filename = "plots/fig2_overall_shifts.jpg", width = 12, height = 5
)

#### Mechanism plots ####

diffs_by_mech <- posterior_pred %>% 
  mutate(drought_diff = drought_N - nodrought_N,
         drought_diff_pct = drought_diff / nodrought_N
  ) %>% 
  group_by(species, cond_type) %>%
  summarize(med_diff = median(drought_diff), 
            med_diff_pct = median(drought_diff_pct),
            diff_pct_025 = quantile(drought_diff_pct, probs = 0.025),
            diff_pct_975 = quantile(drought_diff_pct, probs = 0.975),
            decl_rate = mean(decline)) %>% 
  mutate(cred_type = ifelse(decl_rate >= 0.975, "Credible decline",
                            ifelse(decl_rate <= 0.025, "Credible increase", "No credible diff."))) %>% 
  arrange(med_diff_pct)


mech_summary <- diffs_by_mech %>% 
  select(species, cond_type, cred_type) %>% 
  pivot_wider(names_from = cond_type, values_from = cred_type)

(cred_hist_by_mech <- diffs_by_mech %>% 
    mutate(cred_type = factor(cred_type, levels = names(cred_colors))) %>% 
    filter(cond_type != "All") %>% 
    ggplot() +
    geom_bar(aes(group = cred_type, fill = cred_type, cred_type)) +
    scale_fill_manual("", values = cred_colors) +
    ylab("Number of species") +
    xlab("") +
    theme_minimal() +
    theme(axis.text.x = element_blank()) +
    facet_wrap(~cond_type, nrow = 1)
)

ggsave(cred_hist_by_mech, width = 8, height = 3, 
       filename = "plots/fig_s4_cred_hist_by_mech.jpg")


#### Traits plots ####

# Traits vs. overall shifts
traits <- readxl::read_xlsx("data/AVONET2_eBird.xlsx", 
                            sheet = "AVONET2_eBird")

birdcodes <- left_join(birdcodes, traits[, c("Species2", "Order2",
                                             "Mass", "Habitat", "Migration",
                                             "Trophic.Niche")],
                       by = c("SCIENTIFIC.NAME" = "Species2"))

(cred_hist_niche <- diffs %>% 
    left_join(birdcodes, by = c("species" = "alpha.code")) %>% 
    mutate(cred_type = factor(cred_type, levels = names(cred_colors))) %>% 
    mutate(Trophic.Niche = recode(Trophic.Niche,
                                  "Frugivore" = "Other",
                                  "Nectarivore" = "Other",
                                  "Scavenger" = "Other",
                                  "Vertivore" = "Other",
                                  "Herbivore terrestrial" = "Other",
                                  "Herbivore aquatic" = "Other"
    )) %>% 
    ggplot() +
    geom_bar(aes(group = cred_type, fill = cred_type, cred_type)) +
    scale_fill_manual("", values = cred_colors)+
    ylab("Number of species") +
    xlab("") +
    theme_minimal() + 
    theme(axis.text.x = element_blank()) +
    facet_wrap(~Trophic.Niche, scales = "free_y") +
    ggtitle("(a)")
)

trait_legend <- get_legend(cred_hist_niche)
cred_hist_niche <- cred_hist_niche + theme(legend.position = "none")

(cred_hist_order <- diffs %>% 
    left_join(birdcodes, by = c("species" = "alpha.code")) %>% 
    mutate(cred_type = factor(cred_type, levels = names(cred_colors))) %>% 
    mutate(Order2 = ifelse(Order2 == "Passeriformes", "Passerine", "Non-passerine")) %>% 
    ggplot() +
    geom_bar(aes(group = cred_type, fill = cred_type, cred_type)) +
    scale_fill_manual("", values = cred_colors)+
    ylab("Number of species") +
    xlab("") +
    theme_minimal() + 
    theme(axis.text.x = element_blank()) +
    facet_wrap(~Order2) +
    ggtitle("(b)") +
    theme(legend.position = "none")
)

(count_mag_mass <- diffs %>% 
    left_join(birdcodes, by = c("species" = "alpha.code")) %>% 
    mutate(diff_pct_975 = ifelse(diff_pct_975 > 6, 6, diff_pct_975),
           species = factor(species, levels = diffs$species)) %>% 
    ggplot(aes(log(Mass), med_diff_pct,
               ymin = diff_pct_025, ymax = diff_pct_975,
               col = cred_type, alpha = cred_type)) +
    geom_hline(color = "black", size = 1, yintercept = 0) +
    geom_pointrange() +
    theme_minimal() +
    scale_color_manual("", values = cred_colors) +
    scale_alpha_manual("", values = cred_alphas) +
    theme(panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank()) +
    xlab("Log mass (g)") + ylab("Pct. change in abundance during drought") +
    ggtitle("(c)") +
    theme(legend.position = "none")
)

(cred_hist_migr <- diffs %>% 
    left_join(birdcodes, by = c("species" = "alpha.code")) %>% 
    mutate(cred_type = factor(cred_type, levels = names(cred_colors))) %>% 
    mutate(Migration = factor(recode(Migration, 
                                     "1" = "Sedentary",
                                     "2" = "Partially migratory",
                                     "3" = "Migratory"),
                              levels = c("Sedentary", "Partially migratory",
                                         "Migratory"))) %>% 
    ggplot() +
    geom_bar(aes(group = cred_type, fill = cred_type, cred_type)) +
    scale_fill_manual("", values = cred_colors)+
    ylab("Number of species") +
    xlab("") +
    theme_minimal() + 
    theme(axis.text.x = element_blank()) +
    facet_wrap(~Migration) +
    ggtitle("(d)") +
    theme(legend.position = "none")
)


ggsave(filename = "plots/fig_s2_trait_count.jpg", width = 12, height = 7,
       grid.arrange(cred_hist_niche, cred_hist_order, count_mag_mass,
                    cred_hist_migr, right = trait_legend))





#### Interactions ####
spec_habint_summary <- bind_rows(lapply(
  list.files("intermediate",
             pattern = "spec_habint_draws_1.csv", full.names = TRUE),
  read_csv
)) %>% 
  group_by(cond_type, species, hab, H, Hq) %>% 
  summarize(
    med = median(value),
    Q025 = quantile(value, probs = 0.025),
    Q975 = quantile(value, probs = 0.975),
    .groups = "drop"
  ) %>% 
  mutate(nonzero = sign(Q025) == sign(Q975),
         cred_type = ifelse(nonzero & sign(Q025) == 1, "Credible decrease",
                            ifelse(nonzero & sign(Q025) == -1, "Credible increase", "No credible diff."))
         
  )


species_interactions <- spec_habint_summary %>% 
  group_by(species) %>% 
  summarize(has_interaction = sum(nonzero) > 0,
            num_interactions = sum(nonzero))

hab_interactions <- spec_habint_summary %>% 
  filter(Hq == 0.9) %>% 
  count(hab, cred_type)


spec_w_interactions <- species_interactions$species[species_interactions$has_interaction]

hab_counts <- read_csv("intermediate/spec_habcount_draws_1.csv")



#### Traits vs. hab shifts ####
cred_colors_shift <- cred_colors[1:2]
names(cred_colors_shift) <- c("Habitat shift", "No habitat shift")

(cred_hist_niche <- diffs %>% 
    left_join(birdcodes, by = c("species" = "alpha.code")) %>% 
    mutate(hab_shifts = ifelse(species %in% spec_w_interactions,
                               "Habitat shift", "No habitat shift")) %>% 
    mutate(Trophic.Niche = recode(Trophic.Niche,
                                  "Frugivore" = "Other",
                                  "Nectarivore" = "Other",
                                  "Scavenger" = "Other",
                                  "Vertivore" = "Other",
                                  "Herbivore terrestrial" = "Other",
                                  "Herbivore aquatic" = "Other"
    )) %>% 
    ggplot() +
    geom_bar(aes(group = hab_shifts, fill = hab_shifts, hab_shifts)) +
    scale_fill_manual("", values = cred_colors_shift)+
    ylab("Number of species") +
    xlab("") +
    theme_minimal() + 
    theme(axis.text.x = element_blank()) +
    facet_wrap(~Trophic.Niche, scales = "free_y") +
    ggtitle("(a)")
)

trait_legend <- get_legend(cred_hist_niche)
cred_hist_niche <- cred_hist_niche + theme(legend.position = "none")

(cred_hist_order <- diffs %>% 
    left_join(birdcodes, by = c("species" = "alpha.code")) %>% 
    mutate(hab_shifts = ifelse(species %in% spec_w_interactions,
                               "Habitat shift", "No habitat shift"))  %>% 
    mutate(Order2 = ifelse(Order2 == "Passeriformes", "Passerine", "Non-passerine")) %>% 
    ggplot() +
    geom_bar(aes(group = hab_shifts, fill = hab_shifts, hab_shifts)) +
    scale_fill_manual("", values = cred_colors_shift)+
    ylab("Number of species") +
    xlab("") +
    theme_minimal() + 
    theme(axis.text.x = element_blank()) +
    facet_wrap(~Order2) +
    ggtitle("(b)") +
    theme(legend.position = "none")
)

(count_mag_mass <- diffs %>% 
    left_join(birdcodes, by = c("species" = "alpha.code")) %>% 
    mutate(hab_shifts = ifelse(species %in% spec_w_interactions,
                               "Habitat shift", "No habitat shift")) %>% 
    mutate(diff_pct_975 = ifelse(diff_pct_975 > 6, 6, diff_pct_975),
           species = factor(species, levels = diffs$species)) %>% 
    ggplot(aes(log(Mass), med_diff_pct,
               ymin = diff_pct_025, ymax = diff_pct_975,
               col = hab_shifts, alpha = hab_shifts)) +
    geom_hline(color = "black", size = 1, yintercept = 0) +
    geom_pointrange() +
    theme_minimal() +
    scale_color_manual("", values = cred_colors_shift) +
    scale_alpha_manual("", values = unname(cred_alphas[1:2])) +
    theme(panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank()) +
    xlab("Log mass (g)") + ylab("Pct. change in abundance during drought") +
    ggtitle("(c)") +
    theme(legend.position = "none")
)

(cred_hist_migr <- diffs %>% 
    left_join(birdcodes, by = c("species" = "alpha.code")) %>% 
    mutate(hab_shifts = ifelse(species %in% spec_w_interactions,
                               "Habitat shift", "No habitat shift"))  %>% 
    mutate(Migration = factor(recode(Migration, 
                                     "1" = "Sedentary",
                                     "2" = "Partially migratory",
                                     "3" = "Migratory"),
                              levels = c("Sedentary", "Partially migratory",
                                         "Migratory"))) %>% 
    ggplot() +
    geom_bar(aes(group = hab_shifts, fill = hab_shifts, hab_shifts)) +
    scale_fill_manual("", values = cred_colors_shift)+
    ylab("Number of species") +
    xlab("") +
    theme_minimal() + 
    theme(axis.text.x = element_blank()) +
    facet_wrap(~Migration) +
    ggtitle("(d)") +
    theme(legend.position = "none")
)


ggsave(filename = "plots/fig_s3_trait_shift.jpg", width = 12, height = 7,
       grid.arrange(cred_hist_niche, cred_hist_order, count_mag_mass,
                    cred_hist_migr, right = trait_legend))


hab_props_long <- hab_counts %>%
  pivot_longer(cols = all_of(c(land_covars, "other")),
               names_to = "hab", values_to = "count") %>%
  select(-cond_iter, -cond_type) %>%
  ungroup() %>%
  mutate(count = ifelse(count < 0, 0, count)) %>%
  pivot_wider(id_cols = c("i", "hab", "species"),
              names_from = "type", values_from = "count") %>%
  group_by(species, i) %>%
  mutate(drought_prop = drought / sum(drought),
         nodrought_prop = nodrought / sum(nodrought),
         prop_diff = drought_prop - nodrought_prop) %>%
  ungroup()

hab_props_summary <- hab_props_long %>%
  group_by(species, hab) %>%
  summarize(pd_med = median(prop_diff),
            pd_Q025 = quantile(prop_diff, 0.025),
            pd_Q975 = quantile(prop_diff, 0.975)) %>%
  mutate(cred_type = ifelse(pd_Q975 < 0, "Credible decline",
                            ifelse(pd_Q025 > 0, "Credible increase",
                                   "No credible diff.")))

pd_df <- hab_props_summary %>%
  filter(species %in% spec_w_interactions) %>%
  select(species, hab, pd_med) %>%
  pivot_wider(names_from = "hab", values_from = "pd_med")

pd_mtx <- as.matrix(pd_df[, 2:7])
rownames(pd_mtx) <- pd_df$species

# In which hab is each species most abundant?
spec_hab_maindiff <- data.frame(
  spec = rownames(pd_mtx),
  biggest_inc = colnames(pd_mtx)[apply(pd_mtx, 1, which.max)],
  biggest_dec = colnames(pd_mtx)[apply(pd_mtx, 1, which.min)]
)



hab_props_summary_nodrought <- hab_props_long %>%
  # mutate(habgrp = ifelse(hab %in% c("riparian", "other_natural", "grass_pasture"),
  #                        "Natural", ifelse(hab %in% c("row_field_ag", "perennial_ag"),
  #                                          "Agriculture", "Developed"))) %>%
  group_by(species, hab, i) %>%
  summarize(nodrought_prop = sum(nodrought_prop)) %>%
  summarize(nd_med = median(nodrought_prop),
            nd_Q025 = quantile(nodrought_prop, 0.025),
            nd_Q975 = quantile(nodrought_prop, 0.975)) %>%
  mutate(shifted_dist = species %in% spec_w_interactions) %>%
  group_by(species) %>%
  ungroup() %>%
  filter(species %in% spec_w_interactions)
hab_props_summary_drought <- hab_props_long %>%
  # mutate(habgrp = ifelse(hab %in% c("riparian", "other_natural", "grass_pasture"),
  #                        "Natural", ifelse(hab %in% c("row_field_ag", "perennial_ag"),
  #                                          "Agriculture", "Developed"))) %>%
  group_by(species, hab, i) %>%
  summarize(drought_prop = sum(drought_prop)) %>%
  summarize(nd_med = median(drought_prop),
            nd_Q025 = quantile(drought_prop, 0.025),
            nd_Q975 = quantile(drought_prop, 0.975)) %>%
  mutate(shifted_dist = species %in% spec_w_interactions) %>%
  group_by(species) %>%
  ungroup() %>%
  filter(species %in% spec_w_interactions)

tplist_drought <- list()
tplist_nodrought <- list()
tplist_line <- list()
drought_col <- character(nrow(hab_props_summary_drought) / 6)
nodrought_col <- character(nrow(hab_props_summary_drought) / 6)
line_col <- character(nrow(hab_props_summary_drought) / 6)

for (i in 1:length(unique(hab_props_summary_drought$species))) {
  thisspec <- unique(hab_props_summary_drought$species)[[i]]
  
  temp <- hab_props_summary_drought %>%
    filter(species == unique(hab_props_summary_drought$species)[i]) %>%
    arrange(hab)
  
  drought_col[i] <- ifelse(thisspec %in% spec_w_interactions,
                           "#e31a1c", "#fb9a99")
                           
  tplist_drought[[i]] <- temp$nd_med
  
  temp <- hab_props_summary_nodrought %>%
    filter(species == unique(hab_props_summary_drought$species)[i]) %>%
    arrange(hab)
  
  tplist_nodrought[[i]] <- temp$nd_med
  
  tplist_line[[i]] <- list(tplist_drought[[i]], tplist_nodrought[[i]])
  
  nodrought_col[i] <- ifelse(thisspec %in% spec_w_interactions,
                             "#1f78b4", "#a6cee3")
                             line_col[i] <- ifelse(thisspec %in% spec_w_interactions,
                                                   "#616161", "#969696")
                                                   
}

library(Ternary)

#### Hab proportion shift plots ####
# hab_shift_colors <- c("Decreased" = "#d95f02", "Increased" = "#1b9e77")
hab_shift_colors <- c("Decreased" = "#FD8700", "Increased" = "#444E7E")

riparian_shift_plot <- bind_rows(
  hab_props_summary_drought %>% mutate(condition = "drought"),
  hab_props_summary_nodrought %>% mutate(condition = "non-drought")
) %>% 
  arrange(condition) %>% 
  group_by(species) %>% 
  filter(hab == "riparian") %>% 
  mutate(inc = ifelse(max(nd_med) == nd_med[[1]], "Increased", "Decreased")) %>% 
  ggplot(aes(x = nd_med, y = condition, group = species)) +
  geom_line(aes(col = inc)) +
  geom_point(aes(col = inc, shape = condition), show.legend = FALSE) +
  xlab("") +
  ylab("") +
  scale_shape_manual(values = c(1, 19)) +
  scale_color_manual("Change w. drought",
                     values = hab_shift_colors) +
  theme_minimal() +
  ggtitle("(d) Riparian") +
  theme(plot.title = element_text(size = rel(1)))

legend_hab <- get_legend(riparian_shift_plot)

riparian_shift_plot <- riparian_shift_plot  +
  theme(legend.position = "none")
grass_shift_plot <- bind_rows(
  hab_props_summary_drought %>% mutate(condition = "drought"),
  hab_props_summary_nodrought %>% mutate(condition = "non-drought")
) %>% 
  arrange(condition) %>% 
  group_by(species) %>% 
  filter(hab == "grass_pasture") %>% 
  mutate(inc = ifelse(max(nd_med) == nd_med[[1]], "Increased", "Decreased")) %>% 
  ggplot(aes(x = nd_med, y = condition, group = species)) +
  geom_line(aes(col = inc)) +
  geom_point(aes(col = inc, shape = condition)) +
  xlab("") +
  ylab("") +
  scale_shape_manual(values = c(1, 19)) +
  scale_color_manual("Change w. drought",
                     values = hab_shift_colors) +
  theme_minimal() +
  ggtitle("(e) Grassland") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size = rel(1)))



othernat_shift_plot <- bind_rows(
  hab_props_summary_drought %>% mutate(condition = "drought"),
  hab_props_summary_nodrought %>% mutate(condition = "non-drought")
) %>% 
  arrange(condition) %>% 
  group_by(species) %>% 
  filter(hab == "other_natural") %>% 
  mutate(inc = ifelse(max(nd_med) == nd_med[[1]], "Increased", "Decreased")) %>% 
  ggplot(aes(x = nd_med, y = condition, group = species)) +
  geom_line(aes(col = inc)) +
  geom_point(aes(col = inc, shape = condition)) +
  xlab("") +
  ylab("") +
  scale_shape_manual(values = c(1, 19)) +
  scale_color_manual("Change w. drought",
                     values = hab_shift_colors) +
  theme_minimal() +
  ggtitle("(f) Other natural") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size = rel(1)))


rowfield_shift_plot <- bind_rows(
  hab_props_summary_drought %>% mutate(condition = "drought"),
  hab_props_summary_nodrought %>% mutate(condition = "non-drought")
) %>% 
  arrange(condition) %>% 
  group_by(species) %>% 
  filter(hab == "row_field_ag") %>% 
  mutate(inc = ifelse(max(nd_med) == nd_med[[1]], "Increased", "Decreased")) %>% 
  ggplot(aes(x = nd_med, y = condition, group = species)) +
  geom_line(aes(col = inc)) +
  geom_point(aes(col = inc, shape = condition)) +
  xlab("") +
  ylab("") +
  scale_shape_manual(values = c(1, 19)) +
  scale_color_manual("Change w. drought",
                     values = hab_shift_colors) +
  theme_minimal() +
  ggtitle("(b) Row and field ag") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size = rel(1)))


perennial_shift_plot <- bind_rows(
  hab_props_summary_drought %>% mutate(condition = "drought"),
  hab_props_summary_nodrought %>% mutate(condition = "non-drought")
) %>% 
  arrange(condition) %>% 
  group_by(species) %>% 
  filter(hab == "perennial_ag") %>% 
  mutate(inc = ifelse(max(nd_med) == nd_med[[1]], "Increased", "Decreased")) %>% 
  ggplot(aes(x = nd_med, y = condition, group = species)) +
  geom_line(aes(col = inc)) +
  geom_point(aes(col = inc, shape = condition)) +
  xlab("") +
  ylab("") +
  scale_shape_manual(values = c(1, 19)) +
  scale_color_manual("Change w. drought",
                     values = hab_shift_colors) +
  theme_minimal() +
  ggtitle("(c) Perennial agriculture") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size = rel(1)))


other_shift_plot <- bind_rows(
  hab_props_summary_drought %>% mutate(condition = "drought"),
  hab_props_summary_nodrought %>% mutate(condition = "non-drought")
) %>% 
  arrange(condition) %>% 
  group_by(species) %>% 
  filter(hab == "other") %>% 
  mutate(inc = ifelse(max(nd_med) == nd_med[[1]], "Increased", "Decreased")) %>% 
  ggplot(aes(x = nd_med, y = condition, group = species)) +
  geom_line(aes(col = inc)) +
  geom_point(aes(col = inc, shape = condition)) +
  ylab("") +
  scale_shape_manual(values = c(1, 19)) +
  scale_color_manual("Change w. drought",
                     values = hab_shift_colors) +
  theme_minimal() +
  ggtitle("(g) Developed/other") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size = rel(1)))



ggsave(grid.arrange(rowfield_shift_plot, perennial_shift_plot, riparian_shift_plot,
                    grass_shift_plot, othernat_shift_plot, other_shift_plot,
                    ncol = 1, right = legend_hab), filename = "plots/hab_fracshift_plot.jpg",
       width = 7, height = 10)



library(ggtern)
hab_props_summary_drought <- hab_props_long %>% 
  mutate(habgrp = ifelse(hab %in% c("riparian", "other_natural", "grassland"),
                         "Natural", ifelse(hab %in% c("row_field_ag", "perennial_ag"),
                                           "Agriculture", "Developed"))) %>%
  group_by(species, habgrp, i) %>% 
  summarize(drought_prop = sum(drought_prop)) %>% 
  summarize(nd_med = median(drought_prop),
            nd_Q025 = quantile(drought_prop, 0.025),
            nd_Q975 = quantile(drought_prop, 0.975)) %>% 
  mutate(shifted_dist = species %in% spec_w_interactions) %>% 
  group_by(species) %>% 
  ungroup() %>% 
  filter(species %in% spec_w_interactions)
hab_props_summary_nodrought <- hab_props_long %>% 
  mutate(habgrp = ifelse(hab %in% c("riparian", "other_natural", "grassland"),
                         "Natural", ifelse(hab %in% c("row_field_ag", "perennial_ag"),
                                           "Agriculture", "Developed"))) %>%
  group_by(species, habgrp, i) %>% 
  summarize(nodrought_prop = sum(nodrought_prop)) %>% 
  summarize(nd_med = median(nodrought_prop),
            nd_Q025 = quantile(nodrought_prop, 0.025),
            nd_Q975 = quantile(nodrought_prop, 0.975)) %>% 
  mutate(shifted_dist = species %in% spec_w_interactions) %>% 
  group_by(species) %>% 
  ungroup() %>% 
  filter(species %in% spec_w_interactions)


(tern_plot <- bind_rows(
  hab_props_summary_drought %>% mutate(condition = "Drought"),
  hab_props_summary_nodrought %>% mutate(condition = "Non-drought")
) %>% 
    select(condition, species, nd_med, habgrp) %>% 
    pivot_wider(names_from = habgrp, values_from = nd_med) %>% 
    ggtern(aes(Agriculture, Developed, Natural, group = species)) +
    geom_point(aes(shape = condition )) +
    geom_line() +
    scale_shape_manual("", values = c(1, 19)) +
    xlab("Agriculture") +
    theme_light() + theme_showarrows() +
    theme(legend.position = "bottom",
          tern.axis.arrow = element_line(color = "darkgray", size = 1)) +
    Tlab("", labelarrow = "Developed") +
    Rlab("", labelarrow = "Natural") +
    Llab("", labelarrow = "Agriculture")
)

ly_mtx <- matrix(c(
  1,1,2,2,2,
  1,1,3,3,3,
  1,1,4,4,4,
  1,1,5,5,5,
  1,1,6,6,6,
  1,1,7,7,7
), nrow = 6, byrow = TRUE)

ly_mtx <- matrix(c(
  1,1,1,1,2,2,2,5,5,5,
  1,1,1,1,3,3,3,6,6,6,
  1,1,1,1,4,4,4,7,7,7
), nrow = 3, byrow = TRUE)

ggsave(
  grid.arrange(tern_plot + ggtitle("(a)") + theme(plot.margin = margin(rep(-2, 4), "cm")),
               rowfield_shift_plot, 
               perennial_shift_plot, 
               riparian_shift_plot + xlab("Relative use"), 
               grass_shift_plot, 
               othernat_shift_plot, 
               other_shift_plot + xlab("Relative use"), 
               right = legend_hab, layout_matrix = ly_mtx),
  width = 12, height = 5, filename = "plots/fig3_hab_shift_plot.jpg"
)



#### Summarize LMM results ####

get_hab_slopes <- function(xlist) {
  
  res <- data.frame(
    hab = c(land_covars, "other"),
    b_med = NA,
    b_025 = NA,
    b_975 = NA
  )
  
  for (i in 1:nrow(res)) {
    lc <- res$hab[i]
    
    lc_flag <- c("zero_pct", "onehundred_pct")[as.numeric(land_covars == lc) + 1]
    
    beta_thishab <- xlist$maineff_draws$b_spei + 
      xlist$maineff_draws$`b_spei:row_field_ag` * 
      land_scaling_factors[land_scaling_factors$param == "row_field_ag", 
                           lc_flag[which(land_covars == "row_field_ag")]] +
      xlist$maineff_draws$`b_spei:perennial_ag` * 
      land_scaling_factors[land_scaling_factors$param == "perennial_ag", 
                           lc_flag[which(land_covars == "perennial_ag")]] +
      xlist$maineff_draws$`b_spei:riparian` * 
      land_scaling_factors[land_scaling_factors$param == "riparian", 
                           lc_flag[which(land_covars == "riparian")]] +
      xlist$maineff_draws$`b_spei:grass_pasture` * 
      land_scaling_factors[land_scaling_factors$param == "grass_pasture", 
                           lc_flag[which(land_covars == "grass_pasture")]] +
      xlist$maineff_draws$`b_spei:other_natural` * 
      land_scaling_factors[land_scaling_factors$param == "other_natural", 
                           lc_flag[which(land_covars == "other_natural")]]
    
    res$b_med[i] <- median(beta_thishab)
    res$b_025[i] <- quantile(beta_thishab, probs = c(0.025))
    res$b_975[i] <- quantile(beta_thishab, probs = c(0.975))
  }
  
  res
}


evi_summary <- readRDS("intermediate/EVI_LMM_posterior.RDS") %>% 
  get_hab_slopes()
ndwi_summary <- readRDS("intermediate/NDWI_LMM_posterior.RDS") %>% 
  get_hab_slopes()
tmax_summary <- readRDS("intermediate/tmax_ssn_LMM_posterior.RDS") %>% 
  get_hab_slopes()
precip_summary <- readRDS("intermediate/ppt_LMM_posterior.RDS") %>% 
  get_hab_slopes()

