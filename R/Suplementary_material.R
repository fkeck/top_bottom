SUPP_MAT <- list()


# TAXONOMY RANK 2

FIG_nmds_main <- nmds_rar_list %>% 
  filter(INVENTORY == "Taxonomie Rank 3") %>% 
  .$SITE_SCORES %>% 
  .[[1]] %>% 
  left_join(lakes_meta_nr, by = "LAC_POS") %>% 
  left_join(lakes_bv, by = "LAC") %>% 
  filter(!is.na(LAC)) %>% 
  mutate(LABS_LAC_BOTTOM = ifelse(POSITION == "Bottom", LAC, "")) %>% 
  ggplot() +
  geom_line(aes(x = NMDS1, y = NMDS2, group = LAC), color = "lightgrey") +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, fill = POSITION), geom = "polygon", level = 0.95, alpha = 0.2) +
  geom_point(aes(x = NMDS1, y = NMDS2, shape = POSITION, color = POSITION), size = 2) +
  geom_text(aes(x = NMDS1, y = NMDS2, label = LABS_LAC_BOTTOM), size = 2, nudge_y = -0.03) +
  theme_bw() +
  theme(axis.text.y.right = element_text(size = 6, angle = 45),
        axis.text.x.top = element_text(size = 6, angle = 45, hjust = 0, vjust = 0),
        legend.position = "none")

FIG_nmds_bray <- TB_pw_BC_plots %>% 
  filter(INVENTORY == "Taxonomie Rank 3") %>% .$hist_plots %>% .[[1]] + xlab("Bray-Curtis dissimilarity")

FIG_nmds_intra <- beta_div_TB %>% 
  filter(INVENTORY == "Taxonomie Rank 3") %>% .$TB_intra_plots %>% .[[1]]

dat <- com_rar_list %>% 
  filter(INVENTORY == "Taxonomie Rank 3") %>% .$BRAY_DIST %>%
  .[[1]] %>% 
  filter(LAC_POS_1 != LAC_POS_2, LAC_1 == LAC_2) %>% 
  left_join(select(lakes_bv, LAC, ALTITUDE, LOG_DEPTH_MAX, LOG_SURFACE, LOG_VOLUME_PYR,
                   Artificial_surface, Agricultural_surface, Natural_surface,
                   Anthropised_surface, Natural_Artificial_ratio, Natural_Anthropised_ratio),
            by = c("LAC_1" = "LAC"))

tr_mod_alt <- rpart(BRAY_DIST ~ ALTITUDE, data = dat, model = TRUE)
alt_split <- tr_mod_alt$splits[1, "index"]

bray_split <- dat %>%
  mutate(ALT_SPLIT = ifelse(ALTITUDE > alt_split, "High altitude", "Low altitude")) %>% 
  group_by(ALT_SPLIT) %>% 
  summarise(mean_BRAY = mean(BRAY_DIST),
            sd_BRAY = sd(BRAY_DIST),
            n_BRAY = n()) %>% 
  mutate(CI_BRAY_SUP = mean_BRAY + 1.96 * sd_BRAY / sqrt(n_BRAY),
         CI_BRAY_INF = mean_BRAY - 1.96 * sd_BRAY / sqrt(n_BRAY))


FIG_alt_tree_regression <- dat %>%
  ggplot() +
  geom_rect(aes(xmin = 0, xmax = alt_split,
                ymin = bray_split$CI_BRAY_INF[2], ymax = bray_split$CI_BRAY_SUP[2]), fill = "grey") +
  geom_rect(aes(xmin = alt_split, xmax = max(dat$ALTITUDE),
                ymin = bray_split$CI_BRAY_INF[1], ymax = bray_split$CI_BRAY_SUP[1]), fill = "grey") +
  geom_point(aes(ALTITUDE, BRAY_DIST)) +
  geom_segment(aes(x = 0, y = bray_split[2, "mean_BRAY", drop = TRUE],
                   xend = alt_split, yend = bray_split[2, "mean_BRAY", drop = TRUE]), size = 1) +
  geom_segment(aes(x = alt_split, y = bray_split[1, "mean_BRAY", drop = TRUE],
                   xend = max(dat$ALTITUDE), yend = bray_split[1, "mean_BRAY", drop = TRUE]), size = 1) +
  xlab("Altitude (m)") + ylab("Bray-Curtis dissimilarity")

SUPP_MAT$TAXO <- plot_grid(FIG_nmds_main, FIG_nmds_bray, FIG_nmds_intra, FIG_alt_tree_regression, labels = "AUTO")




##### JACCARD (Switch to jaccard first) ####

FIG_nmds_main <- nmds_rar_list %>% 
  filter(INVENTORY == "OTU") %>% 
  .$SITE_SCORES %>% 
  .[[1]] %>% 
  left_join(lakes_meta_nr, by = "LAC_POS") %>% 
  left_join(lakes_bv, by = "LAC") %>% 
  filter(!is.na(LAC)) %>% 
  mutate(LABS_LAC_BOTTOM = ifelse(POSITION == "Bottom", LAC, "")) %>% 
  ggplot() +
  geom_line(aes(x = NMDS1, y = NMDS2, group = LAC), color = "lightgrey") +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, fill = POSITION), geom = "polygon", level = 0.95, alpha = 0.2) +
  geom_point(aes(x = NMDS1, y = NMDS2, shape = POSITION, color = POSITION), size = 2) +
  geom_text(aes(x = NMDS1, y = NMDS2, label = LABS_LAC_BOTTOM), size = 2, nudge_y = -0.03) +
  theme_bw() +
  theme(axis.text.y.right = element_text(size = 6, angle = 45),
        axis.text.x.top = element_text(size = 6, angle = 45, hjust = 0, vjust = 0),
        legend.position = "none")

FIG_nmds_bray <- TB_pw_BC_plots %>% 
  filter(INVENTORY == "OTU") %>% .$hist_plots %>% .[[1]] + xlab("Jaccard dissimilarity")

FIG_nmds_intra <- beta_div_TB %>% 
  filter(INVENTORY == "OTU") %>% .$TB_intra_plots %>% .[[1]]


dat <- com_rar_list %>% 
  filter(INVENTORY == "OTU") %>% .$BRAY_DIST %>%
  .[[1]] %>% 
  filter(LAC_POS_1 != LAC_POS_2, LAC_1 == LAC_2) %>% 
  left_join(select(lakes_bv, LAC, ALTITUDE, LOG_DEPTH_MAX, LOG_SURFACE, LOG_VOLUME_PYR,
                   Artificial_surface, Agricultural_surface, Natural_surface,
                   Anthropised_surface, Natural_Artificial_ratio, Natural_Anthropised_ratio),
            by = c("LAC_1" = "LAC"))

tr_mod_alt <- rpart(BRAY_DIST ~ ALTITUDE, data = dat, model = TRUE)
alt_split <- tr_mod_alt$splits[1, "index"]

bray_split <- dat %>%
  mutate(ALT_SPLIT = ifelse(ALTITUDE > alt_split, "High altitude", "Low altitude")) %>% 
  group_by(ALT_SPLIT) %>% 
  summarise(mean_BRAY = mean(BRAY_DIST),
            sd_BRAY = sd(BRAY_DIST),
            n_BRAY = n()) %>% 
  mutate(CI_BRAY_SUP = mean_BRAY + 1.96 * sd_BRAY / sqrt(n_BRAY),
         CI_BRAY_INF = mean_BRAY - 1.96 * sd_BRAY / sqrt(n_BRAY))


FIG_alt_tree_regression <- dat %>%
  ggplot() +
  geom_rect(aes(xmin = 0, xmax = alt_split,
                ymin = bray_split$CI_BRAY_INF[2], ymax = bray_split$CI_BRAY_SUP[2]), fill = "grey") +
  geom_rect(aes(xmin = alt_split, xmax = max(dat$ALTITUDE),
                ymin = bray_split$CI_BRAY_INF[1], ymax = bray_split$CI_BRAY_SUP[1]), fill = "grey") +
  geom_point(aes(ALTITUDE, BRAY_DIST)) +
  geom_segment(aes(x = 0, y = bray_split[2, "mean_BRAY", drop = TRUE],
                   xend = alt_split, yend = bray_split[2, "mean_BRAY", drop = TRUE]), size = 1) +
  geom_segment(aes(x = alt_split, y = bray_split[1, "mean_BRAY", drop = TRUE],
                   xend = max(dat$ALTITUDE), yend = bray_split[1, "mean_BRAY", drop = TRUE]), size = 1) +
  xlab("Altitude (m)") + ylab("Jaccard dissimilarity")


SUPP_MAT$JACC <- plot_grid(FIG_nmds_main, FIG_nmds_bray, FIG_nmds_intra, FIG_alt_tree_regression, labels = "AUTO")

save(SUPP_MAT, file = "Supp_mat.Rdata")
