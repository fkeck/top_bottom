

TB_pw_BC_plots <- com_rar_list %>%
  mutate(plots = map2(BRAY_DIST, INVENTORY, function(BRAY_DIST, INVENTORY) {
    BRAY_DIST %>%
      filter(LAC_POS_1 != LAC_POS_2, LAC_1 == LAC_2) %>%
      ggplot() +
      geom_col(aes(forcats::fct_reorder(LAC_1, BRAY_DIST), BRAY_DIST)) +
      xlab("Lacs") + ylab("Top-Bottom BC dissimilarity") +
      labs(title = INVENTORY) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  })) %>%
  mutate(hist_plots = map2(BRAY_DIST, INVENTORY, function(BRAY_DIST, INVENTORY) {
    BRAY_DIST %>%
      filter(LAC_POS_1 != LAC_POS_2, LAC_1 == LAC_2) %>%
      ggplot() +
      geom_histogram(aes(BRAY_DIST), bins = 20) +
      xlab("Bray-Curtis dissimilarity") + ylab("Count") +
      theme_bw()
      #labs(title = INVENTORY)
  }))
# 
# 
# TB_pw_BC_covar_plots <- com_rar_list %>% 
#   unnest(BRAY_DIST) %>% 
#   filter(LAC_POS_1 != LAC_POS_2, LAC_1 == LAC_2) %>% 
#   select(INVENTORY, LAC_1, BRAY_DIST) %>% 
#   left_join(select(lakes_bv, LAC, ALTITUDE),
#             by = c("LAC_1" = "LAC")) %>%
#   left_join(select(filter(lakes_cn, POSITION == "Top"), LAC, C_org_perc_MS, C_tot_perc_MS),
#             by = c("LAC_1" = "LAC")) %>% 
#   dplyr::rename(LAC = LAC_1) %>% 
#   gather(key = "VAR", "VAR_VALUE", -INVENTORY, -LAC, -BRAY_DIST) %>% 
#   group_by(INVENTORY) %>% 
#   nest() %>% 
#   mutate(bray_tb_covar_plots = map2(data, INVENTORY,
#                                     function(data, INVENTORY) {
#                                       mod <- data %>%
#                                         group_by(VAR) %>%
#                                         nest() %>%
#                                         mutate(mod = map(data, ~lm(.x$BRAY_DIST ~ .x$VAR_VALUE))) %>% 
#                                         mutate(R2 = map_dbl(mod, ~summary(.x)$r.squared),
#                                                pv = map_dbl(mod, ~summary(.x)$coefficients[2, 4]),
#                                                stat_txt = paste0("R2 = ", round(R2, 2), "\npval = ", scales::pvalue(pv)))
#                                       
#                                       ggplot(data = data) +
#                                         geom_smooth(aes(VAR_VALUE, BRAY_DIST), se = FALSE, color = "darkgrey", size = 0.8) +
#                                         geom_smooth(aes(VAR_VALUE, BRAY_DIST), se = FALSE, color = "black", size = 0.1, method = "lm") +
#                                         #geom_point(aes(VAR_VALUE, BRAY_DIST)) +
#                                         geom_text(aes(VAR_VALUE, BRAY_DIST, label = LAC), size = 2.3) +
#                                         geom_text(aes(-Inf, -Inf, label = stat_txt), data = mod, hjust = -0.1, vjust = -0.5, size = 3) +
#                                         labs(title = INVENTORY) +
#                                         facet_wrap(~VAR, scales = "free", strip.position = "bottom") +
#                                         ylab("Bray-Curtis dissimilarity Top-Bottom") + xlab("")
#                                     }))
# 
# TB_pw_BC_covar_plots <-
#   TB_pw_BC_covar_plots %>% 
#   mutate(lm_mod = map2(data, INVENTORY,
#                        function(data, INVENTORY) {
#                          dat <- spread(data, key = VAR, value = VAR_VALUE) %>% 
#                            filter(complete.cases(.))
#                          full_mod <- lm(BRAY_DIST ~ ALTITUDE + Anthropised_surface + LOG_VOLUME_PYR + C_tot_perc_MS, data = dat, na.action = "na.fail")
#                          all_mod <- MuMIn::dredge(full_mod, rank ="AICc", extra = c("R^2", "adjR^2"))
#                          best_mod <- MuMIn::get.models(all_mod, 1)[[1]]
#                          
#                          rib <- relaimpo::booteval.relimp(relaimpo::boot.relimp(full_mod, b = 100))
#                          plot <- tibble(var = fct_reorder(as.factor(rib@namen[-1]), rib@lmg, .desc = TRUE),
#                                         lmg = rib@lmg,
#                                         lmg_lower = as.vector(rib@lmg.lower),
#                                         lmg_upper = as.vector(rib@lmg.upper),
#                                         AIC_select = ifelse(rib@namen[-1] %in% attr(best_mod$terms, "term.labels"), "In", "Out")) %>% 
#                            ggplot() +
#                            geom_col(aes(var, lmg, fill = AIC_select)) +
#                            geom_linerange(aes(var, ymin = lmg_lower, ymax = lmg_upper)) +
#                            scale_fill_manual(values = c("In" = "#F8766D", "Out" = "#00BFC4")) +
#                            labs(title = INVENTORY)
#                          tibble(lm_mod_all = list(all_mod), lm_mod_plot = list(plot))
#                        })) %>%
#   mutate(tree_mod = map2(data, INVENTORY,
#                          function(data, INVENTORY) {
#                            dat <- spread(data, key = VAR, value = VAR_VALUE)
#                            tree_mod <- rpart(BRAY_DIST ~ ALTITUDE + Anthropised_surface + LOG_VOLUME_PYR + C_tot_perc_MS,
#                                              data = dat, model = TRUE)
#                            tree_mod <- prune(tree_mod, cp = tree_mod$cptable[which.min(tree_mod$cptable[, 4]), 1])
#                            tree_mod$title <- INVENTORY
#                            varimp_plot <- tree_mod$variable.importance %>% 
#                              enframe() %>% 
#                              mutate(CV_select = ifelse(name %in% unique(tree_mod$frame$var[!tree_mod$frame$var == "<leaf>"]), "In", "Out")) %>% 
#                              ggplot() +
#                              geom_col(aes(fct_reorder(as.factor(name), value, .desc = TRUE), value, fill = CV_select)) +
#                              scale_fill_manual(values = c("In" = "#F8766D", "Out" = "#00BFC4")) +
#                              labs(title = INVENTORY)
#                            tibble(tree_mod = list(tree_mod), tree_mod_varimp_plot = list(varimp_plot))
#                          })) %>%
#   unnest(lm_mod, tree_mod)



# rich <- com_merged_raw_list %>% 
#   filter(INVENTORY == "OTU") %>% 
#   .$data %>% .[[1]] %>% 
#   spread_cdm(LAC_POS, TAXON, COUNT) %>% 
#   rarefy(sample = 10000) %>% 
#   enframe("LAC_POS", "RAR_RICHNESS") %>% 
#   left_join(lakes_meta_nr)
# 
# fun_min <- function(x) {x[2] - x[1]}
# 
# d_rich <- rich %>% 
#   group_by(LAC) %>% 
#   summarise(delta_rich = fun_min(RAR_RICHNESS))
# 
# dhill <- com_rar_list$BRAY_DIST[[1]] %>% 
#   filter(LAC_POS_1 != LAC_POS_2, LAC_1 == LAC_2)
# 
# d_rich %>% 
#   left_join(dhill, by = c("LAC" = "LAC_1")) %>% 
#   ggplot() +
#   geom_point(aes(delta_rich, BRAY_DIST)) +
#   geom_text(aes(delta_rich, BRAY_DIST, label = LAC), nudge_y = -0.01) +
#   xlab("Delta Richness")
# 
# 
# d_rich %>% 
#   left_join(lakes_bv) %>% 
#   ggplot() +
#   geom_point(aes(ALTITUDE, delta_rich))
