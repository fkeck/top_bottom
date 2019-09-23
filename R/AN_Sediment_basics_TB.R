sed_basics_tests <- list()
sed_basics_plots <- list()

theme_boxplot_sed <- theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


dat <- nmds_rar_list %>%
  filter(INVENTORY == "OTU") %>%
  .$SITE_SCORES %>% .[[1]] %>%
  separate(LAC_POS, into = c("LAC", "POSITION"), sep = "_") %>%
  left_join(lakes_cn, by = c("LAC", "POSITION")) %>%
  semi_join(lakes_meta_nr, by = c("LAC", "POSITION")) %>% 
  left_join(lakes_dna) %>% 
    filter(POSITION != "Bottom 2",
           !LAC %in% c("Bourget", "Léman", "Marion", "Chalain", "Vert")) %>% 
  mutate(mug_ADN_BY_g_C_DRY_SED = (ng_ADN_BY_g_DRY_SED/1000) * (1/(C_tot_perc_MS/100)),
         mug_ADN_BY_g_Corg_DRY_SED = (ng_ADN_BY_g_DRY_SED/1000) * (1/(C_org_perc_MS/100))) %>%
  mutate(jitter_val = runif(nrow(.), -0.1, 0.1))

sed_basics_plots$DNA <- ggplot(dat) +
  geom_boxplot(aes(as.numeric(as.factor(POSITION)), (ng_ADN_BY_g_DRY_SED/1000) , group = POSITION, fill = POSITION), outlier.shape = NA)+
  geom_line(aes(as.numeric(as.factor(POSITION)), (ng_ADN_BY_g_DRY_SED/1000), group = LAC), alpha = 0.2) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Past", "Recent")) +
  xlab("") + ylab(expression(paste("DNA (", µg~g^-1, " dry sediment)"))) +
  theme_boxplot_sed

sed_basics_plots$DNA_corrected_Ctot <- ggplot(dat) +
  geom_boxplot(aes(as.numeric(as.factor(POSITION)), mug_ADN_BY_g_C_DRY_SED , group = POSITION, fill = POSITION), outlier.shape = NA)+
  geom_line(aes(as.numeric(as.factor(POSITION)), mug_ADN_BY_g_C_DRY_SED, group = LAC), alpha = 0.2) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Past", "Recent")) +
  xlab("") + ylab(expression(paste("DNA (", µg~g^-1, "TC dry sediment)"))) +
  theme_boxplot_sed

sed_basics_plots$DNA_corrected_Corg <- ggplot(dat) +
  geom_boxplot(aes(as.numeric(as.factor(POSITION)), mug_ADN_BY_g_Corg_DRY_SED , group = POSITION, fill = POSITION), outlier.shape = NA)+
  geom_line(aes(as.numeric(as.factor(POSITION)), mug_ADN_BY_g_Corg_DRY_SED, group = LAC), alpha = 0.2) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Past", "Recent")) +
  xlab("") + ylab(expression(paste("DNA (", µg~g^-1, "TOC dry sediment)"))) +
  theme_boxplot_sed


sed_basics_tests$DNA <- wilcox.test((ng_ADN_BY_g_DRY_SED/1000) ~ POSITION, paired = TRUE, data = dat, exact = FALSE)
sed_basics_tests$DNA_corrected_Ctot <- wilcox.test(mug_ADN_BY_g_C_DRY_SED ~ POSITION, paired = TRUE, data = dat, exact = FALSE)
sed_basics_tests$DNA_corrected_Corg <- wilcox.test(mug_ADN_BY_g_Corg_DRY_SED ~ POSITION, paired = TRUE, data = dat, exact = FALSE)



dat <- nmds_rar_list %>%
  filter(INVENTORY == "OTU") %>%
  .$SITE_SCORES %>% .[[1]] %>%
  separate(LAC_POS, into = c("LAC", "POSITION"), sep = "_") %>%
  left_join(lakes_cn, by = c("LAC", "POSITION")) %>% 
  filter(POSITION != "Bottom 2",
         !LAC %in% c("Bourget", "Léman", "Marion", "Vert")) %>% 
  semi_join(lakes_meta_nr, by = c("LAC", "POSITION")) %>% 
  mutate(jitter_val = runif(nrow(.), -0.1, 0.1))

sed_basics_plots$COT <- ggplot(dat) +
  geom_boxplot(aes(as.numeric(as.factor(POSITION)), C_tot_perc_MS, group = POSITION, fill = POSITION), outlier.shape = NA)+
  geom_line(aes(as.numeric(as.factor(POSITION)), C_tot_perc_MS, group = LAC), alpha = 0.2) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Past", "Recent")) +
  xlab("") + ylab("Total carbon (% of dry mass)") +
  theme_boxplot_sed

sed_basics_plots$Corg <- ggplot(dat) +
  geom_boxplot(aes(as.numeric(as.factor(POSITION)), C_org_perc_MS, group = POSITION, fill = POSITION), outlier.shape = NA)+
  geom_line(aes(as.numeric(as.factor(POSITION)), C_org_perc_MS, group = LAC), alpha = 0.2) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Past", "Recent")) +
  xlab("") + ylab("TOC (% of dry mass)") +
  theme_boxplot_sed

sed_basics_tests$COT <- wilcox.test(C_tot_perc_MS ~ POSITION, paired = TRUE, data = dat, exact = FALSE)
sed_basics_tests$Corg <- wilcox.test(C_org_perc_MS ~ POSITION, paired = TRUE, data = dat, exact = FALSE)


dat <- nmds_rar_list %>%
  filter(INVENTORY == "OTU") %>%
  .$SITE_SCORES %>% .[[1]] %>%
  separate(LAC_POS, into = c("LAC", "POSITION"), sep = "_") %>%
  left_join(lakes_cn, by = c("LAC", "POSITION")) %>% 
  filter(POSITION != "Bottom 2",
         !LAC %in% c("Bourget", "Léman", "Marion", "Mont Coua", "Vert")) %>% 
  semi_join(lakes_meta_nr, by = c("LAC", "POSITION")) %>% 
  mutate(CAROTENOIDS_BY_g_C_DRY_SED = CAROTENOIDS * (1/(C_tot_perc_MS/100)),
         CAROTENOIDS_BY_g_Corg_DRY_SED = CAROTENOIDS * (1/(C_org_perc_MS/100))) %>%
  mutate(jitter_val = runif(nrow(.), -0.1, 0.1))

sed_basics_plots$PIG <- ggplot(dat) +
  geom_boxplot(aes(as.numeric(as.factor(POSITION)), log(CAROTENOIDS), group = POSITION, fill = POSITION), outlier.shape = NA)+
  geom_line(aes(as.numeric(as.factor(POSITION)), log(CAROTENOIDS), group = LAC), alpha = 0.2) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Past", "Recent")) +
  xlab("") + ylab(expression(paste("Total carotenoids (", ng~g^-1, " dry sediment, log-transf.)"))) +
  theme_boxplot_sed

sed_basics_plots$PIG_corrected_Ctot <- ggplot(dat) +
  geom_boxplot(aes(as.numeric(as.factor(POSITION)), log(CAROTENOIDS_BY_g_C_DRY_SED), group = POSITION, fill = POSITION), outlier.shape = NA)+
  geom_line(aes(as.numeric(as.factor(POSITION)), log(CAROTENOIDS_BY_g_C_DRY_SED), group = LAC), alpha = 0.2) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Past", "Recent")) +
  xlab("") + ylab(expression(paste("Total carotenoids (", ng~g^-1, "TC dry sediment, log-transf.)"))) +
  theme_boxplot_sed

sed_basics_plots$PIG_corrected_Corg <- ggplot(dat) +
  geom_boxplot(aes(as.numeric(as.factor(POSITION)), log(CAROTENOIDS_BY_g_Corg_DRY_SED), group = POSITION, fill = POSITION), outlier.shape = NA)+
  geom_line(aes(as.numeric(as.factor(POSITION)), log(CAROTENOIDS_BY_g_Corg_DRY_SED), group = LAC), alpha = 0.2) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Past", "Recent")) +
  xlab("") + ylab(expression(paste("Total carotenoids (", ng~g^-1, "TOC dry sediment, log-transf.)"))) +
  theme_boxplot_sed

sed_basics_tests$PIG <- wilcox.test(CAROTENOIDS ~ POSITION, paired = TRUE, data = dat, exact = FALSE)
sed_basics_tests$PIG_corrected_Ctot <- wilcox.test(CAROTENOIDS_BY_g_C_DRY_SED ~ POSITION, paired = TRUE, data = dat, exact = FALSE)
sed_basics_tests$PIG_corrected_Corg <- wilcox.test(CAROTENOIDS_BY_g_Corg_DRY_SED ~ POSITION, paired = TRUE, data = dat, exact = FALSE)




