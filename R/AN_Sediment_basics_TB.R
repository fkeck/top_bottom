sed_basics_tests <- list()
sed_basics_plots <- list()

theme_boxplot_sed <- theme_grey() +
  theme(legend.position = "none",
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10))


dat <- nmds_rar_list %>%
  filter(INVENTORY == "OTU") %>%
  .$SITE_SCORES %>% .[[1]] %>%
  separate(LAC_POS, into = c("LAC", "POSITION"), sep = "_") %>%
  left_join(lakes_cn, by = c("LAC", "POSITION")) %>% 
  filter(POSITION != "Bottom 2",
         !LAC %in% c("Bourget", "Léman", "Godivelle", "Corbeaux")) %>% 
  semi_join(lakes_meta_nr, by = c("LAC", "POSITION")) %>% 
  mutate(jitter_val = runif(nrow(.), -0.1, 0.1))

sed_basics_plots$DNA <- ggplot(dat) +
  geom_boxplot(aes(as.numeric(as.factor(POSITION)), MASS_SAMPLE_DNA, group = POSITION, fill = POSITION), outlier.shape = NA)+
  geom_point(aes(as.numeric(as.factor(POSITION)) + jitter_val, MASS_SAMPLE_DNA), size = 0.4) +
  geom_line(aes(as.numeric(as.factor(POSITION)) + jitter_val, MASS_SAMPLE_DNA, group = LAC), alpha = 0.2) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Bottom", "Top")) +
  xlab("") + ylab("DNA") +
  theme_boxplot_sed


sed_basics_tests$DNA <- wilcox.test(MASS_SAMPLE_DNA ~ POSITION, paired = TRUE, data = dat, exact = FALSE)


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
  geom_point(aes(as.numeric(as.factor(POSITION)) + jitter_val, C_tot_perc_MS), size = 0.4) +
  geom_line(aes(as.numeric(as.factor(POSITION)) + jitter_val, C_tot_perc_MS, group = LAC), alpha = 0.2) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Bottom", "Top")) +
  xlab("") + ylab("Total carbon") +
  theme_boxplot_sed

sed_basics_tests$COT <- wilcox.test(C_tot_perc_MS ~ POSITION, paired = TRUE, data = dat, exact = FALSE)


dat <- nmds_rar_list %>%
  filter(INVENTORY == "OTU") %>%
  .$SITE_SCORES %>% .[[1]] %>%
  separate(LAC_POS, into = c("LAC", "POSITION"), sep = "_") %>%
  left_join(lakes_cn, by = c("LAC", "POSITION")) %>% 
  filter(POSITION != "Bottom 2",
         LAC != "Mont Coua", !is.na(CAROTENOIDS)) %>% 
  semi_join(lakes_meta_nr, by = c("LAC", "POSITION")) %>% 
  mutate(jitter_val = runif(nrow(.), -0.1, 0.1))

sed_basics_plots$PIG <- ggplot(dat) +
  geom_boxplot(aes(as.numeric(as.factor(POSITION)), log(CAROTENOIDS), group = POSITION, fill = POSITION), outlier.shape = NA)+
  geom_point(aes(as.numeric(as.factor(POSITION)) + jitter_val, log(CAROTENOIDS)), size = 0.4) +
  geom_line(aes(as.numeric(as.factor(POSITION)) + jitter_val, log(CAROTENOIDS), group = LAC), alpha = 0.2) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Bottom", "Top")) +
  xlab("") + ylab("Carotenoids (log-transformed)") +
  theme_boxplot_sed

sed_basics_tests$PIG <- wilcox.test(CAROTENOIDS ~ POSITION, paired = TRUE, data = dat, exact = FALSE)
