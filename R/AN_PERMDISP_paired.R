# Homogeneity of multivariate dispersions (change in dispersion)

beta_div_TB <- com_rar_list %>% 
  #filter(INVENTORY == "OTU") %>% 
  mutate(BRAY_DIST = map(data, ~ .x %>% 
                           spread_cdm(LAC_POS, TAXON, COUNT) %>% 
                           vegdist(method = "jaccard"))) %>%
  mutate(BETA_DISP = map(BRAY_DIST, ~ .x %>% 
                           betadisper(group = ifelse(str_detect(attr(., "Labels"), "_Top"),
                                                     "Top", "Bottom")))
  ) %>% 
  mutate(BETA_DISP_TEST = map(BETA_DISP, function(x) {
    wilcox.test(x = x$distances[x$group == "Bottom"],
                y = x$distances[x$group == "Top"],
                paired = TRUE)
  })) %>% 
  mutate(TB_intra_plots = map2(BETA_DISP, INVENTORY, function(BETA_DISP, INVENTORY){
    BETA_DISP <- tibble(DIST = BETA_DISP$distances, POSITION = BETA_DISP$group)
    ggplot(BETA_DISP, aes(DIST, fill = POSITION)) +
      geom_histogram(alpha = .5, position = "identity", bins = 20) +
      xlab("Distance to group geometric median") + ylab("Count") +
      labs(fill = "Position") +
      theme_bw() +
      theme(legend.position = "none")
    }))
