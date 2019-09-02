nmds_rar_list <-
  com_rar_list %>% 
  mutate(NMDS = map(data, . %>% spread_cdm(LAC_POS, TAXON, COUNT) %>%
                      metaMDS(distance = "bray", trace = FALSE)),
         NMDS_TOP = map(data, . %>% filter(POSITION == "Top") %>%
                          spread_cdm(LAC, TAXON, COUNT) %>%
                          metaMDS(distance = "bray", trace = FALSE)),
         NMDS_BOTTOM = map(data, . %>% filter(POSITION == "Bottom") %>%
                             spread_cdm(LAC, TAXON, COUNT) %>%
                             metaMDS(distance = "bray", trace = FALSE))
  ) %>% 
  mutate(SITE_SCORES = map(NMDS, . %>%
                             scores(display = 'sites') %>%
                             as_tibble(rownames = "LAC_POS")),
         TAXON_SCORES = map(NMDS, . %>%
                              scores(display = 'species') %>%
                              as_tibble(rownames = "TAXON")),
         SITE_SCORES_TOP = map(NMDS_TOP, . %>%
                                 scores(display = 'sites') %>%
                                 as_tibble(rownames = "LAC_POS")),
         TAXON_SCORES_TOP = map(NMDS_TOP, . %>%
                                  scores(display = 'species') %>%
                                  as_tibble(rownames = "TAXON")),
         SITE_SCORES_BOTTOM = map(NMDS_BOTTOM, . %>%
                                    scores(display = 'sites') %>%
                                    as_tibble(rownames = "LAC_POS")),
         TAXON_SCORES_BOTTOM = map(NMDS_BOTTOM, . %>%
                                     scores(display = 'species') %>%
                                     as_tibble(rownames = "TAXON")))

nmds_rar_plots <- nmds_rar_list %>%
  mutate(TOP_TAXON_SCORES_NMDS1 = map(TAXON_SCORES, . %>%
                                        arrange(desc(abs(NMDS1))) %>%
                                        head(10))) %>% 
  mutate(TOP_TAXON_SCORES_NMDS2 = map(TAXON_SCORES, . %>%
                                        arrange(desc(abs(NMDS2))) %>%
                                        head(10))) %>% 
  select(INVENTORY, SITE_SCORES, SITE_SCORES_TOP, SITE_SCORES_BOTTOM,
         TOP_TAXON_SCORES_NMDS1, TOP_TAXON_SCORES_NMDS2) %>% 
  mutate(plots = pmap(select(., INVENTORY, SITE_SCORES, TOP_TAXON_SCORES_NMDS1, TOP_TAXON_SCORES_NMDS2),
                      function(INVENTORY, SITE_SCORES, TOP_TAXON_SCORES_NMDS1, TOP_TAXON_SCORES_NMDS2){
                        TOP_TAXON_SCORES_NMDS1$NMDS1 <- TOP_TAXON_SCORES_NMDS1$NMDS1 * max(SITE_SCORES$NMDS1) / max(abs(TOP_TAXON_SCORES_NMDS1$NMDS1))
                        TOP_TAXON_SCORES_NMDS1$NMDS2 <- TOP_TAXON_SCORES_NMDS2$NMDS2 * max(SITE_SCORES$NMDS2) / max(abs(TOP_TAXON_SCORES_NMDS2$NMDS2))
                        SITE_SCORES %>%
                          left_join(lakes_meta_nr) %>% 
                          left_join(lakes_bv) %>% 
                          filter(!is.na(LAC)) %>% 
                          mutate(LABS_LAC_BOTTOM = ifelse(POSITION == "Bottom", LAC, "")) %>% 
                          ggplot() +
                          geom_line(aes(x = NMDS1, y = NMDS2, group = LAC), color = "lightgrey") +
                          geom_point(aes(x = NMDS1, y = NMDS2, shape = POSITION, color = ALTITUDE), size = 3) +
                          geom_text(aes(x = NMDS1, y = NMDS2, label = LABS_LAC_BOTTOM), size = 3, nudge_y = -0.03) +
                          scale_colour_continuous(type = "viridis") +
                          scale_x_continuous(sec.axis = sec_axis(trans = ~., breaks = TOP_TAXON_SCORES_NMDS1$NMDS1,
                                                                 labels = TOP_TAXON_SCORES_NMDS1$TAXON)) +
                          scale_y_continuous(sec.axis = sec_axis(trans = ~., breaks = TOP_TAXON_SCORES_NMDS2$NMDS2,
                                                                 labels = TOP_TAXON_SCORES_NMDS2$TAXON)) +
                          labs(title = INVENTORY) +
                          theme_bw() +
                          theme(axis.text.y.right = element_text(size = 6, angle = 45),
                                axis.text.x.top = element_text(size = 6, angle = 45, hjust = 0, vjust = 0))
                      }),
         plots_NMDS_TOP = pmap(select(., INVENTORY, SITE_SCORES_TOP),
                               function(INVENTORY, SITE_SCORES_TOP){
                                 SITE_SCORES_TOP %>%
                                   left_join(lakes_meta_nr) %>% 
                                   left_join(lakes_bv) %>% 
                                   filter(!is.na(LAC)) %>% 
                                   ggplot() +
                                   geom_point(aes(x = NMDS1, y = NMDS2, shape = POSITION, color = ALTITUDE), size = 3) +
                                   geom_text(aes(x = NMDS1, y = NMDS2, label = sub("_Top$", "", LAC_POS)), size = 3, nudge_y = -0.03) +
                                   scale_colour_continuous(type = "viridis") +
                                   labs(title = INVENTORY) +
                                   theme_bw()
                               }),
         plots_NMDS_BOTTOM = pmap(select(., INVENTORY, SITE_SCORES_BOTTOM),
                                  function(INVENTORY, SITE_SCORES_BOTTOM){
                                    SITE_SCORES_BOTTOM %>%
                                      left_join(lakes_meta_nr) %>% 
                                      left_join(lakes_bv) %>% 
                                      filter(!is.na(LAC)) %>% 
                                      ggplot() +
                                      geom_point(aes(x = NMDS1, y = NMDS2, shape = POSITION, color = ALTITUDE), size = 3) +
                                      geom_text(aes(x = NMDS1, y = NMDS2, label = sub("_Bottom$", "", LAC_POS)), size = 3, nudge_y = -0.01) +
                                      scale_colour_continuous(type = "viridis") +
                                      labs(title = INVENTORY) +
                                      theme_bw()
                                  }))