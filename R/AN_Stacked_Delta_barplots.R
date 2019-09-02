#### STACKED and DELTA BARPLOTS OF TAXO DIVERSITY ####

nperc <- function(x) scales::number(x, suffix = "%")

plot_stack_taxo <- function(x, title){
  x %>% 
    mutate(FREQ = ifelse(POSITION == "Bottom", FREQ * -1, FREQ)) %>% 
    ggplot() +
    geom_col(aes(x = forcats::fct_rev(LAC), y = FREQ, fill = TAXON), color = "grey", size = 0.3) +
    ggpol::facet_share(~POSITION, dir = "h", scales = "free", reverse_num = FALSE) +
    coord_flip() +
    theme_minimal() +
    theme(legend.position="bottom") +
    xlab("") + labs(title = title) +
    scale_y_continuous(breaks = c(-100, -75, -50, -25, 0, 25, 50, 75, 100),
                       labels = c("100%", "75%", "50%", "25%", "0%", "25%", "50%", "75%", "100%"))
}

stack_taxo_plots <- com_rar_list %>% 
  filter(!str_detect(com_rar_list$INVENTORY, "^OTU")) %>% 
  mutate(plots = map2(data, INVENTORY, plot_stack_taxo)) %>% 
  mutate(tb_stack_plots = map2(data, INVENTORY, function(data, INVENTORY){
    group_by(data, POSITION, TAXON) %>% 
      summarise(FREQ = mean(FREQ)) %>% 
      ggplot() +
      geom_col(aes(POSITION, FREQ, fill = TAXON), color = "grey", size = 0.3) +
      xlab("Position") + ylab("Mean Count") + labs(title = INVENTORY) +
      scale_y_continuous(labels = nperc)
  })) %>% 
  mutate(tb_bar_plots = map2(data, INVENTORY, function(data, INVENTORY){
    group_by(data, LAC, TAXON) %>% 
      arrange(LAC, TAXON, desc(POSITION)) %>% 
      mutate_at("FREQ", list(~ . - lag(., default = 0))) %>% 
      summarise_at("FREQ", tail, 1) %>% 
      mutate(FREQ = FREQ * -1) %>% 
      group_by(TAXON) %>% 
      summarise(MEAN_DELTA = mean(FREQ), SD_DELTA = sd(FREQ), N = n()) %>% 
      ggplot() +
      geom_hline(yintercept = 0) +
      geom_col(aes(fct_reorder(TAXON, MEAN_DELTA), MEAN_DELTA)) +
      geom_linerange(aes(x = fct_reorder(TAXON, MEAN_DELTA),
                         ymin = MEAN_DELTA + 1.96 * SD_DELTA / sqrt(N),
                         ymax = MEAN_DELTA - 1.96 * SD_DELTA / sqrt(N))) +
      xlab("Lacs") + ylab("Count Delta Top-Bottom (mean)") + labs(title = INVENTORY) +
      theme(axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1)) +
      scale_y_continuous(labels = nperc)
  }))

tb_delta_taxo_altitude <- com_rar_list %>% 
  filter(INVENTORY %in% c("Taxonomie Rank 3", "Trophic Mode", "Trophic Type")) %>% 
  mutate(tb_bar_plots = map2(data, INVENTORY, function(data, INVENTORY){
    group_by(data, LAC, TAXON) %>% 
      arrange(LAC, TAXON, desc(POSITION)) %>% 
      mutate_at("FREQ", list(~ . - lag(., default = 0))) %>% 
      summarise_at("FREQ", tail, 1) %>% 
      mutate(FREQ = FREQ * -1) %>% 
      left_join(lakes_bv) %>% 
      ggplot() +
      geom_point(aes(ALTITUDE, FREQ)) +
      labs(title = INVENTORY) +
      facet_wrap(~TAXON, scales = "free", strip.position = "bottom") +
      ylab("Delta Top-Bottom") + xlab("")
  }))
