

rar_richness <- com_merged_raw_list %>% 
  mutate(RAR_RICHNESS = map(data, function(data){
    dat <- spread_cdm(data, LAC_POS, TAXON, COUNT)
    raremax <- dat %>% rowSums() %>% min()
    res <- rarefy(dat, sample = raremax) %>% 
      enframe("LAC_POS", "RAR_RICHNESS") %>% 
      left_join(lakes_meta_nr)
    return(res)
  })) %>% 
  mutate(MEAN_RAR_RICHNESS = map(RAR_RICHNESS, function(RAR_RICHNESS){
    RAR_RICHNESS %>% 
      group_by(POSITION) %>% 
      summarise(MEAN_RAR_RICH = mean(RAR_RICHNESS))
  })) %>% 
  mutate(PLOT_RAR_RICHNESS = map(RAR_RICHNESS, function(RAR_RICHNESS){
    RAR_RICHNESS %>% 
      ggplot() +
      geom_boxplot(aes(POSITION, RAR_RICHNESS))
  }))%>% 
  mutate(TEST_RAR_RICHNESS = map(RAR_RICHNESS, function(RAR_RICHNESS){
    val_top <- RAR_RICHNESS$RAR_RICHNESS[RAR_RICHNESS$POSITION == "Top"]
    val_bot <- RAR_RICHNESS$RAR_RICHNESS[RAR_RICHNESS$POSITION == "Bottom"]
    wilcox.test(val_top, val_bot, paired = TRUE)
  }))





es_richness <- com_merged_raw_list %>% 
  mutate(ES_RICHNESS = map(data, function(data){
    spread_cdm(data, LAC_POS, TAXON, COUNT) %>%
      estimateR() %>%
      t() %>%
      set_colnames(c("S_OBS", "S_CHAO1", "SE_CHAO1", "S_ACE", "SE_ACE")) %>% 
      as_tibble(rownames = "LAC_POS") %>% 
      left_join(lakes_meta_nr)
  })) %>% 
  mutate(MEAN_S_CHAO1 = map(ES_RICHNESS, function(ES_RICHNESS){
    ES_RICHNESS %>% 
      group_by(POSITION) %>% 
      summarise(MEAN_RAR_RICH = mean(S_CHAO1))
  })) %>% 
  mutate(PLOT_S_CHAO1 = map(ES_RICHNESS, function(ES_RICHNESS){
    ES_RICHNESS %>% 
      ggplot() +
      geom_boxplot(aes(POSITION, S_CHAO1))
  }))%>% 
  mutate(TEST_S_CHAO1 = map(ES_RICHNESS, function(ES_RICHNESS){
    val_top <- ES_RICHNESS$S_CHAO1[ES_RICHNESS$POSITION == "Top"]
    val_bot <- ES_RICHNESS$S_CHAO1[ES_RICHNESS$POSITION == "Bottom"]
    wilcox.test(val_top, val_bot, paired = TRUE)
  }))
