general_stats <- list()

general_stats$N_reads <- otu_merge_tbl$OTU_COUNT %>% sum()

general_stats$N_OTU <- OTU_all %>% .$TAXON %>% unique() %>% length()
general_stats$N_OTU_BOT <- OTU_Bottom %>% .$TAXON %>% unique() %>% length()
general_stats$N_OTU_TOP <- OTU_Top %>% .$TAXON %>% unique() %>% length()
general_stats$N_OTU_TOP_INTER_BOT <- intersect(OTU_Top %>% .$TAXON %>% unique(), OTU_Bottom %>% .$TAXON %>% unique()) %>% length()
general_stats$N_OTU_BOT_ONLY <- general_stats$N_OTU - general_stats$N_OTU_TOP
general_stats$N_OTU_TOP_ONLY <- general_stats$N_OTU - general_stats$N_OTU_BOT

general_stats$N_OTU_UNCLASSIF_TAX <- otu_merge_tbl$OTU_ID %>% unique() %>% enframe(value = "OTU_ID") %>% left_join(otu_meta) %>% group_by(ADL_RANK_2) %>% summarise(n = n()) %>% filter(is.na(ADL_RANK_2)) %>% .$n
general_stats$N_OTU_UNCLASSIF_TROPH <- otu_merge_tbl$OTU_ID %>% unique() %>% enframe(value = "OTU_ID") %>% left_join(otu_meta) %>% group_by(TROPHIC_TYPE) %>% summarise(n = n()) %>% filter(is.na(TROPHIC_TYPE)) %>% .$n