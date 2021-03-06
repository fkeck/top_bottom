
#### ANOVA 1 TB ####
deseq_TB_AOV1 <- com_merged_raw_list %>% 
  filter(!str_detect(INVENTORY, "OTU")) %>% 
  mutate(results = map(data, function(data){
    tab <- data %>% 
      mutate(TAXON = fct_explicit_na(TAXON, "Unclassified")) %>% 
      spread_cdm(LAC_POS, TAXON, COUNT) %>% 
      t()
    cdat <- colnames(tab) %>%
      enframe(value = "LAC_POS") %>%
      left_join(lakes_meta_nr) %>%
      left_join(lakes_bv) %>% 
      mutate(POSITION = factor(POSITION))
    deseq_res <- DESeqDataSetFromMatrix(tab, colData = cdat,
                                        design = ~ POSITION) %>%
      DESeq(fitType = "mean") %>%
      results(independentFiltering = FALSE, contrast = c("POSITION", "Top", "Bottom")) %>%
      as_tibble(rownames = "Taxa")
    deseq_res
  }))


deseq_TB_AOV1 <- deseq_TB_AOV1 %>% 
  mutate(barplots = map2(INVENTORY, results, function(INVENTORY, results){
  results %>%
      mutate(Wald_test = ifelse(padj < 0.05, "Significant", "Non-significant")) %>% 
      ggplot() +
      geom_col(aes(fct_reorder(Taxa, log2FoldChange), log2FoldChange, fill = Wald_test)) +
      geom_linerange(aes(x = fct_reorder(Taxa, log2FoldChange),
                         ymin = log2FoldChange + lfcSE,
                         ymax = log2FoldChange - lfcSE)) +
      geom_hline(yintercept = 0) +
      xlab("") + ylab("Past-Recent Log2 Fold Change") + labs(title = INVENTORY) +
      scale_fill_manual(values = c("Significant" = "grey30", "Non-significant" = "grey")) +
      coord_flip() +
      theme_bw() +
      theme(legend.position = "none",
            axis.title = element_text(size = 10))
  }))



