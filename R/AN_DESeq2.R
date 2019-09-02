
#### ANOVA 1 TB ####
deseq_TB_AOV1 <- com_merged_raw_list %>% 
  filter(!str_detect(INVENTORY, "OTU")) %>% 
  mutate(results = map(data, function(data){
    tab <- data %>% 
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
      geom_hline(yintercept = 0) +
      geom_col(aes(fct_reorder(Taxa, log2FoldChange), log2FoldChange, fill = Wald_test)) +
      geom_linerange(aes(x = fct_reorder(Taxa, log2FoldChange),
                         ymin = log2FoldChange + lfcSE,
                         ymax = log2FoldChange - lfcSE)) +
      xlab("") + ylab("Log2 Fold Change") + labs(title = INVENTORY) +
      theme(#axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1),
            legend.position = "none") +
      scale_fill_manual(values = c("Significant" = "grey30", "Non-significant" = "grey")) +
      coord_flip()
  }))



