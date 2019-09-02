

# Fraction of intersect OTU Top and Bottom by lake
frac_taxa_temp <- com_merged_raw_list %>% 
  mutate(fractions = map(data, function(data){
    data %>% 
    filter(COUNT > 0) %>%
    split(.$LAC) %>% 
    map_dfr(.f = function(x){
      bot <- x %>% filter(POSITION == "Bottom") %>% .$TAXON
      top <- x %>% filter(POSITION == "Top") %>% .$TAXON
      c(length(bot), length(top), length(intersect(top, bot)), length(setdiff(top, bot)), length(setdiff(bot, top)))
    }) %>% 
    mutate(FRACTION = c("Bottom all", "Top all", "Intersection", "Top only", "Bottom only")) %>% 
    gather(key = "LAC", value = "N", -FRACTION) %>%
    filter(FRACTION %in% c("Intersection", "Top only", "Bottom only"))
  })) %>% 
  mutate(plots = map(fractions, function(fractions){
    fractions %>% 
      ggplot() +
      geom_col(aes(factor(LAC, levels = rev(levels(factor(LAC)))), N, fill = FRACTION),
               position = "fill") +
      coord_flip() +
      scale_y_continuous(labels = scales::percent) +
      theme_minimal() +
      labs(fill = "Fraction") + ylab("Relative proportion") + xlab("Lac")
  })) %>% 
  mutate(mean_fractions = map(fractions, function(fractions){
    fractions %>% 
      group_by(FRACTION) %>% 
      summarise(mean = mean(N))
  }))



# Fraction of intersect OTU Top vs Bottom among lakes
frac_taxa_spatial <- com_merged_raw_list %>% 
  mutate(OTU_lists = map(data, function(data){
    data %>% 
    filter(COUNT > 0) %>%
    split(.$LAC) %>% 
    map(.f = function(x){
      bot <- x %>% filter(POSITION == "Bottom") %>% .$TAXON
      top <- x %>% filter(POSITION == "Top") %>% .$TAXON
      list(Bottom = bot, Top = top)
    })
  })) %>% 
  mutate(plots = map(OTU_lists, function(OTU_lists){
    crossing(LAC_1 = names(OTU_lists), LAC_2 = names(OTU_lists)) %>% 
      filter(LAC_1 < LAC_2) %>% 
      mutate(RES = map2(LAC_1, LAC_2, function(LAC_1, LAC_2){
        l1b <- OTU_lists[[LAC_1]]$Bottom
        l2b <- OTU_lists[[LAC_2]]$Bottom
        l1t <- OTU_lists[[LAC_1]]$Top
        l2t <- OTU_lists[[LAC_2]]$Top
        tibble(POSITION = c("Bottom", "Top"),
               Intersection = c(length(intersect(l1b, l2b)), length(intersect(l1t, l2t))),
               Only_LAC_1 = c(length(setdiff(l1b, l2b)), length(setdiff(l1t, l2t))),
               Only_Lac_2 = c(length(setdiff(l2b, l1b)), length(setdiff(l2t, l1t)))
        )
      })) %>% 
      unnest(RES) %>% 
      mutate(Freq_Intersection = Intersection / (Intersection + Only_LAC_1 + Only_Lac_2),
             Freq_Unique = (Only_LAC_1 + Only_Lac_2) / (Intersection + Only_LAC_1 + Only_Lac_2)) %>% 
      mutate(jitter_val = runif(nrow(.), -0.1, 0.1)) %>% 
      ggplot() +
      geom_boxplot(aes(as.numeric(as.factor(POSITION)), Freq_Intersection, group = POSITION, fill = POSITION), outlier.shape = NA)+
      geom_point(aes(as.numeric(as.factor(POSITION)) + jitter_val, Freq_Intersection)) +
      scale_x_continuous(breaks = c(1, 2), labels = c("Bottom", "Top")) +
      theme(legend.position = "none") +
      xlab("") + ylab("Taxa intersection (%)")
  }))



OTU_all <- com_merged_raw_list %>% 
  filter(INVENTORY == "OTU") %>% 
  .$data %>% .[[1]] %>% 
  filter(COUNT > 0)

OTU_Bottom <- OTU_all %>% 
  filter(POSITION == "Bottom")

OTU_Top <- OTU_all %>% 
  filter(POSITION == "Top")

OTU_Bottom_only <- OTU_Bottom %>% 
  filter(TAXON %in% setdiff(OTU_Bottom %>% .$TAXON %>% unique(), OTU_Top %>% .$TAXON %>% unique()))

OTU_Top_only <- OTU_Top %>% 
  filter(TAXON %in% setdiff(OTU_Top %>% .$TAXON %>% unique(), OTU_Bottom %>% .$TAXON %>% unique()))

# OTU_Bottom_only %>% 
#   group_by(TAXON) %>% 
#   summarise(Number_of_lakes = n(), Total_count = sum(COUNT)) %>% 
#   left_join(otu_meta[, c("OTU_ID", "ADL_RANK_3", "ADL_RANK_12_UNID")], by = c("TAXON" = "OTU_ID")) %>% 
#   write_csv("OTU_Bottom_only.csv")
# 
# OTU_Top_only %>% 
#   group_by(TAXON) %>% 
#   summarise(Number_of_lakes = n(), Total_count = sum(COUNT)) %>% 
#   left_join(otu_meta[, c("OTU_ID", "ADL_RANK_3", "ADL_RANK_12_UNID")], by = c("TAXON" = "OTU_ID")) %>% 
#   write_csv("OTU_Top_only.csv")
# 
# OTU_all %>% 
#   left_join(otu_meta[, c("OTU_ID", "ADL_RANK_3", "ADL_RANK_12_UNID")], by = c("TAXON" = "OTU_ID")) %>% 
#   group_by(LAC, POSITION) %>% 
#   summarise(list_tax = list(unique(ADL_RANK_3))) %>% 
#   group_by(LAC) %>% 
#   summarise(zz = list((function(x){tibble(Bottom_only = paste(setdiff(x[[1]], x[[2]]), collapse = ", "),
#                                         Top_only = paste(setdiff(x[[2]], x[[1]]), collapse = ", "))})(list_tax))) %>% 
#   unnest(zz) %>% 
#   gather(key = "Position", value = "Taxa", -LAC) %>% 
#   write_csv2("Taxon_TB_only.csv")
