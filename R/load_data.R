
#+ Begins
library(tidyverse)
library(cowplot)
library(magrittr)
library(vegan)
library(rpart)
library(DESeq2)
library(flexitarian)

# nest API was broken with tibble v 1.0. Run the following if >= 1.0
nest <- nest_legacy
unnest <- unnest_legacy

setwd("/home/francois/Google Drive/Sync work/projects/top_bottom/")

lakes_meta <- read_csv("data/Lakes_meta.csv") %>% 
  mutate(LAC_POS = paste(LAC, POSITION, sep = "_"))

lakes_troph <- read_csv("data/lakes_trophic_status.csv") %>% 
  mutate(Trophic_Status = str_to_sentence(Trophic_Status) %>% 
           fct_explicit_na(na_level = "Unknown"))

lakes_cn <- read_csv("data/synthese_cn.csv") %>% 
  select(-ID_CAROTTE, -CAROTTE_PROF_MIN, -CAROTTE_PROF_MAX)

lakes_bv <- read_csv("data/synthese_bv.csv")

lakes_dna <- read_csv("data/synthese_dna.csv") %>% 
  group_by(LAC, POSITION) %>% 
  summarise_if(is.numeric, mean)

lakes_kml <- rgdal::readOGR("data/TOP-BOTTOM_48_FK.kml", verbose = FALSE)
lakes_kml$Name <- str_replace(lakes_kml$Name, pattern = "\n", "")

otu_taxo_pr2 <- read_csv("data/OTU_affil_PR2SOUL_MIX_75.csv")
otu_taxo_adl_conv <- read_csv("data/OTU_affil_PR2SOUL_MIX_75_ADLCONV_EDISA.csv") %>%
  select(-COUNT)

lakes_otu <- read_csv("data/OTU_distribution_tax.csv")

otu_meta <- select(lakes_otu,
                   OTU_ID,
                   OTU_SEED) %>%
  select(-matches("^CLASSIF_[0-9]")) %>% 
  left_join(otu_taxo_pr2, by = "OTU_ID") %>% 
  left_join(otu_taxo_adl_conv)


otu_meta <- otu_meta %>% 
  mutate(ADL_RANK_12_UNID =
           select(otu_meta, starts_with("ADL_")) %>% 
           apply(1, function(x) ifelse(all(is.na(x)), NA, x[max(which(!is.na(x)))]))
  )

otu_meta <- otu_meta %>% 
  mutate(ADL_RANK_3_UNID = case_when(
    !is.na(ADL_RANK_3) ~ ADL_RANK_3,
    is.na(ADL_RANK_3) & !is.na(ADL_RANK_2) ~ paste("unidentified", ADL_RANK_2),
    is.na(ADL_RANK_3) & is.na(ADL_RANK_2) ~ paste("unidentified", ADL_RANK_1)
    ))


otu_tbl <- select(lakes_otu,
                  -OTU_ID,
                  -OTU_SEED) %>%
  as.data.frame() %>% t() %>% 
  magrittr::set_colnames(lakes_otu$OTU_ID) %>%
  tidy_cdm("ID_SAMPLE", "OTU_ID", "OTU_COUNT") %>% 
  left_join(lakes_meta) %>% 
  filter(POSITION != "Deep")

info_otu <- otu_tbl %>% 
  group_by(OTU_ID) %>% 
  summarise(TOTAL_COUNT = sum(OTU_COUNT)) %>% 
  arrange(TOTAL_COUNT) %>% 
  mutate(CUMSUM_COUNT = cumsum(TOTAL_COUNT),
         QUANTILE = 1:nrow(.)/nrow(.) * 100)

otu_tbl_replicates <- otu_tbl

# MERGE DATA (SUM and SUBSAMPLE OF REPLICATES)
samples_depth <- otu_tbl %>%
  group_by(LAC, POSITION) %>%
  summarize(DEPTH = sum(OTU_COUNT)) %>%
  unite(LAC_POS, LAC, POSITION, remove = FALSE) %>%
  left_join(
    lakes_meta %>% group_by(LAC, POSITION) %>% summarise(n = n())
  ) %>%
  mutate(MEAN_DEPTH = round(DEPTH/n))

otu_merge_tbl <- otu_tbl %>%
  group_by(LAC, POSITION, OTU_ID) %>%
  summarize(OTU_COUNT = sum(OTU_COUNT)) %>%
  unite(LAC_POS, LAC, POSITION, remove = FALSE) %>%
  spread_cdm(LAC_POS, OTU_ID, OTU_COUNT) %>%
  rrarefy(sample = samples_depth$MEAN_DEPTH) %>%
  tidy_cdm("LAC_POS", "OTU_ID", "OTU_COUNT") %>%
  separate(LAC_POS, into = c("LAC", "POSITION"), sep = "_", remove = FALSE)

rm(samples_depth)

lakes_meta_nr <- lakes_meta %>% 
  filter(REPLICATE == 1) %>% 
  select(-ID_SAMPLE, -REPLICATE)


# otu_merge_tbl <- otu_tbl_replicates %>% 
#   group_by(LAC, POSITION, OTU_ID) %>% 
#   mutate(INTERSECT = (function(x){
#     ifelse(x[1] > 0 & x[2] > 0, TRUE, FALSE)
#   })(OTU_COUNT)) %>% 
#   filter(INTERSECT == TRUE) %>% 
#   summarize(OTU_COUNT = round(mean(OTU_COUNT))) %>% 
#   unite(LAC_POS, LAC, POSITION, remove = FALSE) %>% 
#   spread_cdm(LAC_POS, OTU_ID, OTU_COUNT, fill.missing = 0) %>% 
#   tidy_cdm("LAC_POS", "OTU_ID", "OTU_COUNT") %>% 
#   separate(LAC_POS, into = c("LAC", "POSITION"), sep = "_", remove = FALSE)


# FILTER OTUs WITH < 10 OCCURENCES IN THE DATASET
otu_merge_tbl <- otu_merge_tbl %>%
  left_join(info_otu) %>% 
  filter(TOTAL_COUNT > 10)

#### TAXONOMIC FILTER: Unid Eukaryota, Metazoa and Embryophyta
otu_meta <- otu_meta %>% 
  filter(!is.na(Supergroup)) %>% 
  filter(ADL_RANK_4 != "Embryophyta" | is.na(ADL_RANK_4)) %>% 
  filter(ADL_RANK_2 != "Metazoa" | is.na(ADL_RANK_2))

otu_merge_tbl <- otu_merge_tbl %>% 
  filter(OTU_ID %in% otu_meta$OTU_ID)

# NORMALIZATION WITH PROPORTIONS
otu_rar_tbl <- otu_merge_tbl %>%
  group_by(LAC_POS) %>% 
  mutate(COUNT = OTU_COUNT/sum(OTU_COUNT)) %>% 
  select(-OTU_COUNT) %>% 
  dplyr::rename(TAXON = OTU_ID)

# TAXONOMY
tax_rar_rank_1 <- otu_rar_tbl %>% 
  left_join(otu_meta, by = c("TAXON" = "OTU_ID")) %>% 
  group_by(LAC_POS, ADL_RANK_1) %>% 
  summarise(COUNT = sum(COUNT)) %>%
  dplyr::rename(TAXON = ADL_RANK_1) %>% 
  separate(LAC_POS, sep = "_", into = c("LAC", "POSITION"), remove = FALSE)  

tax_rar_rank_2 <- otu_rar_tbl %>% 
  left_join(otu_meta, by = c("TAXON" = "OTU_ID")) %>% 
  group_by(LAC_POS, ADL_RANK_2) %>% 
  summarise(COUNT = sum(COUNT)) %>%
  dplyr::rename(TAXON = ADL_RANK_2) %>% 
  separate(LAC_POS, sep = "_", into = c("LAC", "POSITION"), remove = FALSE)  

tax_rar_rank_3 <- otu_rar_tbl %>% 
  left_join(otu_meta, by = c("TAXON" = "OTU_ID")) %>% 
  group_by(LAC_POS, ADL_RANK_3) %>% 
  summarise(COUNT = sum(COUNT)) %>%
  dplyr::rename(TAXON = ADL_RANK_3) %>% 
  separate(LAC_POS, sep = "_", into = c("LAC", "POSITION"), remove = FALSE)  

tax_rar_rank_3_unid <- otu_rar_tbl %>% 
  left_join(otu_meta, by = c("TAXON" = "OTU_ID")) %>% 
  group_by(LAC_POS, ADL_RANK_3_UNID) %>% 
  summarise(COUNT = sum(COUNT)) %>%
  dplyr::rename(TAXON = ADL_RANK_3_UNID) %>% 
  separate(LAC_POS, sep = "_", into = c("LAC", "POSITION"), remove = FALSE)  

tax_rar_rank_12_unid <- otu_rar_tbl %>% 
  left_join(otu_meta, by = c("TAXON" = "OTU_ID")) %>% 
  group_by(LAC_POS, ADL_RANK_12_UNID) %>% 
  summarise(COUNT = sum(COUNT)) %>%
  dplyr::rename(TAXON = ADL_RANK_12_UNID) %>% 
  separate(LAC_POS, sep = "_", into = c("LAC", "POSITION"), remove = FALSE)  


# TROPHIC GROUPS
troph_mode_rar <- otu_rar_tbl %>% 
  left_join(otu_meta, by = c("TAXON" = "OTU_ID")) %>% 
  group_by(LAC_POS, TROPHIC_MODE) %>% 
  summarise(COUNT = sum(COUNT)) %>%
  dplyr::rename(TAXON = TROPHIC_MODE) %>% 
  separate(LAC_POS, sep = "_", into = c("LAC", "POSITION"), remove = FALSE)  

troph_type_rar <- otu_rar_tbl %>% 
  left_join(otu_meta, by = c("TAXON" = "OTU_ID")) %>% 
  separate_rows(TROPHIC_TYPE, sep = ";") %>% 
  group_by(LAC_POS, TROPHIC_TYPE) %>% 
  summarise(COUNT = sum(COUNT)) %>%
  dplyr::rename(TAXON = TROPHIC_TYPE) %>% 
  separate(LAC_POS, sep = "_", into = c("LAC", "POSITION"), remove = FALSE)

# SUB MATRIX OF TROPHIC GROUPS
otu_rar_tbl_heterotrophs <- otu_rar_tbl %>% 
  left_join(otu_meta, by = c("TAXON" = "OTU_ID")) %>% 
  filter(TROPHIC_MODE == "Heterotrophs") %>% 
  select(LAC_POS, LAC, POSITION, TAXON, COUNT)

otu_rar_tbl_autotrophs <- otu_rar_tbl %>% 
  left_join(otu_meta, by = c("TAXON" = "OTU_ID")) %>% 
  filter(TROPHIC_MODE == "Autotrophs") %>% 
  select(LAC_POS, LAC, POSITION, TAXON, COUNT)

otu_rar_tbl_mixotrophs <- otu_rar_tbl %>% 
  left_join(otu_meta, by = c("TAXON" = "OTU_ID")) %>% 
  filter(TROPHIC_MODE == "Mixotrophs") %>% 
  select(LAC_POS, LAC, POSITION, TAXON, COUNT)



# LIST OF INVENTORIES
com_rar_list <- list(otu_rar_tbl = otu_rar_tbl,
                     tax_rar_rank_1 = tax_rar_rank_1, tax_rar_rank_2 = tax_rar_rank_2, tax_rar_rank_3 = tax_rar_rank_3,
                     tax_rar_rank_3_unid = tax_rar_rank_3_unid, tax_rar_rank_12_unid = tax_rar_rank_12_unid,
                     troph_mode_rar = troph_mode_rar, troph_type_rar = troph_type_rar,
                     otu_rar_tbl_heterotrophs = otu_rar_tbl_heterotrophs,
                     otu_rar_tbl_autotrophs = otu_rar_tbl_autotrophs,
                     otu_rar_tbl_mixotrophs = otu_rar_tbl_mixotrophs
                     ) %>% 
  bind_rows(.id = "INVENTORY") %>% 
  group_by(INVENTORY) %>% 
  nest() %>% 
  mutate(INVENTORY = c("OTU",
                       "Taxonomie Rank 1", "Taxonomie Rank 2", "Taxonomie Rank 3",
                       "Taxonomie Rank 3 UNID", "Taxonomie Rank 12 UNID",
                       "Trophic Mode", "Trophic Type",
                       "OTU subset Heterotrophs",
                       "OTU subset Autotrophs",
                       "OTU subset Mixotrophs"
                       ))

rm(otu_rar_tbl,
   tax_rar_rank_1, tax_rar_rank_2, tax_rar_rank_3,
   tax_rar_rank_3_unid, tax_rar_rank_12_unid,
   otu_rar_tbl_heterotrophs, otu_rar_tbl_autotrophs, otu_rar_tbl_mixotrophs,
   troph_mode_rar, troph_type_rar)
rm(lakes_otu)


com_rar_list <- com_rar_list %>% 
  mutate(data = map(data, ~group_by(.x, LAC_POS) %>%
                      mutate(FREQ = COUNT/sum(COUNT) * 100)))

com_rar_list <- com_rar_list %>% 
  mutate(BRAY_DIST = map(data, ~ .x %>% 
                           spread_cdm(LAC_POS, TAXON, COUNT) %>% 
                           vegdist(method = "bray") %>% 
                           broom::tidy() %>% 
                           transmute(LAC_POS_1 = as.character(item1),
                                     LAC_POS_2 = as.character(item2),
                                     BRAY_DIST = distance) %>% 
                           left_join(lakes_meta_nr, by = c("LAC_POS_1" = "LAC_POS")) %>% 
                           left_join(lakes_meta_nr, by = c("LAC_POS_2" = "LAC_POS"),
                                     suffix = c("_1", "_2")) %>% 
                           filter(!is.na(LAC_1), !is.na(LAC_2))
  ))


# TAXONOMY NON RARIFIED

otu_raw_tbl <- otu_merge_tbl %>%
  dplyr::rename(TAXON = OTU_ID , COUNT = OTU_COUNT)

tax_raw_rank_1 <- otu_raw_tbl %>% 
  left_join(otu_meta, by = c("TAXON" = "OTU_ID")) %>% 
  group_by(LAC_POS, ADL_RANK_1) %>% 
  summarise(COUNT = sum(COUNT)) %>%
  dplyr::rename(TAXON = ADL_RANK_1) %>% 
  separate(LAC_POS, sep = "_", into = c("LAC", "POSITION"), remove = FALSE)  

tax_raw_rank_2 <- otu_raw_tbl %>% 
  left_join(otu_meta, by = c("TAXON" = "OTU_ID")) %>% 
  group_by(LAC_POS, ADL_RANK_2) %>% 
  summarise(COUNT = sum(COUNT)) %>%
  dplyr::rename(TAXON = ADL_RANK_2) %>% 
  separate(LAC_POS, sep = "_", into = c("LAC", "POSITION"), remove = FALSE)  

tax_raw_rank_3 <- otu_raw_tbl %>% 
  left_join(otu_meta, by = c("TAXON" = "OTU_ID")) %>% 
  group_by(LAC_POS, ADL_RANK_3) %>% 
  summarise(COUNT = sum(COUNT)) %>%
  dplyr::rename(TAXON = ADL_RANK_3) %>% 
  separate(LAC_POS, sep = "_", into = c("LAC", "POSITION"), remove = FALSE)  

tax_raw_rank_3_unid <- otu_raw_tbl %>% 
  left_join(otu_meta, by = c("TAXON" = "OTU_ID")) %>% 
  group_by(LAC_POS, ADL_RANK_3_UNID) %>% 
  summarise(COUNT = sum(COUNT)) %>%
  dplyr::rename(TAXON = ADL_RANK_3_UNID) %>% 
  separate(LAC_POS, sep = "_", into = c("LAC", "POSITION"), remove = FALSE)  

tax_raw_rank_12_unid <- otu_raw_tbl %>% 
  left_join(otu_meta, by = c("TAXON" = "OTU_ID")) %>% 
  group_by(LAC_POS, ADL_RANK_12_UNID) %>% 
  summarise(COUNT = sum(COUNT)) %>%
  dplyr::rename(TAXON = ADL_RANK_12_UNID) %>% 
  separate(LAC_POS, sep = "_", into = c("LAC", "POSITION"), remove = FALSE)  


# TROPHIC GROUPS
troph_mode_raw <- otu_raw_tbl %>% 
  left_join(otu_meta, by = c("TAXON" = "OTU_ID")) %>% 
  group_by(LAC_POS, TROPHIC_MODE) %>% 
  summarise(COUNT = sum(COUNT)) %>%
  dplyr::rename(TAXON = TROPHIC_MODE) %>% 
  separate(LAC_POS, sep = "_", into = c("LAC", "POSITION"), remove = FALSE)  

troph_type_raw <- otu_raw_tbl %>% 
  left_join(otu_meta, by = c("TAXON" = "OTU_ID")) %>% 
  separate_rows(TROPHIC_TYPE, sep = ";") %>% 
  group_by(LAC_POS, TROPHIC_TYPE) %>% 
  summarise(COUNT = sum(COUNT)) %>%
  dplyr::rename(TAXON = TROPHIC_TYPE) %>% 
  separate(LAC_POS, sep = "_", into = c("LAC", "POSITION"), remove = FALSE)

# SUB MATRIX OF TROPHIC GROUPS
otu_raw_tbl_heterotrophs <- otu_raw_tbl %>% 
  left_join(otu_meta, by = c("TAXON" = "OTU_ID")) %>% 
  filter(TROPHIC_MODE == "Heterotrophs") %>% 
  select(LAC_POS, LAC, POSITION, TAXON, COUNT)

otu_raw_tbl_autotrophs <- otu_raw_tbl %>% 
  left_join(otu_meta, by = c("TAXON" = "OTU_ID")) %>% 
  filter(TROPHIC_MODE == "Autotrophs") %>% 
  select(LAC_POS, LAC, POSITION, TAXON, COUNT)

otu_raw_tbl_mixotrophs <- otu_raw_tbl %>% 
  left_join(otu_meta, by = c("TAXON" = "OTU_ID")) %>% 
  filter(TROPHIC_MODE == "Mixotrophs") %>% 
  select(LAC_POS, LAC, POSITION, TAXON, COUNT)


# LIST OF INVENTORIES
com_merged_raw_list <- list(otu_raw_tbl = otu_raw_tbl,
                     tax_raw_rank_1 = tax_raw_rank_1, tax_raw_rank_2 = tax_raw_rank_2, tax_raw_rank_3 = tax_raw_rank_3,
                     tax_raw_rank_3_unid = tax_raw_rank_3_unid, tax_raw_rank_12_unid = tax_raw_rank_12_unid,
                     troph_mode_raw = troph_mode_raw, troph_type_raw = troph_type_raw,
                     otu_raw_tbl_heterotrophs = otu_raw_tbl_heterotrophs,
                     otu_raw_tbl_autotrophs = otu_raw_tbl_autotrophs,
                     otu_raw_tbl_mixotrophs = otu_raw_tbl_mixotrophs) %>% 
  bind_rows(.id = "INVENTORY") %>% 
  group_by(INVENTORY) %>% 
  nest() %>% 
  mutate(INVENTORY = c("OTU",
                       "Taxonomie Rank 1", "Taxonomie Rank 2", "Taxonomie Rank 3",
                       "Taxonomie Rank 3 UNID", "Taxonomie Rank 12 UNID",
                       "Trophic Mode", "Trophic Type",
                       "OTU subset Heterotrophs",
                       "OTU subset Autotrophs",
                       "OTU subset Mixotrophs"))

rm(otu_raw_tbl,
   tax_raw_rank_1, tax_raw_rank_2, tax_raw_rank_3,
   tax_raw_rank_3_unid, tax_raw_rank_12_unid,
   otu_raw_tbl_heterotrophs, otu_raw_tbl_autotrophs, otu_raw_tbl_mixotrophs,
   troph_mode_raw, troph_type_raw)
