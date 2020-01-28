
set.seed(999)
source("R/load_data.R", echo = TRUE)

source("R/AN_Richness.R")
source("R/AN_NMDS.R")
source("R/AN_BrayCurtis_TB.R")
source("R/AN_Stacked_Delta_barplots.R")
source("R/AN_PERMANOVA_paired.R")
source("R/AN_PERMDISP_paired.R")
source("R/AN_DESeq2.R")
source("R/AN_Sediment_basics_TB.R")
source("R/AN_Fraction_intersect.R")
source("R/AN_General_stats.R")



rmarkdown::render("Rmd/results_report.Rmd", output_dir = "results", output_file = "results_report.html")
browseURL("results/results_report.html")

rmarkdown::render("Rmd/supplementary_info.Rmd", output_dir = "results", output_file = "supplementary_info.pdf")
browseURL("results/supplementary_info.pdf")

# Supplementary Dataset 1
tibble(OTU_ID = com_rar_list$data[[1]]$TAXON %>% unique()) %>% 
  left_join(otu_meta) %>% 
  select(OTU_ID, DNA_SEQ, starts_with("ADL_RANK"), TROPHIC_MODE, TROPHIC_TYPE) %>% 
  select(-ends_with("_UNID")) %>% 
  write_csv("results/SI_Dataset_1_OTU_meta.csv")



save.image("results/session_image.RData")
writeLines(capture.output(sessionInfo()), "results/session_info.txt")

# load("results/session_image.RData")

