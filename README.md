<img src="logo_tb.png" width="220" align="right"/>

# Code and data for the paper

[![DOI](https://zenodo.org/badge/236010392.svg)](https://zenodo.org/badge/latestdoi/236010392)

This repository provides R code and data to reproduce results and figures from the paper:

Keck, F. et al. Assessing the response of micro-eukaryotic diversity to the Great Acceleration using lake sedimentary DNA. Nature Communications 11, 3831 (2020).


## Analyses

Scripts can be found in the `/R` directory. They are organized in separated modules but all the analyses can be ran from the master script `RUN_Main.R`.
Once the analyses are completed, a report of the main results and figures can be generated from `Rmd/results_report.Rmd`.
Supplementary materials can be generated from `Rmd/supplementary_info.Rmd`.

## Dependencies

Install R packages from CRAN:

    install.packages(c("tidyverse", "cowplot", "broom", "ggpol",
                       "vegan", "rpart", "rmarkdown", "rgdal",
                       "knitr", "kableExtra", "ape", "phytools"))

Install R packages from GitHub (you need devtools):

    devtools::install_github("fkeck/flexitarian")

Install R packages from Bioconductor (you need BiocManager):

    BiocManager::install("DESeq2")

