# Code and data for the paper

This repository provides R code and data to reproduce results and figures from the paper:

> XXX [DOI URL](DOI URL)

## Analyses

All the analyses can be ran from the master script `RUN_Main.R`.
A report of the main results and figures can be generated from `results_report.Rmd`. Supplementary materials can be generated from ...

## Dependencies

Install R packages from CRAN:

    install.packages(c("tidyverse", "cowplot", "magrittr", "broom",
                       "vegan", "rpart", rmarkdown))

Install R packages from GitHub (you need devtools):

    devtools::install_github("fkeck/flexitarian")

Install R packages from Bioconductor (you need BiocManager):

    BiocManager::install("DESeq2")


## R Session Info