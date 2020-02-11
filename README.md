<img src="logo_tb.png" width="220" align="right"/>

# Code and data for the paper

This repository provides R code and data to reproduce results and figures from the paper:

**Assessing responses of micro-eukaryotes biodiversity to the Great Acceleration: a paleolimnological view based on sedimentary DNA.**

François Keck, Laurent Millet, Didier Debroas, David Etienne, Didier Galop and Isabelle Domaizon


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

