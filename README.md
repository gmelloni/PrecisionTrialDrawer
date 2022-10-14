<!-- badges: start -->
[![R-CMD-check](https://github.com/gmelloni/PrecisionTrialDrawer/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gmelloni/PrecisionTrialDrawer/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# PrecisionTrialDrawer

A set of tools to help bioinformaticians, clinicians and biostatisticians
to design, analyze and finalize a custom gene panel for clinical trial in cancer genomics.

## Install

You can install **PrecisionTrialDrawer** using devtools

```{r}
install.packages("devtools")
devtools::install_github("gmelloni/PrecisionTrialDrawer")
```
In case some of the dependencies are not automatically downloaded,
you can install the dependencies first and then run *install_github*

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite(c('cBioPortalData','parallel','stringr',
    'reshape2','data.table','RColorBrewer',
    'BiocParallel','magrittr','biomaRt',
    'XML','httr','jsonlite',
    'ggplot2','ggrepel','grid',
    'S4Vectors','IRanges','GenomicRanges',
    'googleVis','shiny',
    'shinyBS','DT','brglm','BiocStyle',
    'knitr','rmarkdown','dplyr','testthat'))
biocLite('LowMACAAnnotation')
devtools::install_github("gmelloni/PrecisionTrialDrawer")
```
## Info

For an explanation of the package functionalities, refer to the official website: https://gmelloni.github.io/ptd/

## Web resource

A demo version of PTD is available as a web resource [here](https://gmelloni.github.io/ptd/shinyapp.html). The web version is meant for demonstration purpose and only a limited group of functionalities are implemented. 
