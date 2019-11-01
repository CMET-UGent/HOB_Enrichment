# HOB_Enrichment
Repository for data analyses scripts of Hu *et al.* (2019)

## Files

Script Filename | Content 
---------|----------
Heatmaps_and_Ordinations.R | main script for drawing multiple heatmaps and ordination plots
DataLoader.R | separate script to preload and wrangle the data for DE_testing.R
DE_testing.R | DESeq2 differential abundance testing, knitr spin compatible

## Dependencies

All data was analyzed in R (v 3.6.1.). The following additional packages should
be installed to run the provided code:

- ggplot2         (CRAN)
- extrafont       (CRAN)
- phyloseq        (BioConductor)
- NMF             (CRAN)
- vegan           (CRAN)
- gplots          (CRAN)   
- readxl          (CRAN)
- DESeq2          (BioConductor)
- data.table      (CRAN)
- splitstackshape (CRAN)
- dplyr           (CRAN)
- scales          (CRAN)
- ape             (CRAN)