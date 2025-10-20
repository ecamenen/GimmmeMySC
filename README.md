
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GimmmeMySC

<!-- badges: start -->

<!-- badges: end -->

#### Version: 1.0.0

#### Author:

Etienne CAMENEN (<etienne.camenen@inserm.fr>)

#### Abstract:

End-to-end workflow for the analysis of single-cell and single-nucleus
omics data, including scRNA-seq and scATAC-seq. Features include
rigorous quality control, ambient RNA correction, and advanced detection
of doublets and duplicated cells. Offers tools for normalization,
feature selection, dimensionality reduction, clustering, and robust
batch effect correction. Supports cell-type annotation using reference
atlases and marker-based methods, sex-specific marker analysis, and cell
cycle scoring. Enables multi-omic integration, pseudo-bulk analysis,
differential expression and accessibility testing, peak calling, motif
enrichment, and gene regulatory network reconstruction. Includes
extensive visualization capabilities such as UMAP/t-SNE embeddings,
violin and dot plots, coverage maps, co-accessibility networks, and
interactive heatmaps.

#### Key-words:

scRNA-seq; scATAC-seq; Multi-omics; Cell annotation; Quality control;
Cell trajectory; Clustering analysis; Epigenetic; Chromatin
accessibility

## Installation

You can install the development version of GimmmeMySC from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("ecamenen/GimmmeMySC")
```

## Example

``` r
library(GimmmeMySC)
#> Le chargement a nécessité le package : tidyverse
#> ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.1.4     ✔ readr     2.1.5
#> ✔ forcats   1.0.0     ✔ stringr   1.5.1
#> ✔ ggplot2   3.5.2     ✔ tibble    3.3.0
#> ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
#> ✔ purrr     1.1.0     
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
```
