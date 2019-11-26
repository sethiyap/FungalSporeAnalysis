
<!-- README.md is generated from README.Rmd. Please edit that file -->
FungalSporeAnalysis
===================

<!-- badges: start -->
[![platform](https://img.shields.io/badge/R-%3E%20v3.5.1-brightgreen)](https://shields.io/category/platform-support)
<!-- badges: end -->
Fungal spore analysis is a repository and package consisting of functions used to analyze genomics data of fungal spores. Following functions have been iteratively used in the analysis;

-   `genelist_specific_profileplot`
-   `profiles_normalized_by_control`
-   `lineplot_from_bw`
-   `ggplot_heatmap`
-   `GO_diamond`
-   `bw_corr`

### Install

    options(repos  = BiocManager::repositories())

    install.packages("devtools")

    devtools::install_github("sethiyap/FungalSporeAnalysis")

For more details on functions and analysis see <https://sethiyap.github.io/FungalSporeAnalysis/>.
