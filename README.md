
<!-- README.md is generated from README.Rmd. Please edit that file -->
FungalSporeAnalysis
===================

Load feature file and data (bw files stored as GRanges object).

``` r
devtools::load_all()
gene_list <- readr::read_delim("data/An_Spore_Pol2.txt", delim="\t", col_names = FALSE)
feature_gr <- AnnotationDbi::loadDb("R/an_feature_file_s10_m04_r07.sqlite")
```

Profiles of pre-initiation complex factors ie. RNAP-II, TBP and TFII-B in spores of *A. nidulans* generated using `genelist_specific_profileplot` function.

``` r
1. RNAP-II 
genelist_specific_profileplot(feature_gr=feature_gr,bw_files = "pol2_veA_wt_spore", genelist=gene_list, output_name="An_Spore_Pol2", ymin=3,max_key = 10, min_key = 0, ymax = 5.5, palette = "white_red")

2. TBP
genelist_specific_profileplot(feature_gr=feature_gr,bw_files = "TBP_veA_wt_spore", genelist=gene_list,max_key=4.5,min_key = 0, output_name="An_Spore_TBP", ymin=3, palette = "white_green", ymax = 5.5)

3. TFII-B
genelist_specific_profileplot(feature_gr=feature_gr,bw_files = "TFIIB_veA_wt_spore", genelist=gene_list,max_key=5,min_key = 1, output_name="An_Spore_TFIIB", ymin=3, palette = "white_blue", ymax = 5.5)
```

<img src="plots/An_Spore_Pol2_1030_hm.png" style="width:50.0%" /> <img src="plots/An_Spore_TBP_1030_hm.png" style="width:50.0%" /> <img src="plots/An_Spore_TFIIB_1030_hm.png" style="width:50.0%" />
