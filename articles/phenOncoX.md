# Getting started

  
  

## Installation

``` r
if (!("remotes" %in% installed.packages()) {
  install.packages("remotes")
}
remotes::install_github('sigven/phenOncoX')
```

  
  

## Get OncoTree terms

This shows how to retrieve cancer phenotype terms in, as defined in
[OncoTree](http://oncotree.mskcc.org/#/home) (*max tree depth* = 2)

``` r

## load the data
download_dir <- tempdir()

oncotree <- phenOncoX::get_tree(
  cache_dir = download_dir, max_tree_depth = 2)

## Number of records
nrow(oncotree$records)
#> [1] 249

## Show metadata for underlying resources
oncotree$metadata
#>     source                                    source_description
#> 1 OncoTree A cancer classification system for precision oncology
#>                         source_url
#> 1 http://oncotree.mskcc.org/#/home
#>                                         source_citation source_version
#> 1 Kundra et al., JCO Clin Cancer Inform, 2021; 33625877     2025_10_03
#>   source_abbreviation source_license
#> 1            oncotree      CC BY 4.0
#>                             source_license_url
#> 1 https://creativecommons.org/licenses/by/4.0/
```

  
  

### Get all (OncoTree-expanded) cancer phenotype terms

``` r

## get all cancer phenotype terms

oncoterms <- phenOncoX::get_terms(
  cache_dir = download_dir)

## Number of records
nrow(oncoterms$records)
#> [1] 25492
```

  
  

### Number of (OncoTree-expanded) cancer phenotype terms per primary tumor type/tissue

``` r

## get all cancer phenotype terms

oncoterms <- phenOncoX::get_terms(
  cache_dir = download_dir)

## Number of records
as.data.frame(oncoterms$records |>
  dplyr::filter(!is.na(primary_site)) |>
  dplyr::group_by(primary_site) |> 
  dplyr::summarise(num_terms = dplyr::n(),
                   .groups = "drop") |>
  dplyr::arrange(dplyr::desc(num_terms)))
#>                 primary_site num_terms
#> 1                   Lymphoid      4113
#> 2                Soft Tissue      2799
#> 3              Head and Neck      2331
#> 4                  CNS/Brain      2277
#> 5                       Skin      1820
#> 6                    Myeloid      1772
#> 7                       Lung      1142
#> 8               Colon/Rectum      1009
#> 9       Ovary/Fallopian Tube       885
#> 10                    Breast       764
#> 11         Esophagus/Stomach       739
#> 12                      Bone       491
#> 13     Bladder/Urinary Tract       488
#> 14                    Uterus       455
#> 15                    Kidney       413
#> 16             Biliary Tract       398
#> 17                     Liver       386
#> 18                  Pancreas       347
#> 19                    Cervix       286
#> 20 Peripheral Nervous System       280
#> 21                       Eye       279
#> 22              Vulva/Vagina       260
#> 23                   Thyroid       241
#> 24                    Testis       236
#> 25                  Prostate       221
#> 26             Other/Unknown       171
#> 27                    Pleura       117
#> 28                Peritoneum       101
#> 29                    Thymus        87
#> 30             Adrenal Gland        72
#> 31                     Penis        60
#> 32          Ampulla of Vater        46
```

### Get all (OncoTree-expanded) cancer phenotype terms relevant for prostate cancer

``` r

## get all oncoterms for prostate cancer
oncoterms_prostate <- phenOncoX::get_terms(
  cache_dir = download_dir, 
  site = "Prostate")

## Make as datatable
prostate_terms_table <- DT::datatable(
  dplyr::select(
    oncoterms_prostate$records, 
    primary_site, cui, cui_name,
    dplyr::everything()),
  escape = FALSE,
  extensions = c("Buttons", "Responsive"), 
  width = "100%",
  options = list(
    buttons = c("csv", "excel"), 
    dom = "Bfrtip"))
```

  

  
  

## Session Info

``` r
# set eval = FALSE if you don't want this info (useful for reproducibility) 
# to appear
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] jsonlite_2.0.0    dplyr_1.1.4       compiler_4.5.2    crayon_1.5.3     
#>  [5] tidyselect_1.2.1  phenOncoX_1.1.2   jquerylib_0.1.4   systemfonts_1.3.1
#>  [9] textshaping_1.0.4 yaml_2.3.12       fastmap_1.2.0     R6_2.6.1         
#> [13] generics_0.1.4    curl_7.0.0        knitr_1.51        htmlwidgets_1.6.4
#> [17] tibble_3.3.1      desc_1.4.3        bslib_0.10.0      pillar_1.11.1    
#> [21] rlang_1.1.7       DT_0.34.0         cachem_1.1.0      lgr_0.5.2        
#> [25] xfun_0.56         fs_1.6.6          sass_0.4.10       otel_0.2.0       
#> [29] cli_3.6.5         withr_3.0.2       pkgdown_2.2.0     magrittr_2.0.4   
#> [33] crosstalk_1.2.2   digest_0.6.39     lifecycle_1.0.5   vctrs_0.7.1      
#> [37] evaluate_1.0.5    gargle_1.6.1      glue_1.8.0        ragg_1.5.0       
#> [41] googledrive_2.1.2 httr_1.4.7        rmarkdown_2.30    purrr_1.2.1      
#> [45] tools_4.5.2       pkgconfig_2.0.3   htmltools_0.5.9
```

  
  

## References

Huang, Kuan-Lin, R Jay Mashl, Yige Wu, Deborah I Ritter, Jiayin Wang,
Clara Oh, Marta Paczkowska, et al. 2018. “Pathogenic Germline Variants
in 10,389 Adult Cancers.” *Cell* 173 (2): 355–370.e14.
<http://dx.doi.org/10.1016/j.cell.2018.03.039>.

Kundra, Ritika, Hongxin Zhang, Robert Sheridan, Sahussapont Joseph
Sirintrapun, Avery Wang, Angelica Ochoa, Manda Wilson, et al. 2021.
“OncoTree: A Cancer Classification System for Precision Oncology.” *JCO
Clin Cancer Inform* 5 (February): 221–30.
<http://dx.doi.org/10.1200/CCI.20.00108>.

Louden, Diana Nelson. 2020. “MedGen: NCBI’s Portal to Information on
Medical Conditions with a Genetic Component.” *Med. Ref. Serv. Q.* 39
(2): 183–91. <https://doi.org/10.1080/02763869.2020.1726152>.

Malone, James, Ele Holloway, Tomasz Adamusiak, Misha Kapushesky, Jie
Zheng, Nikolay Kolesnikov, Anna Zhukova, Alvis Brazma, and Helen
Parkinson. 2010. “Modeling Sample Variables with an Experimental Factor
Ontology.” *Bioinformatics* 26 (8): 1112–18.
<http://dx.doi.org/10.1093/bioinformatics/btq099>.

Schriml, Lynn Marie, Cesar Arze, Suvarna Nadendla, Yu-Wei Wayne Chang,
Mark Mazaitis, Victor Felix, Gang Feng, and Warren Alden Kibbe. 2012.
“Disease Ontology: A Backbone for Disease Semantic Integration.”
*Nucleic Acids Res.* 40 (Database issue): D940–6.
<http://dx.doi.org/10.1093/nar/gkr972>.
