# Getting started

  
  

### Installation

``` r

if (!("remotes" %in% installed.packages())) {
 install.packages("remotes")
}
remotes::install_github('sigven/phenOncoX')
```

  
  

### Core OncoTree terms

This shows how to retrieve cancer phenotype terms in, as defined in
[OncoTree](http://oncotree.mskcc.org/#/home) (*max tree depth* = 2)

``` r

## load the data
download_dir <- tempdir()

oncotree <- phenOncoX::get_tree(
  cache_dir = download_dir, max_tree_depth = 2)

## Number of records
nrow(oncotree$records)
```

    [1] 265

``` r

## Show metadata for underlying resources
oncotree$metadata
```

        source                                    source_description
    1 OncoTree A cancer classification system for precision oncology
                            source_url
    1 http://oncotree.mskcc.org/#/home
                                            source_citation source_version
    1 Kundra et al., JCO Clin Cancer Inform, 2021; 33625877     2025_10_03
      source_abbreviation source_license
    1            oncotree      CC BY 4.0
                                source_license_url
    1 https://creativecommons.org/licenses/by/4.0/

  
  

### Cancer phenotype terms

``` r

## get all cancer phenotype terms - OncoTree-expanded 
oncoterms <- phenOncoX::get_terms(
  cache_dir = download_dir)

## Number of records
nrow(oncoterms$records)
```

    [1] 27694

  
  

### Term statistics per primary tumor type/tissue

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
```

                    primary_site num_terms
    1                   Lymphoid      4362
    2                Soft Tissue      3061
    3              Head and Neck      2608
    4                  CNS/Brain      2475
    5                    Myeloid      1955
    6                       Skin      1885
    7                       Lung      1197
    8               Colon/Rectum      1122
    9       Ovary/Fallopian Tube       938
    10                    Breast       821
    11         Esophagus/Stomach       758
    12                      Bone       563
    13     Bladder/Urinary Tract       497
    14                    Uterus       486
    15                     Liver       478
    16                    Kidney       474
    17 Peripheral Nervous System       469
    18                  Pancreas       368
    19             Biliary Tract       357
    20                    Cervix       288
    21                       Eye       274
    22              Vulva/Vagina       269
    23                   Thyroid       259
    24                    Testis       243
    25                  Prostate       227
    26             Other/Unknown       187
    27                    Pleura       135
    28                Peritoneum       110
    29                    Thymus        95
    30             Adrenal Gland        94
    31                     Penis        64
    32          Ampulla of Vater        47

### Terms relevant for prostate cancer

``` r

## get all oncoterms for prostate cancer
oncoterms_prostate <- phenOncoX::get_terms(
  cache_dir = download_dir, 
  site = "Prostate")

# ## Make as datatable
# prostate_terms_table <- DT::datatable(
#   dplyr::select(
#     oncoterms_prostate$records, 
#     primary_site, cui, cui_name,
#     dplyr::everything()),
#   escape = FALSE,
#   extensions = c("Buttons", "Responsive"), 
#   width = "100%",
#   options = list(
#     buttons = c("csv", "excel"), 
#     dom = "Bfrtip"))

prostate_terms_table <- reactable::reactable(
  dplyr::select(
    oncoterms_prostate$records, 
    primary_site, cui, cui_name,
    dplyr::everything()),
  filterable = TRUE,
  sortable = TRUE,
  striped = TRUE,
  compact = TRUE,
  searchable = TRUE,
  defaultPageSize = 10
  #downloadable = TRUE
)
```

  

  
  

## Session Info

``` r

sessionInfo()
```

    R version 4.6.1 (2026-06-24)
    Platform: x86_64-pc-linux-gnu
    Running under: Ubuntu 24.04.4 LTS

    Matrix products: default
    BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
    LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0

    locale:
     [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8
     [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8
     [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C
    [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C

    time zone: UTC
    tzcode source: system (glibc)

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base

    loaded via a namespace (and not attached):
     [1] jsonlite_2.0.0    dplyr_1.2.1       compiler_4.6.1    crayon_1.5.3
     [5] tidyselect_1.2.1  phenOncoX_1.2.3   yaml_2.3.12       fastmap_1.2.0
     [9] R6_2.6.1          generics_0.1.4    curl_7.1.0        knitr_1.51
    [13] htmlwidgets_1.6.4 tibble_3.3.1      reactable_0.4.5   pillar_1.11.1
    [17] rlang_1.3.0       lgr_0.5.2         reactR_0.6.1      xfun_0.59
    [21] fs_2.1.0          otel_0.2.0        cli_3.6.6         withr_3.0.3
    [25] magrittr_2.0.5    crosstalk_1.2.2   digest_0.6.39     lifecycle_1.0.5
    [29] vctrs_0.7.3       evaluate_1.0.5    gargle_1.6.1      glue_1.8.1
    [33] googledrive_2.1.2 rmarkdown_2.31    purrr_1.2.2       httr_1.4.8
    [37] tools_4.6.1       pkgconfig_2.0.3   htmltools_0.5.9  

  
  

## References

Huang, Kuan-Lin, R Jay Mashl, Yige Wu, et al. 2018. “Pathogenic Germline
Variants in 10,389 Adult Cancers.” *Cell* 173 (2): 355–370.e14.
<http://dx.doi.org/10.1016/j.cell.2018.03.039>.

Kundra, Ritika, Hongxin Zhang, Robert Sheridan, et al. 2021. “OncoTree:
A Cancer Classification System for Precision Oncology.” *JCO Clin Cancer
Inform* 5 (February): 221–30. <http://dx.doi.org/10.1200/CCI.20.00108>.

Louden, Diana Nelson. 2020. “MedGen: NCBI’s Portal to Information on
Medical Conditions with a Genetic Component.” *Med. Ref. Serv. Q.* 39
(2): 183–91. <https://doi.org/10.1080/02763869.2020.1726152>.

Malone, James, Ele Holloway, Tomasz Adamusiak, et al. 2010. “Modeling
Sample Variables with an Experimental Factor Ontology.” *Bioinformatics*
26 (8): 1112–18. <http://dx.doi.org/10.1093/bioinformatics/btq099>.

Schriml, Lynn Marie, Cesar Arze, Suvarna Nadendla, et al. 2012. “Disease
Ontology: A Backbone for Disease Semantic Integration.” *Nucleic Acids
Res.* 40 (Database issue): D940–6.
<http://dx.doi.org/10.1093/nar/gkr972>.
