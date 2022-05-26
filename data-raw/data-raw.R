library(magrittr)
source('data-raw/functions.R')

disease_ontology_release <- 'v2022-04-28'
efo_release <- 'v3.42.0'
oncotree_release <- '2021_11_02'

## Get UMLS / DiseaseOntology / EFO mappings
umls_map <- map_umls(
  update = T,
  basedir = here::here())


do_map <- map_disease_ontology(
  skip_non_cui_mapped = T,
  release = disease_ontology_release,
  umls_map = umls_map,
  basedir = here::here())


efo_map <- map_efo(
  umls_map = umls_map,
  efo_release = efo_release,
  update = T,
  basedir = here::here())

## Use OncoTree as starting point for phenotype cross-map
onco_map <- onco_pheno_map(
  umls_map = umls_map,
  efo_map = efo_map,
  do_map = do_map,
  oncotree_release = oncotree_release,
  efo_release = efo_release,
  do_release = disease_ontology_release)


oncotree_basic <- onco_map$oncotree_basic
oncotree_expanded_full <- onco_map$oncotree_expanded_full
oncotree_expanded_slim <- onco_map$oncotree_expanded_slim
auxiliary_maps <- list()
auxiliary_maps[['umls']] <- umls_map
auxiliary_maps[['do']] <- do_map
auxiliary_maps[['efo']] <- efo_map
auxiliary_maps[['icd10']] <- onco_map$icd10_map


version_date <- format(Sys.Date(), format="%Y%m%d")

usethis::use_data(oncotree_basic, overwrite = T)
usethis::use_data(oncotree_expanded_full, overwrite = T)
usethis::use_data(oncotree_expanded_slim, overwrite = T)
usethis::use_data(auxiliary_maps, overwrite = T)
