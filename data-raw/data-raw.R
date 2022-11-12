source('data-raw/pheno_oncox_utilities.R')

disease_ontology_release <- 'v2022-10'
efo_release <- 'v3.47.0'
oncotree_release <- '2021_11_02'

## get metadata from metadata.xlsx
metadata <- list()
for (elem in c("pheno_oncox")) {
  metadata[[elem]] <- as.data.frame(readxl::read_excel(
    "data-raw/metadata.xlsx",
    sheet = elem, col_names = TRUE
  ) |>
    dplyr::mutate(version = dplyr::if_else(
      is.na(version) &
        abbreviation == "umls",
      as.character(Sys.Date()),
      as.character(version)
    )))
}

## set logging layout
lgr::lgr$appenders$console$set_layout(
  lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T")
)
## Get UMLS / DiseaseOntology / EFO mappings
umls_map <- map_umls(
  update = T,
  basedir = here::here())

icd10_map <- map_icd10(
  basedir = here::here()
)

do_map <- map_disease_ontology(
  skip_non_cui_mapped = T,
  release = metadata$pheno_oncox[
    metadata$pheno_oncox$source == "Disease Ontology",]$version,
  umls_map = umls_map,
  basedir = here::here())


efo_map <- map_efo(
  umls_map = umls_map,
  efo_release = metadata$pheno_oncox[
    metadata$pheno_oncox$source == "Experimental Factor Ontology",]$version,
  update = T,
  basedir = here::here())

## Use OncoTree as starting point for phenotype cross-map
onco_map <- onco_pheno_map(
  umls_map = umls_map,
  efo_map = efo_map,
  do_map = do_map,
  icd10_map = icd10_map)


oncotree_expanded <- list()
oncotree_expanded$records <- onco_map$oncotree_expanded
oncotree_expanded$metadata <- metadata$pheno_oncox

oncotree_core <- list()
oncotree_core$records <- onco_map$oncotree
oncotree_core$metadata <- metadata$pheno_oncox[1,]

auxiliary_maps <- list()
auxiliary_maps$records <- list()

auxiliary_maps$records[['umls']] <- umls_map
auxiliary_maps$records[['do']] <- do_map
auxiliary_maps$records[['efo']] <- efo_map
auxiliary_maps$records[['icd10']] <- icd10_map
auxiliary_maps$metadata <- metadata$pheno_oncox[1:4,]

db <- list()
db[['oncotree_core']] <- oncotree_core
db[['oncotree_expanded']] <- oncotree_expanded
db[['auxiliary_maps']] <- auxiliary_maps

version_bumped <- "0.5.0"
gd_records <- list()
db_id_ref <- data.frame()


for (elem in c("oncotree_core", "auxiliary_maps", "oncotree_expanded")) {
  saveRDS(db[[elem]],
          file = paste0(
            "data-raw/gd_local/", elem, "_v",
            version_bumped, ".rds"
          )
  )
  
  (gd_records[[elem]] <- googledrive::drive_upload(
    paste0("data-raw/gd_local/", 
           elem, "_v", version_bumped, ".rds"),
    paste0("phenoOncoX/", elem, "_v", version_bumped, ".rds")
  ))
  
  google_rec_df <-
    dplyr::select(
      as.data.frame(gd_records[[elem]]), name, id
    ) |>
    dplyr::rename(
      gid = id,
      filename = name
    ) |>
    dplyr::mutate(
      name = stringr::str_replace(filename, "_v\\S+$", ""),
      date = as.character(Sys.Date()),
      pVersion = version_bumped
    ) |>
    dplyr::mutate(
      md5Checksum =
        gd_records[[elem]]$drive_resource[[1]]$md5Checksum
    )
  
  db_id_ref <- db_id_ref |>
    dplyr::bind_rows(google_rec_df)
}

usethis::use_data(db_id_ref, internal = TRUE, overwrite = TRUE)

# usethis::use_data(oncotree_basic, overwrite = T)
# usethis::use_data(oncotree_expanded_full, overwrite = T)
# usethis::use_data(oncotree_expanded_slim, overwrite = T)
# usethis::use_data(auxiliary_maps, overwrite = T)
