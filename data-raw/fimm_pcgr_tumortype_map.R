
# Load all the ICD10 codes
icd10_map <- map_icd10(
  basedir = here::here()
)

# Load phenOncoX icd10 to site mappings
icd10_to_site <- readRDS(
    "data-raw/gd_local/oncotree_expanded_v0.9.2.rds")$records |>
  dplyr::select(primary_site, icd10_code) |>
  dplyr::rename(identifier = icd10_code) |>
  dplyr::mutate(identifier = stringr::str_replace(
    identifier, "\\.[0-9]$",""
  )) |>
  dplyr::distinct() |>
  dplyr::filter(!is.na(primary_site)) |>
  dplyr::filter(!is.na(identifier))

# Load FIMM's ICD10 terms
fimm_icd10 <- readr::read_tsv(
  file="data-raw/fimm_icd10/fimm_icd10.tsv", 
  show_col_types = F) |>
  janitor::clean_names() |>
  
  dplyr::filter(
    !stringr::str_detect(
      .data$identifier, "^D(5|6|7|8)"
    )
  ) |>
  dplyr::mutate(identifier = stringr::str_replace(
    identifier, "&$", ""
  )) |>
  dplyr::left_join(
    icd10_to_site, by = "identifier",
    relationship = "many-to-many"
  ) |>
  dplyr::mutate(primary_site = dplyr::case_when(
    is.na(primary_site) & (
      identifier == "C00-C14" |
        identifier == "C30" |
        identifier == "C31" |
        identifier == "C32" |
        stringr::str_detect(parent_terms, "^C00-C14")) ~ "Head and Neck",
    is.na(primary_site) & (
        identifier == "C40" |
        identifier == "C41")  ~ "Bone",
    is.na(primary_site) & (
      identifier == "C17")  ~ "Esophagus/Stomach",
    is.na(primary_site) & (
      identifier == "C64")  ~ "Kidney",
    is.na(primary_site) & (
      identifier == "C70" |
        identifier == "D42" |
        identifier == "D43")  ~ "CNS/Brain",
    is.na(primary_site) & (
      identifier == "C49" |
        identifier == "D21")  ~ "Soft Tissue",
    is.na(primary_site) & (
      identifier == "C47")  ~ "Peripheral Nervous System",
    is.na(primary_site) & (
      identifier == "C63")  ~ "Testis",
    is.na(primary_site) & (
      identifier == "D19")  ~ "Pleura",
    is.na(primary_site) & (
      identifier == "C23")  ~ "Biliary Tract",
    is.na(primary_site) & (
      identifier == "C58" |
        identifier == "D26")  ~ "Uterus",
    is.na(primary_site) & (
      identifier == "D04" |
        identifier == "D22" |
        identifier == "D23")  ~ "Skin",
    is.na(primary_site) & (
      identifier == "D06")  ~ "Cervix",
    is.na(primary_site) & (
      identifier == "D34")  ~ "Thyroid",
    is.na(primary_site) & (
      identifier == "C39" |
        identifier == "D15" |
        identifier == "C33")  ~ "Lung",
    is.na(primary_site) & (
      identifier == "C55" |
      identifier == "C65" |
        identifier == "C68" |
        identifier == "D41")  ~ "Bladder/Urinary Tract",
    is.na(primary_site) & (
      identifier == "C19" |
        identifier == "D12" |
        identifier == "C20" |
        identifier == "C21")  ~ "Colon/Rectum",
    TRUE ~ as.character(primary_site)
  )) |>
  dplyr::filter(
    !stringr::str_detect(identifier, "-")
  )

fimm_icd10_final <- 
  fimm_icd10 |> 
  dplyr::group_by(reference_database, identifier) |> 
  dplyr::summarise(
    name = paste(unique(name), collapse="; "), 
    primary_site = paste(unique(primary_site), collapse=";"),
    .groups = "drop") |>
  ## Fix some of the primary sites (misassignments)
  dplyr::mutate(primary_site = dplyr::case_when(
    identifier == "C92" ~ "Myeloid",
    identifier == "C96" ~ "Lymphoid;Myeloid",
    identifier == "C78" ~ "Esophagus/Stomach;Lung",
    identifier == "C46" ~ "Skin;Soft Tissue",
    identifier == "C84" ~ "Lymphoid",
    identifier == "D20" ~ "Peritoneum;Soft Tissue",
    identifier == "C79" | identifier == "D35" ~ as.character(NA),
    TRUE ~ as.character(primary_site)
  )) |>
  dplyr::mutate(primary_site = dplyr::case_when(
    primary_site == "NA" & stringr::str_detect(name, "breast") ~ "Breast",
    primary_site == "NA" & stringr::str_detect(name, "ovary") ~ "Ovary/Fallopian Tube",
    primary_site == "NA" & stringr::str_detect(name, " eye") ~ "Eye",
    primary_site == "NA" & stringr::str_detect(name, "bone") ~ "Bone",
    primary_site == "NA" & stringr::str_detect(name, "digestive") ~ "Esophagus/Stomach",
    primary_site == "NA" & stringr::str_detect(name, "lipomatous") ~ "Soft Tissue",
    primary_site == "NA" & stringr::str_detect(name, "meninges") ~ "CNS/Brain",
    primary_site == "NA" & stringr::str_detect(name, "mouth|middle ear") ~ "Head and Neck",
    TRUE ~ as.character(primary_site)
  )) |>
  dplyr::mutate(primary_site = dplyr::case_when(
    stringr::str_detect(name,"uterus") & 
      !stringr::str_detect(primary_site,"Uterus") ~ "Uterus",
    identifier == "C76" | identifier == "D09" ~ as.character(NA),
    TRUE ~ as.character(primary_site)
  )) |>
  dplyr::rename(tumor_site_pcgr = primary_site,
                icd10_code = identifier) |>
  dplyr::distinct()

readr::write_tsv(
  fimm_icd10_final, "data-raw/fimm_icd10/fimm_icd10_pcgr_20250501.tsv")
