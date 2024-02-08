#library(magrittr)
library(jsonlite)
library(dplyr)
library(stringr)
library(ontologyIndex)

map_icd10 <- function(basedir = NULL) {
  
  umls_icd10_mapping <- read.table(
    gzfile(
      file.path(basedir, "data-raw", 
                "disgenet", "disease_mappings.tsv.gz")), 
    sep = "\t", header = T, stringsAsFactors = F, quote = "") |> 
    dplyr::rename(cui = diseaseId, 
                  cui_name = name, 
                  icd10_code = code) |> 
    dplyr::filter(vocabulary == "ICD10") |>
    dplyr::select(cui, icd10_code) |>
    dplyr::distinct()
  
  return(umls_icd10_mapping)
  
}
map_efo <- function(umls_map,
                    efo_release = "v3.47.0",
                    update = T,
                    basedir = NULL) {
  
  if (!dir.exists(basedir)) {
    lgr::lgr$info("Base directory does not exist")
  } else{
    if (!file.exists(
      file.path(basedir,"data-raw", "efo", paste0(
        "efo.", efo_release, ".obo"))) |
      update == T) {
      download.file(
        paste0("https://github.com/EBISPOT/efo/releases/download/",
               efo_release,"/efo.obo"),
        destfile =
          file.path(basedir, "data-raw", "efo", paste0(
            "efo.", efo_release, ".obo"))
      )
    }
  }
  con <- file(
    description =
      file.path(
        basedir,"data-raw","efo",
        paste0("efo.",
               efo_release,".obo")),
    open = "r")
  lines <- readLines(con)
  close(con)
  
  i <- 1
  efo_id <- NA
  nci_t <- NA
  msh_id <- NA
  icd10_id <- NA
  snomed_id <- NA
  cui_all <- c()
  cui_close_all <- c()
  cui_exact <- NA
  cui_exact_all <- c()
  msh_all <- c()
  nci_all <- c()
  icd10_all <- c()
  snomed_all <- c()
  efo_map <- NULL
  name <- NA
  cui <- NA
  cui_close <- NA
  obsolete <- FALSE
  ancestor <- NA
  ancestors <- c()
  while (i <= length(lines)) {
    line <- lines[i]
    #cat(line,'\n')
    ## only include major ontologies (phenotype-related)
    if (stringr::str_detect(line, "^id: (EFO|HP|DOID|MONDO|GO|Orphanet):[0-9]{1,}$")) {
      if (!is.na(efo_id) &
         !is.na(name) &
         !stringr::str_detect(name, "measurement|^CS") &
         obsolete == FALSE) {
        df <-
          data.frame("efo_id" = efo_id,
                     "nci_t" = paste(unique(nci_all),collapse = ","),
                     "msh" = paste(unique(msh_all),collapse = ","),
                     "cui" = paste(unique(cui_all),collapse = ","),
                     "icd10" = paste(unique(icd10_all), collapse = ","),
                     "snomed" = paste(unique(snomed_all), collapse = ","),
                     "cui_exact" = paste(unique(cui_exact_all), collapse = ","),
                     "efo_name" = name,
                     "cui_close" = paste(unique(cui_close_all), collapse = ","),
                     "ancestors" = paste(ancestors,collapse = ","),
                     stringsAsFactors = F)
        efo_map <- rbind(efo_map, df)
      }
      do_id <- NA
      name <- NA
      nci_t <- NA
      msh_id <- NA
      cui <- NA
      icd10_id <- NA
      cui_close <- NA
      cui_exact <- NA
      snomed_id <- NA
      icd10_all <- c()
      cui_all <- c()
      cui_close_all <- c()
      cui_exact_all <- c()
      snomed_all <- c()
      msh_all <- c()
      nci_all <- c()
      cell_line <- FALSE
      obsolete <- FALSE
      ancestor <- NA
      ancestors <- c()
      efo_id <- stringr::str_replace(line,"^id: ","")
    }
    if (stringr::str_detect(line,"name: .+$")) {
      name <- stringr::str_replace(line, "name: ","")
      if ("http" %in% name) {
        name <- NA
      }
    }
    if (stringr::str_detect(line,"is_obsolete: true$")) {
      obsolete <- TRUE
    }
    if (stringr::str_detect(line,"is_a: EFO:[0-9]{1,} ")) {
      ancestor <- stringr::str_replace(
        stringr::str_match(line, "EFO:[0-9]{1,} ")[[1]]," ","")
      ancestors <- c(ancestors,ancestor)
    }
    
    if (stringr::str_detect(line,"xref: MSH:[A-Z]{1,2}[0-9]{1,}")) {
      msh_id <- stringr::str_replace_all(line,"xref: MSH:", "")
      msh_all <- c(msh_all, msh_id)
    }
    if (stringr::str_detect(line,"xref: SNOMEDCT:[0-9]{1,}")) {
      snomed_id <- stringr::str_replace_all(line,"xref: SNOMEDCT:", "")
      snomed_all <- c(snomed_all, snomed_id)
    }
    if (stringr::str_detect(line,"xref: ICD10:[A-Z]{1}[0-9]{1,2}\\.[0-9]{1,2}")) {
      icd10_id <- stringr::str_trim(stringr::str_replace_all(line, 
                                                             "xref: ICD10:",""))
      icd10_id <- stringr::str_replace_all(icd10_id," \\{.+\\}$","")
      icd10_all <- c(icd10_all, icd10_id)
    }
    if (stringr::str_detect(line,"xref: UMLS:[A-Z]{1,2}[0-9]{1,}")) {
      cui <- stringr::str_replace_all(
        stringr::str_match(line,"xref: UMLS:[A-Z]{1,2}[0-9]{1,}")[[1]],
        "xref: UMLS:","")
      cui_all <- c(cui_all, cui)
    }
    if (stringr::str_detect(
      line,
      "property_value: \"closeMatch\" http://linkedlifedata.com/resource/umls/id/")) {
      cui_close <- stringr::str_replace_all(
        line,
        "property_value: \"closeMatch\" http://linkedlifedata.com/resource/umls/id/","")
      cui_close_all <- c(cui_close_all, cui_close)
    }
    if (stringr::str_detect(
      line,
      "property_value: exactMatch http://linkedlifedata.com/resource/umls/id/")) {
      cui_exact <- stringr::str_replace_all(
        line,
        "property_value: exactMatch http://linkedlifedata.com/resource/umls/id/","")
      cui_exact_all <- c(cui_exact_all, cui_exact)
    }
    if (stringr::str_detect(
      line,
      "^xref: NCIt:")) {
      nci_t <- stringr::str_replace_all(
        line,
        "xref: NCIt:","")
      nci_all <- c(nci_all, nci_t)
    }
    if(i %% 10000 == 0){
      cat(paste0("Processing line: ", i), sep="",'\n')
    }
    i <- i + 1
  }
  
  efo_map <- as.data.frame(
    efo_map |>
      dplyr::mutate(cui_close =
                      dplyr::if_else(cui_close == "",
                                     as.character(NA),cui_close)) |>
      dplyr::mutate(nci_t = dplyr::if_else(nci_t == "",
                                           as.character(NA),nci_t)) |>
      dplyr::mutate(msh = dplyr::if_else(msh == "",
                                         as.character(NA),msh)) |>
      dplyr::mutate(cui = dplyr::if_else(cui == "",
                                         as.character(NA),cui)) |>
      tidyr::separate_rows(cui, sep = ",") |>
      tidyr::separate_rows(cui_close, sep = ",") |>
      tidyr::separate_rows(nci_t, sep = ",") |>
      tidyr::separate_rows(msh, sep = ",") |>
      dplyr::mutate(
        efo_name =
          stringr::str_replace(
            efo_name,
            " \\{http://www.co-ode.org/patterns#createdBy=","")) |>
      dplyr::mutate(
        efo_name =
          stringr::str_replace(
            efo_name,'"http://www.ebi.ac.uk/ontology/webulous#OPPL_pattern"\\}',
            '')) |>
      dplyr::distinct() 
  )
  
  efo2name <- efo_map |>
    dplyr::select(efo_id, efo_name) |>
    dplyr::distinct() |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(
          tolower(efo_name), "carcinoma|cancer|neoplasm|tumor") &
          stringr::str_detect(
            tolower(efo_name), 
            "colon|anal|anus|cecum|rectum|intestinal|small intestine|rectal|colorectal"),
        "Colon/Rectum",
        as.character(NA)
      )
    )  |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(
          tolower(efo_name), 
          "bladder|urinary|ureter|urethra|urothelial|transitional cell ") &
          stringr::str_detect(
            tolower(efo_name), "carcinoma|cancer|neoplasm|tumor"),
        "Bladder/Urinary Tract",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "ampulla of vater"),
        "Ampulla of Vater",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(
          tolower(efo_name), "pancreas|acinar cell|pancreatic") &
          stringr::str_detect(
            tolower(efo_name), "carcinoma|cancer|neoplasm|tumor"),
        "Pancreas",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "thymoma") |
        (stringr::str_detect(tolower(efo_name), "thymic|thymus") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplasm|tumor")),
        "Thymus",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name),"schwannoma") |
          (stringr::str_detect(tolower(efo_name), 
                               "neuroectodermal|peripheral nervous|peripheral nerve sheath") &
             stringr::str_detect(tolower(efo_name), 
                                 "carcinoma|cancer|neoplasm|tumor")),
        "Peripheral Nervous System",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "thyroid") &
          stringr::str_detect(tolower(efo_name), 
                              "carcinoma|cancer|neoplasm|tumor"),
        "Thyroid",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "^myelo") |
        stringr::str_detect(tolower(efo_name), "leukemia|myeloma") |
          (stringr::str_detect(tolower(efo_name), "t-cell|b-cell") &
             stringr::str_detect(tolower(efo_name), "cancer|tumor|neoplasm")) |
          (stringr::str_detect(tolower(efo_name), "myeloid") &
             stringr::str_detect(tolower(efo_name), "leukemia|neoplasm")),
        "Myeloid",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), 
                            "gastric|stomach|esophagus|esophageal") &
          stringr::str_detect(tolower(efo_name), 
                              "carcinoma|cancer|neoplas|tumor"),
        "Esophagus/Stomach",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "retinoblastoma") |
        (stringr::str_detect(
          tolower(efo_name), 
          "eye |ocular |orbit(al)? |retinal |((uveal|ocular) melanoma)|lacrimal gland|palbrepal ") &
           stringr::str_detect(tolower(efo_name), 
                               "carcinoma|cancer|neoplas|tumor")),
        "Eye",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), 
                            "gallbladder|biliary tract|bile duct|hepatobiliary|cholangio|biliary intraepithelial") &
          stringr::str_detect(tolower(efo_name), 
                              "carcinoma|cancer|neoplas|tumor"),
        "Biliary Tract",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "chondroblastom|chordoma|giant cell") |
          (stringr::str_detect(tolower(efo_name), "bone") &
             stringr::str_detect(tolower(efo_name), "carcinoma|sarcom|cancer|neoplas|tumor")),
        "Bone",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        !stringr::str_detect(tolower(efo_name), "bone") &
          stringr::str_detect(tolower(efo_name), 
                              "sarcom|soft tissue|gastrointestinal stromal|leiomy|connective tissue"),
        "Soft Tissue",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "prostat") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor"),
        "Prostate",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "ovary|ovarian|fallopian tube") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|blastom|tumor"),
        "Ovary/Fallopian Tube",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "uterine|endometri|female reproductive") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor"),
        "Uterus",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "cervix|cervical") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor"),
        "Cervix",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "tonsil |tonsillar |head and neck|mouth|neck|glottis|larynx|pharynx |pharyngeal |gum |lip |parotid gland|salivary gland|oral squamous|tongue|nasal cavity|nasopharyngeal|laryngeal|sinus |hypopharyn|oral cavity|oropharynx") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor"),
        "Head and Neck",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "liver|hepatocellular") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor"),
        "Liver",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "peritone") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor|mesotheliom"),
        "Peritoneum",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "penis|penile") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor"),
        "Penis",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name),"mesotheliom|mesothelial") |
          (stringr::str_detect(tolower(efo_name), "pleura") &
             stringr::str_detect(tolower(efo_name), "carcinoma|cancer|mesotheliom|neoplas|tumor")),
        "Pleura",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "thoracic|lung|bronchi|bronchogenic|bronchus|bronchoalveolar|respiratory system|large cell neuroendocrine") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor"),
        "Lung",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "neuroblastom|meningiom|glioblastom|medulloblastom|glioma|cerebellum cancer|hemangioblastom|cerebellar neoplasm|astrocytom|scwhannom") |
          (stringr::str_detect(tolower(efo_name), "glioneuronal|central nervous|skull |pituitary |nervous system|neuronal |brain|cerebellar|cerebral") &
             stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor|teratoma")),
        "CNS/Brain",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "breast|nipple") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor"),
        "Breast",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "kidney|renal|wilms|clear cell|nephroblastom") &
          !stringr::str_detect(tolower(efo_name),"ovarian|cervical") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor"),
        "Kidney",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name),"seminoma") |
          (stringr::str_detect(tolower(efo_name), "testis|testicular|embryonal|male reproductive|germ cell") &
             !stringr::str_detect(tolower(efo_name),"female|endometri|ovarian") &
             stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor")),
        "Testis",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name),"melanoma|ear cancer") |
          (stringr::str_detect(tolower(efo_name), "skin|keratinocyte|eccrine") &
             stringr::str_detect(tolower(efo_name), "squamous cell|carcinoma|cancer|neoplas|tumor")),
        "Skin",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        (stringr::str_detect(tolower(efo_name), "vulva|vagina|bartholin gland") &
           stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor")),
        "Vulva/Vagina",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "lynch|li-fraumeni|cowden|von hippel|myelodysplastic syndrome|familial adenomatous ") |
          (stringr::str_detect(tolower(efo_name), "syndrome") &
             stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor")) |
          (stringr::str_detect(tolower(efo_name), "hereditary|familial|susceptibility") &
             stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neuroblast|sarcoma|glioma|melanoma|myeloma|lymphoma|neoplas|tumor")),
        "Other/Unknown",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        #stringr::str_detect(tolower(cui_name), "waldenstr") |
        stringr::str_detect(tolower(efo_name), "lymphom") |
          (stringr::str_detect(tolower(efo_name), "cancer|neoplasm") &
             stringr::str_detect(tolower(efo_name), "lympho|hematopoietic"))
        ,
        "Lymphoid",
        as.character(primary_site)
      )
    ) |>
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "adrenal gland|adrenal") &
          stringr::str_detect(tolower(efo_name), "carcinoma|pheochromocytoma|cancer|neoplas|tumor|blastom|lipom"),
        "Adrenal Gland",
        as.character(primary_site)
      )
    )
  
  cui_map <- efo_map |>
    dplyr::filter(!is.na(cui)) |>
    dplyr::select(efo_id, cui) |>
    dplyr::left_join(
      dplyr::filter(umls_map$concept, main_term == T),
      by = "cui", multiple = "all", relationship = "many-to-many") |>
    dplyr::mutate(xref_source = "UMLS") |>
    dplyr::select(efo_id,cui,cui_name, xref_source) |>
    dplyr::distinct()
  
  nci_map <- efo_map |>
    dplyr::filter(!is.na(nci_t)) |>
    dplyr::select(efo_id, nci_t) |>
    dplyr::left_join(
      umls_map$nciXref, 
      by = "nci_t", 
      multiple = "all",
      relationship = "many-to-many") |>
    dplyr::left_join(
      dplyr::filter(
        umls_map$concept,
        main_term == T), 
      by = "cui", multiple = "all",
      relationship = "many-to-many") |>
    dplyr::filter(!is.na(cui)) |>
    dplyr::mutate(xref_source = "NCI") |>
    dplyr::select(efo_id, cui, 
                  cui_name, xref_source) |>
    dplyr::distinct()
  
  msh_map <- efo_map |>
    dplyr::filter(!is.na(msh)) |>
    dplyr::select(efo_id,msh) |>
    dplyr::left_join(
      umls_map$mshXref, 
      by = "msh", multiple = "all",
      relationship = "many-to-many") |>
    dplyr::left_join(
      dplyr::filter(
        umls_map$concept, main_term == T),
      by = "cui", multiple = "all",
      relationship = "many-to-many") |>
    dplyr::filter(!is.na(cui)) |>
    dplyr::mutate(xref_source = "MESH") |>
    dplyr::select(efo_id, cui, 
                  cui_name, xref_source) |>
    # dplyr::mutate(cui_name =
    #                 stringi::stri_enc_toascii(cui_name)) |>
    dplyr::distinct()
  
  snomed_map <- efo_map |>
    dplyr::filter(!is.na(snomed)) |>
    dplyr::select(efo_id, snomed) |>
    tidyr::separate_rows(snomed) |>
    dplyr::left_join(
      umls_map$snomedXref, 
      by = "snomed", multiple = "all",
      relationship = "many-to-many") |>
    dplyr::left_join(
      dplyr::filter(
        umls_map$concept, main_term == T),
      by = "cui", multiple = "all",
      relationship = "many-to-many") |>
    dplyr::filter(!is.na(cui)) |>
    dplyr::mutate(xref_source = "SNOMED") |>
    # dplyr::mutate(cui_name =
    #                 stringi::stri_enc_toascii(cui_name)) |>
    dplyr::select(efo_id, cui, 
                  cui_name, xref_source) |>
    dplyr::distinct()
  
  efo2xref <- dplyr::bind_rows(
    cui_map, 
    #cui_map_close, 
    msh_map, 
    nci_map, 
    snomed_map) |>
    dplyr::distinct() |>
    dplyr::left_join(
      efo2name, by = "efo_id",
      multiple = "all",
      relationship = "many-to-many") |>
    dplyr::filter(
      is.na(cui) | 
        (!is.na(cui) & !is.na(cui_name)))
  
  return(list("efo2name" = efo2name,
              "efo2xref" = efo2xref))
}


map_disease_ontology <- function(
  skip_non_cui_mapped = T,
  release = "2022-10",
  umls_map = NULL,
  basedir = NULL) {
  
  do_github_raw_url <-
    paste0("https://raw.githubusercontent.com/DiseaseOntology/",
           "HumanDiseaseOntology/main/src/ontology/releases")
  
  release_dest <- stringr::str_replace_all(release,"-","")
  if (!file.exists(
    file.path(basedir,
              "data-raw",
              "do",
              paste0("TopNodes_Docancerslim.", release_dest,".obo")))) {
    download.file(
      url =
        paste0(do_github_raw_url,
               "/subsets/TopNodes_DOcancerslim.obo"),
      destfile = file.path(
        basedir, "data-raw","do",
        paste0("TopNodes_Docancerslim.", release_dest,".obo"))
    )
  }
  do_cancer_top <-
    ontologyIndex::get_ontology(
      file.path(basedir, "data-raw","do",
                paste0(
                  "TopNodes_DOcancerslim.", release_dest, ".obo")
      ),
      extract_tags = "everything")
  do_cancer_top_id <- do_cancer_top$id
  names(do_cancer_top_id) <- NULL
  do_cancer_top_name <- do_cancer_top$name
  names(do_cancer_top_name) <- NULL
  do_cancer_top_df <-
    data.frame("do_id" = do_cancer_top_id,
               "do_name" = do_cancer_top_name, stringsAsFactors = F)
  do_cancer_top_df <-
    dplyr::filter(do_cancer_top_df, stringr::str_detect(do_id,"DOID:"))
  if (!file.exists(
    file.path(basedir, "data-raw","do",
              paste0("doid.",release_dest,".obo")))) {
    download.file(
      url =
        paste0(do_github_raw_url, "/doid.obo"),
      destfile =
        file.path(basedir, "data-raw", "do",
                  paste0("doid.", release_dest, ".obo")))
  }
  do_index <-
    ontologyIndex::get_ontology(
      file.path(basedir,"data-raw","do",
                paste0("doid.", release_dest, ".obo")),
      extract_tags = "everything")
  i <- 1
  
  umls_map_lower <- umls_map$concept |>
    dplyr::mutate(cui_name_lc = tolower(cui_name))
  
  do_ids <- do_index$id
  do_names <- do_index$name
  names(do_ids) <- NULL
  names(do_names) <- NULL
  do_map <- data.frame()
  do_map_top <- data.frame()
  while(i < length(do_ids)) {
    do_id <- do_ids[i]
    if (stringr::str_detect(do_id,"DOID:")) {
      do_name <- do_names[i]
      all_ancestors <-
        ontologyIndex::get_term_property(
          ontology = do_index,
          property_name ="ancestors",
          term = do_id)
      all_xref <-
        ontologyIndex::get_term_property(
          ontology = do_index,
          property_name = "xref",
          term = do_id)
      all_subset <-
        ontologyIndex::get_term_property(
          ontology = do_index,
          property_name = "subset",
          term = do_id)
      cui <- NA
      do_cancer_slim <- FALSE
      do_rare_slim <- FALSE
      do_cancer_slim_top <- FALSE
      top_node_cancer_do_id <- NA
      top_node_cancer <- data.frame()
      cui_ids <- c()
      if (length(all_xref) > 0) {
        j <- 1
        while(j <= length(all_xref)) {
          if (stringr::str_detect(all_xref[j],"^UMLS_CUI:C[0-9]{1,}$")) {
            cui <- stringr::str_replace(all_xref[j],"^UMLS_CUI:","")
            cui_ids <- c(cui_ids,cui)
          }
          j <- j + 1
        }
      }
      if (length(all_ancestors) > 0) {
        all_ancestor_df <-
          data.frame("do_id" = all_ancestors,
                     stringsAsFactors = F)
        top_node_cancer <-
          as.data.frame(
            dplyr::inner_join(
              all_ancestor_df,
              do_cancer_top_df,by="do_id") |>
              dplyr::distinct())
      }
      if (length(all_subset) > 0) {
        if ("TopNodes_DOcancerslim" %in% all_subset) {
          do_cancer_slim_top <- T
        }
        if ("DO_cancer_slim" %in% all_subset) {
          do_cancer_slim <- T
        }
      }
      
      
      if (length(cui_ids) == 0) {
        do_entry <-
          data.frame("do_id" = do_id,
                     "do_name" = do_name,
                     "cui" = NA,
                     "do_cancer_slim" = do_cancer_slim,
                     stringsAsFactors = F)
        do_map <- rbind(do_map, do_entry)
      }
      else{
        u <- 1
        while(u <= length(cui_ids)) {
          do_entry <-
            data.frame("do_id" = do_id,
                       "do_name" = do_name,
                       "cui" = cui_ids[u],
                       "do_cancer_slim" = do_cancer_slim,
                       stringsAsFactors = F)
          do_map <- rbind(do_map, do_entry)
          u <- u + 1
        }
      }
      
      if (nrow(top_node_cancer) > 0) {
        k <- 1
        while(k <= nrow(top_node_cancer)) {
          if (length(cui_ids) == 0) {
            do_entry <-
              data.frame("do_id" = do_id, "do_name" = do_name,
                         "cui" = NA,
                         "do_cancer_slim_top" = do_cancer_slim_top,
                         "top_node_cancer_slim" = top_node_cancer[k,]$do_id,
                         stringsAsFactors = F)
            do_map_top <- rbind(do_map_top, do_entry)
          }
          else{
            u <- 1
            while(u <= length(cui_ids)) {
              do_entry <-
                data.frame("do_id" = do_id,
                           "do_name" = do_name, "cui" = cui_ids[u],
                           "do_cancer_slim_top" = do_cancer_slim_top,
                           "top_node_cancer_slim" = top_node_cancer[k,]$do_id,
                           stringsAsFactors = F)
              do_map_top <- rbind(do_map_top, do_entry)
              u <- u + 1
            }
          }
          k <- k + 1
        }
      }
      else{
        if (length(cui_ids) == 0) {
          do_entry <-
            data.frame("do_id" = do_id,
                       "do_name" = do_name,
                       "cui" = NA,
                       "do_cancer_slim_top" = do_cancer_slim_top,
                       "top_node_cancer_slim" = top_node_cancer[k,]$do_id,
                       stringsAsFactors = F)
          do_map_top <- rbind(do_map_top, do_entry)
        }
        else{
          u <- 1
          while(u <= length(cui_ids)) {
            do_entry <-
              data.frame("do_id" = do_id,
                         "do_name" = do_name,
                         "cui" = cui_ids[u],
                         "do_cancer_slim_top" = do_cancer_slim_top,
                         "top_node_cancer_slim" = NA, stringsAsFactors = F)
            do_map_top <- rbind(do_map_top, do_entry)
            u <- u + 1
          }
        }
      }
    }
    i <- i + 1
  }
  
  ## manual correction of erroneous or missing UMLS cross-references
  do_map <- do_map |>
    # dplyr::mutate(do_name =
    #                 stringi::stri_enc_toascii(do_name)) |>
    ## sarcoma (non-existent UMLS_CUI = C0153519)
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:1115",
                                       "C1261473",
                                       as.character(cui))) |>
    ## colorectal cancer (UMLS_CUI = C1261473, not part of tree)
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:9256",
                                       "C0009404", as.character(cui))) |>
    ## hepatocellular carcinoma (missing UMLS_CUI)
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:684",
                                       "C2239176", as.character(cui))) |>
    ## chronic leukemia (missing UMLS_CUI)
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:1036",
                                       "C1279296", as.character(cui))) |>
    ## lymphoid leukemia (missing UMLS_CUI)
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:10747",
                                       "C0023448", as.character(cui))) |>
    ## brain glioma (missing UMLS_CUI)
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:0060108",
                                       "C0349661", as.character(cui))) |>
    ## ovarian serious carcinoma (missing UMLS_CUI)
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:0050933",
                                       "C0279663", as.character(cui))) |>
    ## mucosal melanoma (missing UMLS_CUI)
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:0050929",
                                       "C3898222", as.character(cui))) |>
    ## acral lentiginous melanoma (missing UMLS_CUI)
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:6367",
                                       "C0346037", as.character(cui))) |>
    ## brain glioma (missing UMLS_CUI)
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:0080146",
                                       "C0279584",as.character(cui))) |>
    ## correct urothelial carcinoma (transitional cell carcinoma)
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:2671",
                                       "C2145472", as.character(cui))) |>
    ## colon cancer
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:219",
                                       "C0699790", as.character(cui))) |>
    ## stomach cancer
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:10534",
                                       "C0699791", as.character(cui))) |>
    ## skip myelodysplastic syndrome
    dplyr::filter(do_id != "DOID:0050908") |>
    ## B-lymphoblastic Leukemia/lymphoma, BCR-ABL1–like
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:0080650",
                                       "C4329382", as.character(cui))) |>
    ## Childhood B-cell Acute Lymphoblastic Leukemia
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:0080146",
                                       "C0279584", as.character(cui))) |>
    ## Childhood Low-grade Glioma
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:0080830",
                                       "C1997217", as.character(cui))) |>
    ## Estrogen-receptor positive breast cancer
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:0060075",
                                       "C4745240", as.character(cui))) |>
    ## HER2-receptor negative breast cancer
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:0060080",
                                       "C5238910", as.character(cui))) |>
    ## Diffuse Midline Glioma, H3 K27M-mutant
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:0080684",
                                       "C4289688", as.character(cui))) |>
    ## Diffuse large B-cell lymphoma activated B-cell type
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:0080996",
                                       "C1333296", as.character(cui))) |>
    ## Oligodendroglioma
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:3181",
                                       "C0028945", as.character(cui))) |>
    ## B-cell Acute Lymphoblastic Leukemia
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:0080638",
                                       "C1292769", as.character(cui))) |>
    ## B-cell Adult Acute Lymphocytic Leukemia
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:0060592",
                                       "C0279593", as.character(cui))) |>
    ## Diffuse Glioma, H3 G34 Mutant
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:0080880",
                                       "C5669880", as.character(cui))) |>
    
    dplyr::bind_rows(
      data.frame(do_id = "DOID:1788",
                 do_name = "peritoneal mesothelioma",
                 cui = "C1377610", do_cancer_slim = T,
                 stringsAsFactors = F)) |>
    dplyr::bind_rows(data.frame(do_id = "DOID:1790",
                                do_name = "mesothelioma",
                                cui = "C0025500", do_cancer_slim = T,
                                stringsAsFactors = F)) |>
    dplyr::bind_rows(
      data.frame(do_id = "DOID:1612", do_name = "breast cancer",
                 cui = "C0678222",  do_cancer_slim = T,
                 stringsAsFactors = F)) |>
    dplyr::bind_rows(
      data.frame(do_id = "DOID:12603", do_name = "acute leukemia",
                 cui = "C0085669",  do_cancer_slim = T,
                 stringsAsFactors = F)) |>
    dplyr::bind_rows(
      data.frame(do_id = "DOID:1793", do_name = "pancreatic cancer",
                 cui = "C0030297",  do_cancer_slim = T,
                 stringsAsFactors = F)) |>
    dplyr::bind_rows(
      data.frame(do_id = "DOID:1793", do_name = "pancreatic cancer",
                 cui = "C0887833",  do_cancer_slim = T,
                 stringsAsFactors = F)) |>
    dplyr::bind_rows(
      data.frame(do_id = "DOID:1793", do_name = "pancreatic cancer",
                 cui = "C0235974",  do_cancer_slim = T,
                 stringsAsFactors = F)) |>
    dplyr::bind_rows(
      data.frame(do_id = "DOID:8923", do_name = "skin melanoma",
                 cui = "C0025202",  do_cancer_slim = T,
                 stringsAsFactors = F)) |>
    dplyr::bind_rows(
      data.frame(do_id = "DOID:1324", do_name = "lung cancer",
                 cui = "C0684249",  do_cancer_slim = T,
                 stringsAsFactors = F)) |>
    dplyr::bind_rows(
      data.frame(do_id = "DOID:0050908", do_name = "myelodysplastic syndrome",
                 cui = "C3463824",  do_cancer_slim = T,
                 stringsAsFactors = F)) |>
    dplyr::bind_rows(
      data.frame(do_id = "DOID:9256", do_name = "colorectal cancer",
                 cui = "C0699790", do_cancer_slim = T,
                 stringsAsFactors = F)) |>
    dplyr::bind_rows(
      data.frame(do_id = "DOID:0060081",
                 do_name = "triple-receptor negative breast cancer",
                 cui = "C3539878", do_cancer_slim = T,
                 stringsAsFactors = F)) |>
    dplyr::bind_rows(
      data.frame(do_id = "DOID:0060079",
                 do_name = "her2-receptor positive breast cancer",
                 cui = "C1960398", do_cancer_slim = T,
                 stringsAsFactors = F)) |>
    dplyr::bind_rows(
      data.frame(do_id = "DOID:0050745",
                 do_name = "diffuse large b-cell lymphoma",
                 cui = "C0079744", do_cancer_slim = T,
                 stringsAsFactors = F)) |>
    dplyr::bind_rows(
      data.frame(do_id = "DOID:0050861",
                 do_name = "colorectal adenocarcinoma",
                 cui = "C0338106", do_cancer_slim = T,
                 stringsAsFactors = F)) |>
    dplyr::bind_rows(
      data.frame(do_id = "DOID:0060105",
                 do_name = "brain medulloblastoma",
                 cui = "C1332188", do_cancer_slim = T,
                 stringsAsFactors = F)) |>
    dplyr::bind_rows(
      data.frame(do_id = "DOID:0050744",
                 do_name = "anaplastic large cell lymphoma",
                 cui = "C0206180", do_cancer_slim = T,
                 stringsAsFactors = F)) |>
    dplyr::bind_rows(
      data.frame(do_id = "DOID:0050746",
                 do_name = "mantle cell lymphoma",
                 cui = "C0334634", do_cancer_slim = T,
                 stringsAsFactors = F)) |>
    dplyr::bind_rows(
      data.frame(do_id = "DOID:0050749",
                 do_name = "peripheral t-cell lymphoma",
                 cui = "C0079774", do_cancer_slim = T,
                 stringsAsFactors = F)) |>
    dplyr::bind_rows(
      data.frame(do_id = "DOID:0050873",
                 do_name = "follicular lymphoma",
                 cui = "C0024301", do_cancer_slim = T,
                 stringsAsFactors = F)) |>
    dplyr::bind_rows(
      data.frame(do_id = "DOID:0070004", do_name = "myeloma",
                 cui = "C0026764", do_cancer_slim = T,
                 stringsAsFactors = F)) |>
    dplyr::bind_rows(
      data.frame(do_id = "DOID:3965",
                 do_name = "merkel cell carcinoma",
                 cui = "C0007129", do_cancer_slim = T,
                 stringsAsFactors = F))
  
  
  do_map_cui_name_matched <- do_map |>
    dplyr::select(do_id, do_name) |>
    dplyr::left_join(
      umls_map_lower, 
      by = c("do_name" = "cui_name_lc"),
      multiple = "all", relationship = "many-to-many") |>
    dplyr::filter(!is.na(cui)) |>
    dplyr::select(do_id, cui) |>
    dplyr::distinct()
  
  #do_map <- dplyr::left_join(do_map, do_map_cui_name_matched)
  
  if (skip_non_cui_mapped == T) {
    do_map1 <- do_map |> 
      dplyr::filter(!is.na(cui)) |> 
      dplyr::distinct()
    do_map2 <- do_map |> 
      dplyr::filter(is.na(cui)) |>
      dplyr::select(-c(cui)) |>
      dplyr::inner_join(
        do_map_cui_name_matched, 
        by = c("do_id"), 
        multiple = "all", relationship = "many-to-many")
    
    do_map <- dplyr::bind_rows(do_map1, do_map2)
  }
  
  return(do_map)
}

map_umls <- function(
  update = T,
  basedir = NULL) {
  
  if (is.null(basedir)) {
    lgr::lgr$info("Please specifiy valid base directory")
  }
  if (!dir.exists(basedir)) {
    lgr::lgr$info("Base directory does not exist")
  }
  
  for (fn in c('MGCONSO','NAMES','MGREL_1','MGREL_2')) {
    if (!file.exists(
      file.path(basedir, "data-raw", "umls", paste0(fn,".csv.gz"))) | update == T) {
      download.file(
        paste0("ftp://ftp.ncbi.nlm.nih.gov/pub/medgen/csv/",fn,".csv.gz"),
        method = "curl",
        destfile = file.path(
          basedir, "data-raw", "umls", paste0(fn,".csv.gz"))
      )
    }
  }
  
  
  if (!file.exists(file.path(
    basedir, "data-raw", "umls", "Neoplasm_Core.txt"))) {
    download.file(
      "https://evs.nci.nih.gov/ftp1/NCI_Thesaurus/Neoplasm/Neoplasm_Core.txt",
      destfile = file.path(basedir, "data-raw", "umls", "Neoplasm_Core.txt"))
  }
  
  
  neoplasm_core <-
    read.table(
      file.path(
        basedir, "data-raw", "umls", "Neoplasm_Core.txt"),
      sep = "\t",
      header = T, quote = "",
      na.strings = c(""), stringsAsFactors = F) |>
    janitor::clean_names()
  
  umls_rel_1 <- read.csv(
    gzfile(file.path(
      basedir, "data-raw","umls","MGREL_1.csv.gz")
    ),
    stringsAsFactors = F)
  umls_rel_2 <- read.csv(
    gzfile(file.path(
      basedir, "data-raw" ,"umls", "MGREL_2.csv.gz")
    ),
    stringsAsFactors = F)
  umls_rel <- rbind(umls_rel_1, umls_rel_2) |>
    dplyr::filter(REL == "CHD" & SUPPRESS == "N") |>
    dplyr::select(CUI1,CUI2) |>
    dplyr::distinct()
  umls_rel <- dplyr::rename(umls_rel, cui = CUI1, cui_child = CUI2)
  
  umls <- read.csv(
    gzfile(
      file.path(basedir, "data-raw","umls","MGCONSO.csv.gz")
    ),
    stringsAsFactors = F)
  
  concept_names_main <-
    read.csv(
      gzfile(
        file.path(basedir, "data-raw","umls","NAMES.csv.gz")
      ),
      stringsAsFactors = F) |>
    dplyr::select(CUI, name) |>
    dplyr::rename(STR = name) |>
    # dplyr::mutate(STR =
    #                 stringi::stri_enc_toascii(STR))
    dplyr::mutate(main_term = TRUE)
  
  concept_summary_data <-
    read.csv(
      gzfile(
        file.path(basedir,"data-raw","umls","MGCONSO.csv.gz")),
      stringsAsFactors = F) |>
    dplyr::select(CUI,SAB,STR) |>
    dplyr::distinct() |>
    dplyr::left_join(concept_names_main, by = c("CUI","STR")) |>
    dplyr::rename(cui = CUI, source = SAB, cui_name = STR) |>
    # dplyr::mutate(cui_name =
    #                 stringi::stri_enc_toascii(cui_name))
    dplyr::filter(nchar(cui) > 2) |>
    dplyr::mutate(main_term =
                    dplyr::if_else(is.na(main_term),
                                   FALSE,TRUE,TRUE))
  
  umls_nci <- dplyr::filter(umls, SAB == "NCI") |>
    dplyr::select(CUI,SCUI) |>
    dplyr::distinct() |>
    dplyr::rename(cui = CUI, nci_t = SCUI)
  
  umls_msh <- dplyr::filter(umls, SAB == "MSH") |>
    dplyr::select(CUI,SDUI) |>
    dplyr::distinct() |>
    dplyr::rename(cui = CUI, msh = SDUI)
  
  umls_snomed <- dplyr::filter(umls, SAB == "SNOMEDCT_US") |>
    dplyr::select(CUI, SCUI) |>
    dplyr::distinct() |>
    dplyr::rename(cui = CUI, snomed = SCUI)
  
  umls <- list("concept" = concept_summary_data,
               "relation" = umls_rel,
               "snomedXref" = umls_snomed,
               "nciXref" = umls_nci,
               "mshXref" = umls_msh)
  return(umls)
}


get_oncotree_entry_df <- function(
  tree_entry) {
  
  main_type <- tree_entry$mainType
  tissue_name <- tree_entry$tissue
  name <- tree_entry$name
  code <- tree_entry$code
  level <- tree_entry$level
  cui <- NA
  nci_t <- NA
  if ("UMLS" %in% names(tree_entry$externalReferences)) {
    cui <- tree_entry$externalReferences$UMLS
  }
  if ("NCI" %in% names(tree_entry$externalReferences)) {
    nci_t <- tree_entry$externalReferences$NCI
  }
  
  df <- data.frame(
    'tissue' = tissue_name,
    'main_type' = main_type,
    'name' = name,
    'code' = code,
    'level' = level,
    'cui' = cui,
    'nci_t' = nci_t,
    stringsAsFactors = F)
  
  return(df)
  
}

onco_pheno_map <- function(
  umls_map = NULL, 
  efo_map = NULL,
  do_map = NULL,
  icd10_map = NULL,
  oncotree_release = "2021_11_02",
  efo_release = NA,
  do_release = NA) {
  
  main_types_minor <-
    read.table(file = "data-raw/oncotree/oncotree.main_types_minor.txt",
               header = F,stringsAsFactors = F,
               quote = "", sep = "\t") |>
    dplyr::rename(main_type = V1) |>
    dplyr::mutate(minor_type = T)
  
  nci_map <- umls_map$nciXref |>
    dplyr::rename(cui_2 = cui)
  
  cui_name_map <- umls_map$concept |>
    dplyr::filter(main_term == T) |>
    dplyr::select(cui, cui_name) |>
    dplyr::distinct()
  
  ot_tumor_types_tree <-
    jsonlite::fromJSON("http://oncotree.mskcc.org/api/tumorTypes/tree")
  
  
  oncotree_entries <- data.frame()
  
  for (tissue in names(ot_tumor_types_tree$TISSUE$children)) {
    
    ## LEVEL 1 - PRIMARY TUMORS, NOS (TISSUE)
    tree_entry_level1 <- ot_tumor_types_tree$TISSUE$children[[tissue]]
    
    level1_df <- get_oncotree_entry_df(tree_entry = tree_entry_level1)
    code_path_level1 <- level1_df$code
    level1_df$code_path <- code_path_level1
    
    oncotree_entries <- oncotree_entries |>
      dplyr::bind_rows(level1_df)
    
    for (subtype in names(tree_entry_level1$children)) {
      
      ## LEVEL 2 - TUMOR SUBTYPES
      tree_entry_level2 <- tree_entry_level1$children[[subtype]]
      
      level2_df <- get_oncotree_entry_df(tree_entry = tree_entry_level2)
      code_path_level2 <- paste(code_path_level1, 
                                level2_df$code, sep = "-")
      level2_df$code_path <- code_path_level2
      
      oncotree_entries <- oncotree_entries |>
        dplyr::bind_rows(level2_df)
      
      for (subtype2 in names(tree_entry_level2$children)) {
        
        ## LEVEL 3 - TUMOR SUBTYPES
        tree_entry_level3 <- tree_entry_level2$children[[subtype2]]
        
        level3_df <- get_oncotree_entry_df(
          tree_entry = tree_entry_level3)
        
        code_path_level3 <- paste(code_path_level2, 
                                  level3_df$code, sep = "-")
        
        level3_df$code_path <- code_path_level3
        
        oncotree_entries <- oncotree_entries |>
          dplyr::bind_rows(level3_df)
        
        for (subtype3 in names(tree_entry_level3$children)) {
          
          ## LEVEL 4 - TUMOR SUBTYPES
          tree_entry_level4 <- tree_entry_level3$children[[subtype3]]
          
          level4_df <- get_oncotree_entry_df(
            tree_entry = tree_entry_level4
          )
          code_path_level4 <- paste(code_path_level3, 
                                    level4_df$code, sep = "-")
          
          level4_df$code_path <- code_path_level4
          
          oncotree_entries <- oncotree_entries |>
            dplyr::bind_rows(level4_df)
          
          for (subtype4 in names(tree_entry_level4$children)) {
            
            ## LEVEL 4 - TUMOR SUBTYPES
            tree_entry_level5 <- tree_entry_level4$children[[subtype4]]
            
            level5_df <- get_oncotree_entry_df(
              tree_entry = tree_entry_level5)
            code_path_level5 <- paste(code_path_level4, 
                                      level5_df$code, sep = "-")
            
            level5_df$code_path <- code_path_level5
            
            oncotree_entries <- oncotree_entries |>
              dplyr::bind_rows(level5_df)
            
            for (subtype5 in names(tree_entry_level5$children)) {
              
              ## LEVEL 4 - TUMOR SUBTYPES
              tree_entry_level6 <- tree_entry_level5$children[[subtype5]]
              
              level6_df <- get_oncotree_entry_df(
                tree_entry = tree_entry_level6)
              
              code_path_level6 <- paste(code_path_level5, 
                                        level6_df$code, sep = "-")
              
              level6_df$code_path <- code_path_level6
              
              oncotree_entries <- oncotree_entries |>
                dplyr::bind_rows(level6_df)
              
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
  onco_tree <- oncotree_entries |>
    dplyr::arrange(tissue, level) |>
    dplyr::filter(code != "SRCCR") |>
    dplyr::mutate(tissue =
                    dplyr::if_else(tissue == "Bowel",
                                   "Colon/Rectum",as.character(tissue))) |>
    dplyr::mutate(main_type =
                    dplyr::if_else(main_type == "Non-Hodgkin Lymphoma",
                                   "Lymphoma Non_Hodgkin, NOS",
                                   as.character(main_type))) |>
    dplyr::mutate(main_type =
                    dplyr::if_else(main_type == "Hodgkin Lymphoma",
                                   "Lymphoma Hodgkin, NOS",
                                   as.character(main_type))) |>
    
    dplyr::mutate(name =
                    dplyr::if_else(name == "Bowel","Colon/Rectum",
                                   as.character(name))) |>
    dplyr::mutate(main_type =
                    dplyr::if_else(main_type == "Bowel Cancer" |
                                     main_type == "Small Bowel Cancer" |
                                   main_type == "Colorectal Cancer",
                                   "Colorectal Cancer, NOS",
                                   as.character(main_type))) |>
    dplyr::mutate(cui = dplyr::if_else(
      !is.na(nci_t) & nci_t == "C9305",
      "C1334557",
      as.character(cui)
    ))
  
  cui_name_map_lower <- cui_name_map |>
    dplyr::mutate(cui_name_lc = tolower(cui_name))
  
  umls_map_lower <- umls_map$concept |>
    dplyr::mutate(cui_name_lc = tolower(cui_name))
  
  main_type_without_NOS_regex <- paste0(
    "^(Prostate Cancer|Ampullary Carcinoma|Peripheral Nervous System Cancer|",
    "Penile Cancer|Thyroid Cancer|Biliary Tract Cancer|",
    "Blood Cancer|Pancreatic Cancer|Soft Tissue Cancer|Esophageal/Stomach Cancer|",
    "Kidney Cancer|Liver Cancer|Lung Cancer|Other Cancer|Pleural Cancer|Peritoneal Cancer|",
    "Head and Neck Cancer|Hepatobiliary Cancer|Hodgkin Lymphoma|Medulloblastoma|",
    "Leukemia|Melanoma|Mesothelioma|Glioma|Germ Cell Tumor|Testicular Cancer|Uterine Cancer|",
    "Endometrial Cancer|Colorectal Cancer|CNS/Brain Cancer|Skin Cancer|",
    "Esophagaeal/Stomach Cancer|Eye Cancer|Cervical Cancer|Vulvar/Vaginal Cancer|",
    "Ovarian/Fallopian Tube Cancer|Thymic Cancer|Breast Cancer|Bone Cancer)$")

  
  onco_tree_custom_cui_mappings <- 
    read.table(file = "data-raw/oncotree/oncotree.manual_curation.txt",
               header = F, stringsAsFactors = F, sep = "\t") |>
    magrittr::set_names(c("code","cui"))
  
  onco_tree_cui_matched <- onco_tree |>
    dplyr::filter(!is.na(cui))
  
  onco_tree_cui_matched_custom <- onco_tree |>
    dplyr::select(-cui) |>
    dplyr::inner_join(
      onco_tree_custom_cui_mappings, 
      by = "code", multiple = "all",
      relationship = "many-to-many")
  
  onco_tree_cui_matched <- onco_tree_cui_matched |>
    dplyr::anti_join(onco_tree_cui_matched_custom, 
                     by = "code") |>
    dplyr::bind_rows(onco_tree_cui_matched_custom) |>
    dplyr::left_join(cui_name_map, 
                     by = c("cui"), multiple = "all",
                     relationship = "many-to-many")
  
  onco_tree_name_matched <- onco_tree |>
    dplyr::anti_join(
      onco_tree_cui_matched, 
      by = "code") |>
    dplyr::select(-cui) |>
    dplyr::mutate(name_lc = tolower(
      stringr::str_replace(name,"(, NOS)$|(, Other)$",""))) |>
    dplyr::mutate(name_lc = stringr::str_replace_all(
        name_lc, "^mds with", "myelodysplastic syndrome with")) |>
    dplyr::mutate(name_lc = stringr::str_replace_all(
        name_lc, "^aml with",
        "acute myeloid leukemia with")) |>
    dplyr::mutate(name_lc = stringr::str_replace_all(
      name_lc,"-grade "," grade ")) |>
    dplyr::inner_join(
      cui_name_map_lower, 
      by = c("name_lc" = "cui_name_lc"), 
      multiple = "all", relationship = "many-to-many") |>
    dplyr::select(-c(name_lc))
  
  
  remain <- onco_tree |>
    dplyr::anti_join(
      dplyr::bind_rows(onco_tree_cui_matched,
                       onco_tree_name_matched),
      by = "code") |>
    dplyr::select(-cui) |>
    dplyr::mutate(lc_name = tolower(name))
  
  umls_map_fuzzyjoin_candidates <- 
    umls_map$concept |> 
    dplyr::select(cui, cui_name, main_term) |> 
    dplyr::mutate(lc_name = tolower(cui_name)) |> 
    dplyr::filter(nchar(lc_name) >= min(nchar(remain$lc_name)) - 2 &
                    nchar(lc_name) <= max(nchar(remain$lc_name)) + 2)
  
  
  
  onco_tree_fuzzyjoin_matched <- onco_tree |>
    dplyr::anti_join(
      dplyr::bind_rows(onco_tree_cui_matched,
                       onco_tree_name_matched),
      by = "code") |>
    dplyr::select(-cui) |>
    dplyr::mutate(lc_name = tolower(name)) |>
    dplyr::filter(nchar(lc_name) > 3) |>
    fuzzyjoin::stringdist_inner_join(
      umls_map_fuzzyjoin_candidates, by = "lc_name", method = "lv", 
      max_dist = 2, distance_col = "distance") |> 
    dplyr::select(-lc_name.y) |>
    dplyr::distinct() |>
    dplyr::filter(!stringr::str_detect(lc_name.x,"^renal")) |>
    dplyr::filter(!stringr::str_detect(lc_name.x,"lymphoblastic") |
                    (stringr::str_detect(lc_name.x, "lymphoblastic") &
                       stringr::str_sub(tolower(lc_name.x),1,1) == 
                       stringr::str_sub(tolower(cui_name),1,1))) |>
    dplyr::select(-c(lc_name.x,main_term,cui_name,distance)) |>
    dplyr::filter(!(code_path == "MYELOID-MNM-ALAL-MPALBNOS" & 
                      (cui == "C2826055" | cui == "C4726606"))) |>
    dplyr::filter(!(code_path == "MYELOID-MNM-ALAL-MPALTNOS" & 
                      (cui == "C3472616" | cui == "C4726606"))) |>
    dplyr::filter(!(code_path == "MYELOID-MNM-MPN-MPNU" & 
                      (cui == "C1333046" | cui == "CN294221"))) |>
    dplyr::filter(!(code_path == "LYMPH-LNM-NHL-MBN-AHCD" & 
                      cui == "C5209276")) |>
    dplyr::filter(!(code_path == "LYMPH-LNM-BLL" & 
                      cui == "C1292769")) |>
    dplyr::filter(!(code_path == "TESTIS-NSGCT" &
                      cui == "C1336708")) |>
    dplyr::filter(!(code_path == "TESTIS-TSCST" &
                      cui == "C0600113")) |>
    dplyr::filter(!(code_path == "LYMPH-LNM-NHL-MTNN-SEBVTLC" &
                      cui == "C4303422")) |>
    dplyr::filter(!(code == "GHCD" & cui == "C5209277")) |>
    dplyr::filter(!(code == "ADMA" & cui == "C0002448")) |>
    dplyr::filter(!(code == "AMLNOS" & cui == "CN200094")) |>
    dplyr::filter(!(code == "MLNPDGFRA" & cui == "C2827361")) |>
    dplyr::filter(!(code == "MLNPDGFRB" & cui == "C2827360")) |>
    dplyr::filter(!(code == "PMFPES" & cui == "C1516553")) |>
    dplyr::filter(!(code == "SPB" & cui == "C0032131")) |>
    dplyr::filter(!(code == "TMN" & cui == "CN294567")) |>
    dplyr::distinct() |>
    dplyr::left_join(
      cui_name_map, by = "cui", multiple = "all",
      relationship = "many-to-many")
  
  all_oncotree_entries <- onco_tree_cui_matched |>
    dplyr::bind_rows(
      onco_tree_name_matched,
      onco_tree_fuzzyjoin_matched) |>
    dplyr::distinct() |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Adrenal Gland","C0750887",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Ampulla of Vater","C0262401",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Biliary Tract","C0005426",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Bladder/Urinary Tract","C0042076",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Bone","C0005967",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Colon/Rectum","C0009402",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Breast","C1458155",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Cervix","C0007847",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "CNS/Brain","C0153633",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Esophagus/Stomach","C0152018",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Eye","C0496836",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Head and Neck","C0018671",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Kidney","C0022665",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Liver","C0023903",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Lung","C0024121",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Lymphoid","C0024299",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Ovary/Fallopian Tube","C0919267",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Other","C0220647",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Myeloid","C0023418",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Pancreas","C0235974",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Penis","C0153601",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Peripheral Nervous System","C0031118",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Peritoneum","C0153467",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Pleura","C0153494",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Prostate","C0376358",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Skin","C0007114",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Soft Tissue","C4551686",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Testis","C0153594",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Thymus","C0751552",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Thyroid","C0007115",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Uterus","C0153574",as.character(cui))) |>
    dplyr::mutate(cui = dplyr::if_else(
      name == "Vulva/Vagina","C0375071",as.character(cui)))
  
  
  level1_mappings <- all_oncotree_entries |>
    dplyr::filter(level == 1) |>
    dplyr::select(-cui_name) |>
    dplyr::left_join(
      cui_name_map, by = c("cui"), multiple = "all",
      relationship = "many-to-many")
  
  all_oncotree_entries_final <- all_oncotree_entries |>
    dplyr::filter(level > 1) |>
    dplyr::bind_rows(level1_mappings) |> 
    dplyr::arrange(code_path) |>
    dplyr::rename(ot_level = level, ot_code = code) |>
    dplyr::left_join(
      main_types_minor, 
      by = c("main_type"), multiple = "all",
      relationship = "many-to-many") |>
    dplyr::mutate(minor_type = dplyr::if_else(
      is.na(minor_type),
      FALSE,
      as.logical(minor_type)
    )) |>
    dplyr::distinct() |>
    dplyr::mutate(main_type = dplyr::if_else(
      stringr::str_detect(main_type, main_type_without_NOS_regex),
      paste0(main_type,", NOS"),
      as.character(main_type))) |>
    dplyr::mutate(main_type = dplyr::if_else(
      main_type == "Bladder/Urinary Tract Cancer" | 
        main_type == "Bladder Cancer",
      "Bladder/Urinary Tract Cancer, NOS",
      as.character(main_type))) |>
    dplyr::mutate(main_type = dplyr::if_else(
      main_type == "CNS Cancer" | main_type == "",
      "CNS/Brain Cancer, NOS",
      as.character(main_type))) |>
    dplyr::mutate(main_type = dplyr::if_else(
      main_type == "Miscellaneous Neuroepithelial Tumor",
      "CNS/Brain Cancer, NOS",
      as.character(main_type))) |>
    dplyr::mutate(main_type = dplyr::if_else(
      main_type == "Esophagogastric Cancer",
      "Esophageal/Stomach Cancer, NOS",
      as.character(main_type))) |>
    dplyr::mutate(main_type = dplyr::if_else(
      main_type == "Ovarian Cancer",
      "Ovarian/Fallopian Tube Cancer, NOS",
      as.character(main_type))) |>
    dplyr::mutate(main_type = dplyr::if_else(
      main_type == "Thymic Tumor",
      "Thymic Cancer, NOS",
      as.character(main_type))) |>
    dplyr::mutate(main_type = dplyr::if_else(
      main_type == "Medulloblastoma with Extensive Nodularity" | 
        main_type == "Desmoplastic/Nodular Medulloblastoma" | 
        main_type == "Large Cell/Anaplastic Medulloblastoma",
      "Medulloblastoma, NOS",
      as.character(main_type))) |>
    dplyr::mutate(main_type = dplyr::if_else(
      main_type == "Soft Tissue Sarcoma",
      "Soft Tissue Cancer, NOS",
      as.character(main_type))) |>
    dplyr::mutate(main_type = dplyr::if_else(
      main_type == "Vulvar Carcinoma" | main_type == "Vaginal Cancer",
      "Vulvar/Vaginal Cancer, NOS",
      as.character(main_type))) |>
    dplyr::mutate(main_type = dplyr::if_else(
      main_type == "Peripheral Nervous System" | 
        main_type == "Nerve Sheath Tumor",
      "Peripheral Nervous System Cancer, NOS",
      as.character(main_type))) |>
    dplyr::mutate(main_type = dplyr::if_else(
      main_type == "Renal Cell Carcinoma",
      "Kidney Cancer, NOS",
      as.character(main_type))) |>
    dplyr::mutate(main_type = dplyr::if_else(
      main_type == "Retinoblastoma",
      "Eye Cancer, NOS",
      as.character(main_type))) |>
    dplyr::mutate(main_type = dplyr::if_else(
      main_type == "Ampullary Cancer",
      "Ampullary Carcinoma, NOS",
      as.character(main_type))) |>
    dplyr::mutate(main_type = dplyr::if_else(
      main_type == "Lymphatic Cancer",
      "Lymphatic Cancer, NOS",
      as.character(main_type))) |>
    dplyr::mutate(main_type = stringr::str_replace(
      main_type,",NOS",", NOS")) |>
    dplyr::mutate(cui = dplyr::if_else(
      ot_code == "CUP","C0027667",as.character(cui))) |>
    dplyr::arrange(main_type) |>
    dplyr::filter(!is.na(cui_name)) |>
    dplyr::filter(ot_code != "PANET") |>
    dplyr::filter(ot_code != "UCCC") |>
    dplyr::filter(ot_code != "BRCNOS") |>
    dplyr::filter(ot_code != "UASC") |>
    dplyr::filter(ot_code != "ACN") |>
    dplyr::select(-nci_t) |>
    dplyr::distinct() |>
    dplyr::rename(primary_site = tissue, ot_main_type = main_type) |>
    dplyr::distinct() |>
    dplyr::rename(ot_code_path = code_path) |>
    dplyr::rename(ot_name = name) |>
    dplyr::mutate(source = "oncotree_basic") |>
    dplyr::select(primary_site, ot_main_type, ot_name,
                  ot_level, ot_code, ot_code_path,
                  cui, cui_name, dplyr::everything()) |>
    dplyr::arrange(primary_site, ot_main_type, ot_level)
  
  remain <- onco_tree |>
    dplyr::anti_join(
      all_oncotree_entries_final, by = c("code" = "ot_code"))
  
  n_mapped_externally <- 
    nrow(onco_tree_name_matched) + nrow(onco_tree_fuzzyjoin_matched) + 
    nrow(onco_tree_cui_matched_custom)
  n_mapped_internally <- nrow(oncotree_entries) - n_mapped_externally
  
  n_nonmapped <- nrow(onco_tree) - nrow(all_oncotree_entries_final)
  cat("OncoTree entries - UMLS-mapped internally (OncoTree): ", 
      n_mapped_internally, '\n')
  cat("OncoTree entries - UMLS-mapped by phenOncoX (by name matching (exact/fuzzy) and manual curation): ", n_mapped_externally, '\n')
  cat("OncoTree entries - NOT UMLS mapped: ", 
      n_nonmapped, '\n')
  
  oncotree_curated <- all_oncotree_entries_final
  
  umls_icd10_mapping <- read.table(
    gzfile("data-raw/disgenet/disease_mappings.tsv.gz"), 
    sep = "\t", header = T, stringsAsFactors = F, quote="") |> 
    dplyr::rename(cui = diseaseId, 
                  cui_name = name, 
                  icd10_code = code) |> 
    dplyr::filter(vocabulary == "ICD10") |>
    dplyr::select(cui, icd10_code) |>
    dplyr::distinct()
  
  
  syndromes_ding <-
    openxlsx::read.xlsx("data-raw/other/ding_pancancer_cell_2018.xlsx",
                        sheet = 4,startRow = 1) |>
    tidyr::separate_rows(cui_syndrome,sep = ";") |>
    dplyr::select(cui_syndrome) |>
    dplyr::filter(!is.na(cui_syndrome)) |>
    dplyr::rename(cui = cui_syndrome) |>
    dplyr::mutate(cui = stringr::str_trim(cui), 
                  name = "Hereditary Cancer Syndrome, NOS") |>
    dplyr::distinct()
  
  susceptibility_ding <-
    openxlsx::read.xlsx("data-raw/other/ding_pancancer_cell_2018.xlsx",
                        sheet = 4,startRow = 1) |>
    tidyr::separate_rows(cui_susceptibility, sep = ";") |>
    dplyr::select(cui_susceptibility) |>
    dplyr::filter(!is.na(cui_susceptibility)) |>
    dplyr::rename(cui = cui_susceptibility) |>
    dplyr::mutate(cui = stringr::str_trim(cui),
                  name = "Hereditary Cancer Susceptibility, NOS") |>
    dplyr::filter(cui != "C0238198" & 
                    cui != "C0027819" & 
                    cui != "C1333989") |>
    dplyr::distinct()
  
  hereditary_cancers <- openxlsx::read.xlsx(
    "data-raw/other/tumor_nos_umls.xlsx",sheet = 1) |>
    dplyr::mutate(cui = stringr::str_trim(cui)) |>
    dplyr::filter(name == "Hereditary_Cancer_Susceptibility_NOS") |>
    dplyr::mutate(name = "Hereditary Cancer Susceptibility, NOS") |>
    dplyr::bind_rows(syndromes_ding, susceptibility_ding) |>
    dplyr::mutate(ot_main_type = name) |>
    dplyr::distinct() |>
    dplyr::select(-c(tcga_cohort, tmb_high, name, cosmic_mutational_signatures))
  
  other_cuis <- as.data.frame(readr::read_tsv(
    file = "data-raw/other/missing_cui.tsv", na = c("."),
    show_col_types = F
  ))
  
  oncotree_plus_hereditary <- dplyr::bind_rows(
    all_oncotree_entries_final, hereditary_cancers, other_cuis) |>
    dplyr::mutate(minor_type = dplyr::if_else(
      is.na(minor_type),FALSE,as.logical(minor_type))) |>
    dplyr::select(primary_site, ot_main_type, ot_name,
                   ot_level, ot_code, ot_code_path, cui, 
                   minor_type, primary_site, 
                   source) |>
    dplyr::filter(cui != "C0008479" | cui == "C0008479" & 
                    ot_main_type == "Bone Cancer, NOS") |>
    dplyr::filter(cui != "C0279672" | cui == "C0279672" & 
                    ot_main_type == "Soft Tissue Cancer, NOS") |>
    dplyr::filter(cui != "C0029463" | cui == "C0029463" & 
                    ot_main_type == "Bone Cancer, NOS") |>
    dplyr::filter(cui != "C0023467" | cui == "C0023467" & 
                    ot_main_type == "Leukemia, NOS") |>
    dplyr::filter(cui != "C1266144" | cui == "C1266144" & 
                    ot_main_type == "Hereditary Cancer Susceptibility, NOS") |>
    dplyr::filter(cui != "C0008497" | cui == "C0008497" & 
                    ot_main_type == "Germ Cell Tumor, NOS") |>
    dplyr::filter(cui != "C0031511" | cui == "C0031511" & 
                    ot_main_type == "Hereditary Cancer Susceptibility, NOS") |>
    dplyr::filter(cui != "C0035335" | cui == "C0035335" & 
                    ot_main_type == "Eye Cancer, NOS") |>
    dplyr::filter(cui != "C0334663" | cui == "C0334663" & 
                    ot_main_type == "Soft Tissue Cancer, NOS") |>
    dplyr::filter(cui != "C1266111" | cui == "C1266111" & 
                    ot_main_type == "Soft Tissue Cancer, NOS") |>
    dplyr::filter(cui != "C1332564" | cui == "C1332564" & 
                    ot_main_type == "Bladder/Urinary Tract Cancer, NOS") |>
    dplyr::filter(cui != "C0007134" | cui == "C0007134" & 
                    ot_main_type == "Kidney Cancer, NOS") |>
    dplyr::filter(cui != "C0334524" | cui == "C0334524" & 
                    primary_site  == "Vulva/Vagina") |>
    dplyr::filter(cui != "C0346185" | cui == "C0346185" & 
                    primary_site  == "Ovary/Fallopian Tube") |>
    dplyr::filter(cui != "C0751291" | (cui == "C0751291" & 
                                         ot_code == "DMBLNOS")) |>
    dplyr::filter(cui != "C1334970" | 
                    (cui == "C1334970" & ot_code == "MBENNOS")) |>
    dplyr::filter(cui != "C0008497" | 
                    (cui == "C0008497" & ot_code == "BCCA")) |>
    dplyr::filter(cui != "C0023467" | 
                    (cui == "C0023467" & ot_code == "AML")) |>
    dplyr::filter(cui != "C0025149" | 
                    (cui == "C0025149" & ot_code == "MBLNOS")) |>
    dplyr::filter(cui != "C0238196" | 
                    (cui == "C0238196" & ot_code == "SIC")) |>
    dplyr::filter(cui != "C0262401" | 
                    (cui == "C0262401" & ot_code == "AMPULLA_OF_VATER")) |>
    dplyr::filter(cui != "C0279392" | 
                    (cui == "C0279392" & ot_code == "OSMCA")) |>
    dplyr::filter(cui != "C0346185" | 
                    (cui == "C0346185" & ot_code == "ODYS")) |>
    dplyr::filter(cui != "C0687150" | 
                    (cui == "C0687150" & ot_code == "PTH")) |>
    dplyr::filter(cui != "C0948750" | 
                    (cui == "C0948750" & ot_code == "SACA")) |>
    dplyr::filter(cui != "C1266111" | 
                    (cui == "C1266111" & ot_code == "MGST")) |>
    dplyr::filter(cui != "C1458155" | 
                    (cui == "C1458155" & ot_code == "BREAST")) |>
    #dplyr::filter(cui_name != "urachal carcinoma") |>
    dplyr::filter(cui != "C1334970" | 
                    (cui == "C1334970" & ot_code == "MBENNOS")) |>
    dplyr::filter(cui != "C1518872" | 
                    (cui == "C1518872" & ot_code == "PACT")) |>
    dplyr::filter(cui != "C1332166" | 
                    (cui == "C1332166" & ot_code == "EGC")) |>
    dplyr::filter(cui != "C0338113" | 
                    (cui == "C0338113" & ot_code == "USARC")) |>
    dplyr::filter(cui != "C2239246" | 
                    (cui == "C2239246" & ot_code == "UUS")) |>
    dplyr::mutate(source = dplyr::if_else(
      is.na(source), 
      "hereditary_and_extra",
      as.character(source))) |>
    dplyr::distinct() |>
    dplyr::arrange(ot_main_type)
  
  oncotree_plus_hereditary_expanded <- data.frame()
  i <- 1
  for (c in oncotree_plus_hereditary$cui) {
    if (is.na(c)) {
      next
    }
    
    minor_type <- unique(
      oncotree_plus_hereditary[!is.na(oncotree_plus_hereditary$cui) & 
                                 oncotree_plus_hereditary$cui == c,]$minor_type)
    ot_main_type <- unique(
      oncotree_plus_hereditary[!is.na(oncotree_plus_hereditary$cui) & 
                                 oncotree_plus_hereditary$cui == c,]$ot_main_type)
    primary_site <- unique(
      oncotree_plus_hereditary[!is.na(oncotree_plus_hereditary$cui) & 
                                 oncotree_plus_hereditary$cui == c,]$primary_site)
    ot_name_c <- unique(
      oncotree_plus_hereditary[!is.na(oncotree_plus_hereditary$cui) & 
                                 oncotree_plus_hereditary$cui == c,]$ot_name)
    ot_code_c <- unique(
      oncotree_plus_hereditary[!is.na(oncotree_plus_hereditary$cui) & 
                                 oncotree_plus_hereditary$cui == c,]$ot_code)
    ot_level_c <- unique(
      oncotree_plus_hereditary[!is.na(oncotree_plus_hereditary$cui) & 
                                 oncotree_plus_hereditary$cui == c,]$ot_level)
    ot_code_path_c <- unique(
      oncotree_plus_hereditary[!is.na(oncotree_plus_hereditary$cui) & 
                                 oncotree_plus_hereditary$cui == c,]$ot_code_path)
    
   
    
    cat(i, c, ot_code_c, 
        ot_main_type, 
        nrow(oncotree_plus_hereditary_expanded), 
        sep = " - ")
    cat("\n")
    
    if (length(ot_code_c) > 1) {
      ot_code_c <- ot_code_c[1]
      ot_code_path_c <- ot_code_path_c[1]
      ot_level_c <- ot_level_c[1]
      ot_name_c <- ot_name_c[1]
    }
    
    cancer_subtypes <- data.frame()
    if(!is.na(primary_site)){
        
      cancer_subtypes <- 
        get_umls_children(c, 
                          umls_map = umls_map, 
                          umls_cui2name = cui_name_map) |>
        dplyr::mutate(primary_site = primary_site,
                      ot_main_type = ot_main_type,
                      minor_type = minor_type,
                      ot_name = NA,
                      ot_level = NA,
                      ot_code = NA,
                      ot_code_path = NA) |>
        dplyr::mutate(source = "oncotree_basic_hereditary_extra_expanded") |>
        dplyr::select(primary_site, 
                      ot_main_type,
                      ot_name,
                      ot_level,
                      ot_code,
                      ot_code_path,
                      cui,
                      cui_name,
                      minor_type,
                      source) |>
        dplyr::mutate(ot_name = dplyr::if_else(
          cui == c,
          as.character(ot_name_c),
          as.character(ot_name)
        )) |>
        dplyr::mutate(ot_level = dplyr::if_else(
          cui == c,
          as.character(ot_level_c),
          as.character(ot_level)
        )) |>
        dplyr::mutate(ot_code = dplyr::if_else(
          cui == c,
          as.character(ot_code_c),
          as.character(ot_code)
        )) |>
        dplyr::mutate(ot_code_path = dplyr::if_else(
          cui == c,
          as.character(ot_code_path_c),
          as.character(ot_code_path)
        ))
    }else{
      cancer_subtypes <- 
        data.frame(
          'primary_site' = NA,
          'ot_main_type' = ot_main_type,
          'ot_name' = ot_name_c,
          'ot_level' = as.character(ot_level_c),
          'ot_code' = ot_code_c,
          'ot_code_path' = ot_code_path_c,
          'cui' = c,
          'minor_type' = FALSE,
          'source' = 'oncotree_basic_hereditary_extra_expanded'
        ) |>
        dplyr::left_join(
          cui_name_map, by = "cui"
        )
    }
      
    oncotree_plus_hereditary_expanded <- 
      dplyr::bind_rows(
        oncotree_plus_hereditary_expanded, 
        cancer_subtypes) |>
      dplyr::distinct()
    i <- i + 1
  }
  
  cancer_phenotypes_regex_lc <- paste0(
    "tumor|cancer|carcinoma|leukemia|teratoma|seminoma|lymphoma|melanoma|",
    "myeloma|sarcoma|glioma|cholangiocar|medulloblastom|",
    "neoplasm|mesotheliom|(metastasis (from|to))|brca")
  
  non_cancer_terms_regex_lc <- paste0(
    "screening|declined|suspect|history|^rat |response of|progression of|",
    "extends|limited to|recurrence|hamster|porcine|horse|baboon|rhesus|rabbit|",
    "bovine|marmoset|xiphophoru|medaka|fish |tumor size|seen by|infiltration|",
    "category|related to|chicken| pig |canine|protection|poor risk|good risk|",
    "affecting|complicating|configuration|involving|involves|confined|",
    "number of|percent of|cancer diagnosis|therapy|assessed|surgery|",
    "surgical|( in (remission|relapse))|cancer-(related|associated)| the rat |",
    "benign |uncertain behavior|surgery|mouse|invades|",
    "finding|^t[0-9]|^p(t|m|(x|[0-9]{1,}))")
  
  cancer_terms <- cui_name_map_lower |>
    dplyr::filter(
      stringr::str_detect(cui_name_lc, cancer_phenotypes_regex_lc)) |>
    dplyr::filter(
      !stringr::str_detect(cui_name_lc, non_cancer_terms_regex_lc))
  
  missing_cancer_terms <- cancer_terms |>
    dplyr::anti_join(
      oncotree_plus_hereditary_expanded, by = "cui") |>
    dplyr::mutate(ot_main_type = NA, 
                  source = "oncotree_basic_hereditary_extra_expanded",
                  primary_site = NA,
                  minor_type = F)
  
  cancer_term_map <- read.table(
    file = "data-raw/term2tissue.tsv", sep = "\t",
    stringsAsFactors = F, header = T)
  
  i <- 1
  while (i <= nrow(cancer_term_map)) {
    term_regex <- cancer_term_map[i,]$tissue_regex
    ot_main_type_i <- cancer_term_map[i,]$ot_main_type
    primary_site_i <- cancer_term_map[i,]$primary_site
    
    if (cancer_term_map[i,]$extra_cancer_term_needed == 1) {
      missing_cancer_terms <- missing_cancer_terms |>
        dplyr::mutate(ot_main_type = dplyr::if_else(
          is.na(ot_main_type) & stringr::str_detect(cui_name_lc, term_regex) &
            stringr::str_detect(cui_name_lc, "neoplasm|tumor|cancer|carcinoma"),
          as.character(ot_main_type_i),
          as.character(ot_main_type)
        ))
      missing_cancer_terms <- missing_cancer_terms |>
        dplyr::mutate(primary_site = dplyr::if_else(
          is.na(primary_site) & stringr::str_detect(cui_name_lc, term_regex) &
            stringr::str_detect(cui_name_lc, "neoplasm|tumor|cancer|carcinoma"),
          as.character(primary_site_i),
          as.character(primary_site)
        ))
    } else{
      missing_cancer_terms <- missing_cancer_terms |>
        dplyr::mutate(ot_main_type = dplyr::if_else(
          is.na(ot_main_type) & stringr::str_detect(cui_name_lc, term_regex),
          as.character(ot_main_type_i),
          as.character(ot_main_type)
        ))
      missing_cancer_terms <- missing_cancer_terms |>
        dplyr::mutate(primary_site = dplyr::if_else(
          is.na(primary_site) & stringr::str_detect(cui_name_lc, term_regex),
          as.character(primary_site_i),
          as.character(primary_site)
        ))
    }
    i <- i + 1
  }
  
  missing_cancer_terms <- missing_cancer_terms |> 
    dplyr::filter(!is.na(primary_site)) |>
    dplyr::select(-cui_name_lc)
  
  oncotree_plus_hereditary_expanded <- oncotree_plus_hereditary_expanded |>
    dplyr::bind_rows(missing_cancer_terms)
  
  oncotree_plus_hereditary_expanded_final <- oncotree_plus_hereditary_expanded |>
    dplyr::filter(cui != "C0920349") |> # Short Limb Dwarfism-Saddle Nose-Spinal Alterations-Metaphyseal Striation Syndrome
    dplyr::filter(cui != "C0009324") |> # Ulcerative colitis
    dplyr::filter(cui != "C0566602") |> # Primary sclerosing cholangitis
    dplyr::filter(cui != "C0010346") |> # Crohn disease
    dplyr::filter(cui != "C0341332") |> # Indeterminate colitis
    dplyr::filter(cui != "C0021390") |> # Inflammatory bowel disease
    dplyr::filter(cui != "C0014527") |> # Epidermolysis bullosa
    dplyr::filter(cui != "C0008029") |> # Fibrous dysplasia of jaw
    dplyr::filter(cui != "C3854181") |> # Nevus sebaceous
    dplyr::filter(cui != "C1839840") |> # 46,XY sex reversal 8
    dplyr::filter(cui != "C0221026") |> # X-linked agammaglobulinemia
    dplyr::filter(cui != "C1846545") |> # Autoimmune lymphoproliferative syndrome type 2B
    dplyr::filter(cui != "C0238339") |> # Hereditary pancreatitis
    dplyr::filter(
      ot_main_type != "Peripheral Nervous System Cancer, NOS" |
        (ot_main_type == "Peripheral Nervous System Cancer, NOS" & 
           stringr::str_detect(
             tolower(cui_name),
             "schwannom|neuroblastom|neurofibr|neuroendocrine|sheath|nerve|neuro|pheochromo|paraganglio"))) |>
    dplyr::filter(!stringr::str_detect(
      tolower(cui_name),
      "diamond-blackfan|emochromatosis|hrombocytopenia| mononucleosis|performance status")) |>
    dplyr::filter(!stringr::str_detect(
      tolower(cui_name),
      "sarcoidosis|myelofibrosis|neuromyelitis|myelodysplasia")) |>
    dplyr::filter(!stringr::str_detect(
      tolower(cui_name),
      "familial pterygium|polycystic kidney| bullosa|nodular goiter")) |>
    dplyr::filter(!stringr::str_detect(
      tolower(cui_name),
      "chalazion|neutropenia|nevus| lipidosis| myelitis| sclerosus")) |>
    dplyr::filter(!stringr::str_detect(
      tolower(cui_name),
      "minimum residual disease|myasthenia gravis|alpha-1-antitrypsin")) |>
    dplyr::filter(!stringr::str_detect(
      tolower(cui_name),
      "aplastic anemia|atrophic gastritis|leukoplakia|agammaglobulinemia|dyskeratosis")) |>
    dplyr::filter(!stringr::str_detect(
      tolower(cui_name),
      "hypertension|erythroplakia|heel spur|grannuloma annulare| goiter")) |>
    dplyr::filter(!stringr::str_detect(
      tolower(cui_name),
      "parapsoriasis | encephalitis | progression| heart disease| hyperplasia")) |>
    dplyr::filter(!(stringr::str_detect(
      cui_name,"^Leio|Angioleiomyoma| Leio| leio") & 
        stringr::str_detect(primary_site,"Esophagus"))) |>
    dplyr::filter(!(stringr::str_detect(
      tolower(cui_name),"prostate|prostatic") & 
        stringr::str_detect(primary_site,"Head"))) |>
    dplyr::filter(!(stringr::str_detect(
      cui_name,"(O|o)varian|(T|t)esticular|seminoma") & 
        stringr::str_detect(primary_site,"CNS/Brain"))) |>
    dplyr::filter(!(stringr::str_detect(
      tolower(cui_name),"breast|extrahepatic|ampulla of|bile duct") & 
        stringr::str_detect(primary_site,"Colon/Rectum"))) |>
    dplyr::filter(!(stringr::str_detect(
      cui_name,"Wilms") & 
        stringr::str_detect(primary_site,"CNS|Bladder"))) |>
    dplyr::filter(!(stringr::str_detect(
      cui_name,"(K|k)idney") & 
        stringr::str_detect(primary_site,"Bladder"))) |>
    dplyr::filter(!(stringr::str_detect(
      cui_name,"^Lung| lung |^Vulva|^Vaginal|^Verruc") & 
        stringr::str_detect(primary_site,"Skin"))) |>
    dplyr::filter(!(stringr::str_detect(
      tolower(cui_name), "lympho|lympha") & 
        !stringr::str_detect(tolower(cui_name),"skin|papulosis") & 
        (!is.na(primary_site) & primary_site == "Skin"))) |>
    dplyr::filter(!(stringr::str_detect(
      tolower(cui_name), "ewing sarcom") & 
        !stringr::str_detect(tolower(cui_name),"neuro") & 
        (!is.na(primary_site) & primary_site == "CNS/Brain"))) |>
    dplyr::filter(!(stringr::str_detect(
      tolower(cui_name), "sarcoma") & 
        (!is.na(primary_site) & 
           primary_site == "Myeloid" | primary_site == "Lymphoid"))) |>
    dplyr::filter(!(stringr::str_detect(
      tolower(cui_name), "renal cell") & 
        (!is.na(primary_site) & primary_site == "Bladder/Urinary Tract"))) |>
    dplyr::filter(!(stringr::str_detect(
      tolower(cui_name), " sinus ") & 
        (!is.na(primary_site) & primary_site == "Bone"))) |>
    dplyr::filter(!(stringr::str_detect(
      tolower(cui_name), "head and neck") & 
        (!is.na(primary_site) & primary_site == "Skin"))) |>
    dplyr::filter(!(stringr::str_detect(
      tolower(cui_name), "hepatocellular") & 
        (!is.na(primary_site) & primary_site == "CNS/Brain"))) |>
    dplyr::filter(!(stringr::str_detect(
      tolower(cui_name), "pharyngeal|nasopharyn|oropharyng|glotti| sinus | sinus|laryngeal|paranasal|olfactory|thymic|mesothelioma|pharynx") & 
        (!is.na(primary_site) & primary_site == "Lung"))) |>
    dplyr::arrange(ot_main_type) |>
    dplyr::mutate(ot_main_type = stringr::str_replace_all(ot_main_type," |/|, ","_")) |>
    dplyr::mutate(ot_main_type = stringr::str_replace_all(ot_main_type,"_and_","_And_")) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      primary_site == "Other","Other/Unknown",
      as.character(primary_site))) |>
    #dplyr::mutate(primary_site = ot_tissue) |>
    dplyr::mutate(ot_main_type =
                    dplyr::if_else(stringr::str_detect(tolower(cui_name),
                                                       "thyroid") &
                                     !is.na(primary_site) &
                                     primary_site == "Head and Neck",
                                   "Thyroid_Cancer_NOS",
                                   as.character(ot_main_type))) |>
    dplyr::mutate(
      primary_site =
        dplyr::if_else(stringr::str_detect(tolower(cui_name),
                                           "thyroid") &
                         !is.na(primary_site) &
                         primary_site == "Head and Neck",
                       "Thyroid",
                       as.character(primary_site))) |>
    dplyr::mutate(
      primary_site =
        dplyr::if_else(stringr::str_detect(tolower(cui_name),
                                           "thyroid") &
                         primary_site == "Head and Neck" &
                         !is.na(primary_site),
                       "Thyroid",
                       as.character(primary_site))) |>
    dplyr::mutate(
      primary_site =
        dplyr::if_else(stringr::str_detect(tolower(cui_name),
                                           "ovarian|ovary") &
                         !is.na(primary_site) &
                         (primary_site == "Testis" | 
                            primary_site == "CNS/Brain"),
                       "Ovary/Fallopian Tube",
                       as.character(primary_site))) |>
    dplyr::mutate(
      primary_site =
        dplyr::if_else(stringr::str_detect(tolower(cui_name),
                                           "ovarian|ovary") &
                         !is.na(primary_site) &
                         (primary_site == "Testis" | 
                            primary_site == "CNS/Brain"),
                       "Ovary/Fallopian Tube",
                       as.character(primary_site))) |>
    
    dplyr::mutate(
      primary_site =
        dplyr::if_else(!stringr::str_detect(tolower(cui_name),
                                           "testicular|testis|seminom") &
                         !is.na(primary_site) &
                         primary_site == "Testis",
                       "Other/Unknown",
                       as.character(primary_site))) |>
    dplyr::mutate(
      primary_site =
        dplyr::if_else(!stringr::str_detect(tolower(cui_name),
                                           "testicular|testis|seminom") &
                         !is.na(primary_site) &
                         primary_site == "Testis",
                       "Other/Unknown",
                       as.character(primary_site))) |>
    
    
    dplyr::mutate(
      ot_main_type =
        dplyr::if_else(stringr::str_detect(tolower(cui_name),
                                           "barrett") &
                         primary_site == "Head and Neck" &
                         !is.na(primary_site),
                       "Esophageal_Stomach_Cancer_NOS",
                       as.character(ot_main_type))) |>
    dplyr::mutate(
      primary_site =
        dplyr::if_else(stringr::str_detect(tolower(cui_name),
                                           "barrett") &
                         primary_site == "Head and Neck"&
                         !is.na(primary_site),
                       "Esophagus/Stomach",
                       as.character(primary_site))) |>
    dplyr::mutate(
      primary_site =
        dplyr::if_else(stringr::str_detect(tolower(cui_name),
                                           "waldenstr") &
                         is.na(primary_site),
                       "Lymphoid",
                       as.character(primary_site))) |>
    dplyr::mutate(
      primary_site =
        dplyr::if_else(stringr::str_detect(tolower(cui_name),
                                           "barrett") &
                         primary_site == "Head and Neck" &
                         !is.na(primary_site),
                       "Esophagus/Stomach",
                       as.character(primary_site))) |>
    dplyr::filter(ot_main_type != "Adenocarcinoma_In_Situ") |>
    dplyr::filter(!(stringr::str_detect(tolower(cui_name),"liver") &
                      primary_site != "Liver")) |>
    dplyr::filter(!(stringr::str_detect(tolower(cui_name),"central nervous system") &
                      primary_site != "CNS/Brain")) |>
    dplyr::filter(!(stringr::str_detect(tolower(cui_name),"gastric|esophag") &
                      primary_site != "Esophagus/Stomach")) |>
    dplyr::filter(!(stringr::str_detect(tolower(cui_name),"lung") &
                      primary_site != "Lung")) |>
    dplyr::filter(!(stringr::str_detect(tolower(cui_name),"pancreatic|pancreas") &
                      primary_site != "Pancreas")) |>
    dplyr::filter(!(stringr::str_detect(tolower(cui_name),"testicular|testis") &
                      primary_site != "Testis")) |>
    dplyr::filter(!(stringr::str_detect(tolower(cui_name),"ovarian|ovary|fallopian") &
                      primary_site != "Ovary/Fallopian Tube")) |>
    dplyr::filter(!(stringr::str_detect(tolower(cui_name),"thyroid") &
                      primary_site != "Thyroid")) |>
    dplyr::filter(!(ot_main_type != "CNS_Brain_Cancer_NOS" &
                      stringr::str_detect(tolower(cui_name),"meningioma"))) |>
    dplyr::filter(!(ot_main_type == "Skin_Cancer_NOS" &
                      stringr::str_detect(tolower(cui_name),"esophag"))) |>
    dplyr::filter(!((ot_main_type == "Bone_Cancer_NOS" | 
                       ot_main_type == "Soft_Tissue_Cancer_NOS") &
                      stringr::str_detect(tolower(cui_name), 
                                          "neuroblastom|glioblastom"))) |>
    dplyr::filter(!(ot_main_type == "Head_And_Neck_Cancer_NOS" &
                      stringr::str_detect(tolower(cui_name),"esophag"))) |>
    dplyr::filter(!(ot_main_type == "Miscellaneous_Brain_Tumor" &
                      !stringr::str_detect(
                        tolower(cui_name),
                        "brain|glioblastom|medulloblastom|neuroblastom|retinoblastom|pineoblastom"))) |>
    dplyr::filter(!(ot_main_type == "Hereditary_Cancer_Susceptibility_NOS" &
                      !stringr::str_detect(
                        tolower(cui_name),
                        "familial|hereditary|xeroderma|increased risk|fanconi|susceptibility|polyposis|brca|colorectal cancer [1-9]|neuroblastoma [1-9]|pancreatic cancer [1-9]|syndrome"))) |>
    dplyr::filter(!(ot_main_type == "Cancer_of_Unknown_Primary" & 
                      !stringr::str_detect(tolower(cui_name), 
                                           "mixed|(unknown (primary|origin))"))) |>
    dplyr::mutate(ot_main_type = dplyr::if_else(
      stringr::str_detect(tolower(cui_name), 
                          "hereditary|susceptibility|familial"),
      "Hereditary_Cancer_Susceptibility_NOS",
      as.character(ot_main_type)
    )) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      stringr::str_detect(
        tolower(cui_name),"hereditary|susceptibility|familial"),
      as.character(NA),
      as.character(primary_site)
    )) |>
    dplyr::filter(!stringr::str_detect(
      tolower(cui_name),
      "(((high|intermediate|low|increased) risk)|presence of|anatomic location)")) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      stringr::str_detect(tolower(cui_name),"lymphoma") &
        !is.na(primary_site) &
        primary_site != "Lymphoid",
      "Lymphoid",
      as.character(primary_site)
    )) |>
   
    dplyr::mutate(ot_main_type = dplyr::if_else(
      stringr::str_detect(tolower(cui_name),"lymphoma") &
        !is.na(primary_site) &
        primary_site != "Lymphoid",
      "Lymphatic_Cancer_NOS",
      as.character(ot_main_type)
    )) |>
    dplyr::mutate(primary_site = dplyr::if_else(
      stringr::str_detect(tolower(cui_name),"leukemia") &
        !is.na(primary_site) &
        primary_site != "Myeloid",
      "Myeloid",
      as.character(primary_site)
    )) |>
    dplyr::mutate(ot_main_type = dplyr::if_else(
      stringr::str_detect(tolower(cui_name),"leukemia") &
        !is.na(primary_site) &
        primary_site != "Myeloid",
      "Blood_Cancer_NOS",
      as.character(ot_main_type)
    )) |>
    dplyr::distinct() 
    
  missing_hereditary <- dplyr::bind_rows(
    data.frame("cui" = "C0027672", 
               "ot_main_type" = "Hereditary_Cancer_Syndrome_NOS",
               "source" = "oncotree_basic_hereditary_extra_expanded", 
               minor_type = F),
    data.frame("cui" = "C1333600", 
               "ot_main_type" = "Hereditary_Cancer_Susceptibility_NOS",
               "source" = "oncotree_basic_hereditary_extra_expanded", 
               minor_type = F),
    data.frame("cui" = "C1333990", 
               "ot_main_type" = "Hereditary_Cancer_Susceptibility_NOS",
               "source" = "oncotree_basic_hereditary_extra_expanded", 
               minor_type = F),
    data.frame("cui" = "C1846758", 
               "ot_main_type" = "Hereditary_Cancer_Susceptibility_NOS",
               "source" = "oncotree_basic_hereditary_extra_expanded",
               minor_type = F),
    data.frame("cui" = "C0042138", 
               "ot_main_type" = "Uterine_Cancer_NOS", 
                primary_site = "Uterus",
               "source" = "oncotree_basic_hereditary_extra_expanded", 
               minor_type = F),
    data.frame("cui" = "C0014170", 
               "ot_main_type" = "Uterine_Cancer_NOS", 
               primary_site = "Uterus",
               "source" = "oncotree_basic_hereditary_extra_expanded", 
               minor_type = F),
    data.frame("cui" = "C0376544", 
               "ot_main_type" = "Lymphatic_Cancer_NOS", 
               primary_site = "Lymphoid",
               "source" = "oncotree_basic_hereditary_extra_expanded", 
               minor_type = F),
    data.frame("cui" = "C0039590", 
               "ot_main_type" = "Testicular_Cancer_NOS", 
               primary_site = "Testis",
               "source" = "oncotree_basic_hereditary_extra_expanded", 
               minor_type = F),
    data.frame("cui" = "C0031149", 
               "ot_main_type" = "Peritoneal_Cancer_NOS", 
               primary_site = "Peritoneum",
               "source" = "oncotree_basic_hereditary_extra_expanded", 
               minor_type = F),
    data.frame("cui" = "C0238301", 
               "ot_main_type" = "Head_And_Neck_Cancer_NOS", 
               primary_site = "Head and Neck",
               "source" = "oncotree_basic_hereditary_extra_expanded", 
               minor_type = F)) |>
    dplyr::left_join(cui_name_map, by = "cui", multiple = "all",
                     relationship = "many-to-many")
    
  oncotree_plus_hereditary_expanded_final <- 
    oncotree_plus_hereditary_expanded_final |>
    dplyr::bind_rows(missing_hereditary)
  
  
  # remove overlapping entries in hereditary cancer syndrome group and 
  # cancer susceptibility group
  hereditary_condition_entries <-
    dplyr::filter(oncotree_plus_hereditary_expanded_final, 
                  ot_main_type == "Hereditary_Cancer_Susceptibility_NOS")
  hereditary_syndrome_entries <-
    dplyr::filter(oncotree_plus_hereditary_expanded_final, 
                  ot_main_type == "Hereditary_Cancer_Syndrome_NOS")
  hereditary_condition_entries <-
    dplyr::anti_join(hereditary_condition_entries,
                     dplyr::select(hereditary_syndrome_entries, cui),
                     by = "cui")
  sporadic_cancer_entries <-
    dplyr::filter(oncotree_plus_hereditary_expanded_final, 
                  ot_main_type != "Hereditary_Cancer_Susceptibility_NOS" &
                    ot_main_type != "Hereditary_Cancer_Syndrome_NOS")
  
  oncoterms_final <- dplyr::bind_rows(
    hereditary_condition_entries,
    hereditary_syndrome_entries,
    sporadic_cancer_entries) |>
    dplyr::filter(!is.na(cui)) |>
    dplyr::arrange(
      primary_site, 
      ot_main_type, 
      ot_name,
      ot_level, 
      ot_code, 
      ot_code_path,
    )
  
  # final_tree_slim <- final_tree |> 
  #   dplyr::filter(minor_type == F)
  onco_pheno_map <- 
    list("oncotree_expanded" = oncoterms_final)
  
  ## Cross-reference OncoTree with EFO and DO
  for (m in c('oncotree_expanded')) {
    onco_pheno_map[[m]] <- onco_pheno_map[[m]] |>
      dplyr::filter(cui_name != "urachal carcinoma") |>
      dplyr::left_join(
        dplyr::select(efo_map$efo2xref, cui, 
                      efo_id,  efo_name),
        by = "cui", multiple = "all", 
        relationship = "many-to-many") |>
      dplyr::left_join(
        do_map, by = "cui", multiple = "all", 
        relationship = "many-to-many") |>
      dplyr::left_join(
        icd10_map, by = "cui", multiple = "all", 
        relationship = "many-to-many") |>
      dplyr::filter(!(cui == "C0006826" & efo_id == "EFO:0000311")) |>
      dplyr::filter(!(cui == "C0023467" & efo_id == "MONDO:0015667")) |>
      dplyr::filter(!(cui == "C0017075" & efo_id == "MONDO:0016730")) |>
      dplyr::filter(!(cui == "C0042133" & efo_id == "EFO:0000731")) |>
      dplyr::filter(!(cui == "C0021071" & efo_id == "EFO:1001798")) |>
      dplyr::filter(!(cui == "C0027708" & efo_id == "Orphanet:654")) |>
      dplyr::distinct() |>
      dplyr::mutate(
        efo_id =
          dplyr::if_else(
            ot_main_type == "Colorectal_Cancer_NOS" &
              stringr::str_detect(tolower(cui_name),
                                  "metasta"),
            "EFO:1001480",as.character(efo_id))) |>
      dplyr::mutate(
        efo_name =
          dplyr::if_else(
            ot_main_type == "Colorectal_Cancer_NOS" &
              stringr::str_detect(tolower(cui_name),
                                  "metasta"),
            "metastatic colorectal cancer",
            as.character(efo_name)))
    
    ## MAP NON-MAPPED EFO IDENTIFIERS BY NAME
    efo_map_lc <- efo_map$efo2name |>
      dplyr::mutate(efo_name_lc = tolower(efo_name))
    
    tmp1 <- onco_pheno_map[[m]] |>
      dplyr::filter(!is.na(primary_site) & is.na(efo_id)) |>
      dplyr::select(-c(efo_name, efo_id)) |>
      dplyr::mutate(cui_name_lc = tolower(cui_name)) |>
      dplyr::left_join(
        dplyr::filter(efo_map_lc, !is.na(primary_site)),
        by = c("cui_name_lc" = "efo_name_lc", 
               "primary_site" = "primary_site"), 
        multiple = "all", relationship = "many-to-many") |>
      dplyr::select(-cui_name_lc)
    
    tmp2 <- onco_pheno_map[[m]] |>
      dplyr::filter(is.na(primary_site) & is.na(efo_id)) |>
      dplyr::select(-c(efo_name, efo_id, primary_site)) |>
      dplyr::mutate(cui_name_lc = tolower(cui_name)) |>
      dplyr::left_join(
        dplyr::filter(efo_map_lc, is.na(primary_site)),
        by = c("cui_name_lc" = "efo_name_lc"), 
        multiple = "all", relationship = "many-to-many") |>
      dplyr::select(-cui_name_lc)
    
    tmp3 <- onco_pheno_map[[m]] |>
      dplyr::filter(
        (!is.na(primary_site) & !is.na(efo_id)) |
          (is.na(primary_site) & !is.na(efo_id)))
    
    onco_pheno_map[[m]] <- 
      dplyr::bind_rows(tmp1,
                       tmp2,
                       tmp3) |>
      dplyr::filter(
        !(!is.na(efo_name) & 
            stringr::str_detect(
              tolower(efo_name),"wilms tumor") &
            (!is.na(primary_site) & 
               stringr::str_detect(
              primary_site, "CNS/Brain|Bladder/Urinary Tract"
            ))
        )) |>
      dplyr::filter(
        !(!is.na(efo_name) & 
            stringr::str_detect(
              tolower(efo_name),"ewing sarcoma") &
            (!is.na(primary_site) & 
               stringr::str_detect(
                 primary_site, "CNS/Brain|Peripheral Nervous System"
               ))
        )) |>
      dplyr::filter(
        !(!is.na(efo_name) & 
            stringr::str_detect(
              tolower(efo_name),"hidradenoma|papilloma|ocular |uveal |oropharynx|laryngeal|hypopharyngeal|oral squamous") &
            (!is.na(primary_site) & 
               stringr::str_detect(
                 primary_site, "Skin"
               ))
        )) |>
      dplyr::filter(
        !(!is.na(efo_name) & 
            stringr::str_detect(
              tolower(efo_name),"retinoblastoma") &
            (!is.na(primary_site) & 
               stringr::str_detect(
                 primary_site, "Bone|CNS/Brain|Soft Tissue"
               ))
        )) |>
      
      dplyr::filter(
        !(!is.na(cui_name) & 
            stringr::str_detect(
              tolower(cui_name),
              paste0(
                "(colon|colorectal|anal|anus|rectum|rectal|ileal|",
                "jejunal|polyp|peritoneum|pancreat|hepato|",
                "lymphoma|appendiceal|intestinal|sigmoid|cecum|",
                "fibrolamellar|alpha-heavy|klatskin|hepatic|gardner|cecal|",
                "biliary|hepatocellular|bile duct|duodenal|",
                "cholangio|lymphatic|",
                "gallbladder|cholangiocarcinoma|hepato|hepatoblastom|peritoneal|",
                "ampulla of vater|intestine|mesotheliom|polyposis|appendix)")) &
            (!is.na(primary_site) &
               primary_site == "Esophagus/Stomach"))) |>
      
      dplyr::filter(
        !(!is.na(efo_name) & 
            stringr::str_detect(
              tolower(efo_name),"squamous cell") &
            (!is.na(primary_site) &
               primary_site == "Skin" &
               !(stringr::str_detect(tolower(efo_name), 
                                      "keratin|skin ") | 
                   tolower(efo_name) == "squamous cell carcinoma" |
                 tolower(efo_name) == "squamous cell carcinoma in situ"
                )
            )
        )) |> 
      dplyr::filter(
        !(!is.na(efo_name) & 
            stringr::str_detect(
              tolower(efo_name),"glottis|tonsil|mesothelioma|tracheal|laryngeal|thymus|thymic|thymoma|pharynx|pharyngeal") &
            (!is.na(primary_site) &
               primary_site == "Lung"
            )
        )) |>
      dplyr::filter(
        !(((!is.na(efo_name) & 
            stringr::str_detect(
              tolower(efo_name),"prostate")) |
            (!is.na(cui_name) &
               stringr::str_detect(
                 tolower(cui_name),"colorectal"))) &
            (stringr::str_detect(
              tolower(efo_name),"acinar ") &
            (!is.na(primary_site) &
               (primary_site == "Head and Neck" | 
                  primary_site == "Cervix")
            ))
        )) |>
      dplyr::filter(
        !(!is.na(efo_name) & 
            stringr::str_detect(
              tolower(efo_name)," sinus ") &
            (!is.na(primary_site) &
               (primary_site == "Soft Tissue" | 
                  primary_site == "Lung" |
                  primary_site == "Bone" |
                  primary_site == "Lung")
            )
        )) |>
      dplyr::filter(
        !(!is.na(efo_name) & 
            stringr::str_detect(
              tolower(efo_name),"pituitary gland|pituitary tumor") &
            (!is.na(primary_site) &
               (primary_site == "Head Neck")
            )
        )) |>
      dplyr::filter(
        !(!is.na(efo_name) & 
            stringr::str_detect(
              tolower(efo_name),"peritoneal") &
            (!is.na(primary_site) &
               (primary_site == "Pleura")
            )
        )) |>
      dplyr::filter(
        !(!is.na(efo_name) & 
            stringr::str_detect(
              tolower(efo_name),"^(malignant )?renal |collecting duct|hepatoblastoma|pleuropulmonary|kidney|seminoma|nephr") &
            (!is.na(primary_site) &
               (primary_site == "Bladder/Urinary Tract" |
                  primary_site == "CNS/Brain")
            )
        )) |>
      dplyr::mutate(primary_site = dplyr::if_else(
          !is.na(efo_name) & 
            stringr::str_detect(
              tolower(efo_name),"heart (cancer|neoplasm)|neoplasm of heart") &
            (!is.na(primary_site) & primary_site == "Lung"),
          "Other/Unknown",
          as.character(primary_site)
        )) |>
      dplyr::mutate(primary_site = dplyr::if_else(
        !is.na(efo_name) & 
          stringr::str_detect(
            tolower(efo_name)," syndrome$") &
         !is.na(primary_site),
        as.character(NA),
        as.character(primary_site)
      )) |>
      dplyr::mutate(primary_site = dplyr::if_else(
        !is.na(efo_name) & 
          stringr::str_detect(
            tolower(efo_name),"heart (cancer|neoplasm)|neoplasm of heart") &
          (!is.na(primary_site) & 
             (primary_site == "Other/Unknown" | primary_site == "Lung")),
        as.character(NA),
        as.character(primary_site)
      )) |>
      dplyr::mutate(primary_site = dplyr::if_else(
        is.na(primary_site) &
          ot_main_type != "Hereditary_Cancer_Susceptibility_NOS" & 
          stringr::str_detect(
            tolower(cui_name), "myelodysplastic syndrome"
          ),
        "Myeloid",
        as.character(primary_site)
      )) |>
      dplyr::filter(!(
        stringr::str_detect(tolower(cui_name), "heart") &
          (!is.na(ot_main_type) & ot_main_type == "Lung_Cancer_NOS"))
      ) |>
      dplyr::filter(
        !(stringr::str_detect(tolower(cui_name), "cardiac"))) |>
      dplyr::filter(
        !(stringr::str_detect(tolower(cui_name), "non-hereditary"))
      ) |>
      dplyr::select(
        primary_site,
        ot_main_type,
        ot_name,
        ot_level,
        ot_code,
        ot_code_path,
        cui, 
        cui_name,
        efo_id,
        efo_name,
        do_id,
        do_name,
        do_cancer_slim,
        minor_type,
        icd10_code,
        source,
        dplyr::everything()
      ) |>
      dplyr::arrange(
        primary_site,
        ot_main_type,
        ot_level,
        ot_code_path,
        cui_name
      )
  }
  
  onco_pheno_map$oncotree <- oncotree_curated

  return(onco_pheno_map)
  
}

get_umls_children <- function(
  c, 
  umls_map = NULL, 
  umls_cui2name = NULL) {
  children_level0 <- data.frame("cui" = c, stringsAsFactors = F) |>
    dplyr::inner_join(umls_cui2name, by = c("cui"))
  cui_children <- dplyr::filter(umls_map$relation, cui == c) |>
    dplyr::select(cui_child)
  
  ## LEVEL 1
  children_level1 <- data.frame()
  if (nrow(cui_children) > 0) {
    children_level1 <- data.frame("cui" = cui_children$cui_child,
                                  "level" = 1, stringsAsFactors = F) |>
      dplyr::inner_join(umls_cui2name, by = c("cui")) |>
      dplyr::select(-level)
  }
  
  
  ## LEVEL 2
  children_level2 <- data.frame()
  for (c1 in unique(children_level1$cui)) {
    children <- dplyr::filter(umls_map$relation, cui == c1)
    if (nrow(children) > 0) {
      children <- children |>
        dplyr::select(cui_child) |>
        dplyr::rename(cui = cui_child)
      children_level2 <- dplyr::bind_rows(children_level2, children) |>
        dplyr::mutate(level = 2)
    }
  }
  if (nrow(children_level2) > 0) {
    children_level2 <- children_level2 |>
      dplyr::select(-level) |>
      #dplyr::distinct() |>
      dplyr::inner_join(umls_cui2name, by = c("cui"))
  }
  
  ## LEVEL 3
  children_level3 <- data.frame()
  for (c2 in unique(children_level2$cui)) {
    children <- dplyr::filter(umls_map$relation, cui == c2)
    if (nrow(children) > 0) {
      children <- children |>
        dplyr::select(cui_child) |>
        dplyr::rename(cui = cui_child)
      children_level3 <- dplyr::bind_rows(children_level3, children) |>
        dplyr::mutate(level = 3)
    }
  }
  if (nrow(children_level3) > 0) {
    children_level3 <- children_level3 |>
      dplyr::select(-level) |>
      dplyr::distinct() |>
      dplyr::inner_join(umls_cui2name, by = c("cui"))
  }
  
  
  ## LEVEL 4
  children_level4 <- data.frame()
  for (c3 in unique(children_level3$cui)) {
    children <- dplyr::filter(umls_map$relation, cui == c3)
    if (nrow(children) > 0) {
      children <- children |>
        dplyr::select(cui_child) |>
        dplyr::rename(cui = cui_child)
      children_level4 <- dplyr::bind_rows(children_level4, children) |>
        dplyr::mutate(level = 4)
    }
  }
  if (nrow(children_level4) > 0) {
    children_level4 <- children_level4 |>
      dplyr::select(-level) |>
      dplyr::distinct() |>
      dplyr::inner_join(umls_cui2name, by = c("cui"))
  }
  
  ## LEVEL 5
  children_level5 <- data.frame()
  for (c4 in unique(children_level4$cui)) {
    children <- dplyr::filter(umls_map$relation, cui == c4)
    if (nrow(children) > 0) {
      children <- children |>
        dplyr::select(cui_child) |>
        dplyr::rename(cui = cui_child)
      children_level5 <- dplyr::bind_rows(children_level5, children) |>
        dplyr::mutate(level = 5)
    }
  }
  if (nrow(children_level5) > 0) {
    children_level5 <- children_level5 |>
      dplyr::select(-level) |>
      dplyr::distinct() |>
      dplyr::inner_join(umls_cui2name, by = c("cui"))
  }
  
  
  ## LEVEL 5
  children_level6 <- data.frame()
  for (c5 in unique(children_level5$cui)) {
    children <- dplyr::filter(umls_map$relation, cui == c5)
    if (nrow(children) > 0) {
      children <- children |>
        dplyr::select(cui_child) |>
        dplyr::rename(cui = cui_child)
      children_level6 <- dplyr::bind_rows(children_level6, children) |>
        dplyr::mutate(level = 6)
    }
  }
  if (nrow(children_level6) > 0) {
    children_level6 <- children_level6 |>
      dplyr::select(-level) |>
      dplyr::distinct() |>
      dplyr::inner_join(umls_cui2name, by = c("cui"))
  }
  
  all_children <- children_level0 |>
    dplyr::bind_rows(children_level1) |>
    dplyr::bind_rows(children_level2) |>
    dplyr::bind_rows(children_level3) |>
    dplyr::bind_rows(children_level4) |>
    dplyr::bind_rows(children_level5) |>
    dplyr::bind_rows(children_level6) |>
    dplyr::distinct()
  
  return(all_children)
}

