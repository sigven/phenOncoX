library(magrittr)
library(jsonlite)
library(dplyr)
library(stringr)
library(ontologyIndex)

map_efo <- function(umls_map,
                    efo_release = "v3.33.0",
                    update = T,
                    basedir = NULL){

  if(!dir.exists(basedir)){
    rlogging::message("Base directory does not exist")
  }else{
    if(!file.exists(
      file.path(basedir,"data-raw","efo",paste0("efo.",efo_release,".obo"))) |
      update == T){
      download.file(
        paste0("https://github.com/EBISPOT/efo/releases/download/",
               efo_release,"/efo.obo"),
        destfile =
          file.path(basedir, "data-raw","efo",paste0("efo.",efo_release,".obo"))
      )
    }
  }
  con <- file(
    description =
      file.path(
        basedir,"data-raw","efo",
        paste0("efo.",
               efo_release,".obo")),
    open="r")
  lines <- readLines(con)
  close(con)

  i <- 1
  efo_id <- NA
  nci_t <- NA
  msh_id <- NA
  cui_all <- c()
  cui_close_all <- c()
  cui_exact <- NA
  cui_exact_all <- c()
  msh_all <- c()
  nci_all <- c()
  efo_map <- NULL
  name <- NA
  cui <- NA
  cui_close <- NA
  obsolete <- FALSE
  ancestor <- NA
  ancestors <- c()
  while(i <= length(lines)){
    line <- lines[i]
    #cat(line,'\n')
    ## only include major ontologies (phenotype-related)
    if(stringr::str_detect(line, "^id: (EFO|HP|DOID|MONDO|GO|Orphanet):[0-9]{1,}$")){
      if(!is.na(efo_id) &
         !is.na(name) &
         !stringr::str_detect(name, "measurement|^CS") &
         obsolete == FALSE){
        df <-
          data.frame("efo_id" = efo_id,
                     "nci_t" = paste(unique(nci_all),collapse = ","),
                     "msh" = paste(unique(msh_all),collapse = ","),
                     "cui" = paste(unique(cui_all),collapse = ","),
                     "cui_exact" = paste(unique(cui_exact_all), collapse=","),
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
      cui_close <- NA
      cui_exact <- NA
      cui_all <- c()
      cui_close_all <- c()
      cui_exact_all <- c()
      msh_all <- c()
      nci_all <- c()
      cell_line <- FALSE
      obsolete <- FALSE
      ancestor <- NA
      ancestors <- c()
      efo_id <- stringr::str_replace(line,"^id: ","")
    }
    if(stringr::str_detect(line,"name: .+$")){
      name <- stringr::str_replace(line, "name: ","")
      if("http" %in% name){
        name <- NA
      }
    }
    if(stringr::str_detect(line,"is_obsolete: true$")){
      obsolete <- TRUE
    }
    if(stringr::str_detect(line,"is_a: EFO:[0-9]{1,} ")){
      ancestor <- stringr::str_replace(
        stringr::str_match(line, "EFO:[0-9]{1,} ")[[1]]," ","")
      ancestors <- c(ancestors,ancestor)
    }

    # if(stringr::str_detect(line,"xref: NCI(T|t):[A-Z]{1,2}[0-9]{1,}")){
    #   nci_t <- stringr::str_replace(
    #     stringr::str_match(line,
    #                        "xref: NCI(t|T)[A-Z]{1,2}[0-9]{1,}")[[1]],
    #     "xref: NCI(T|t):","")
    #   nci_all <- c(nci_all,nci_t)
    # }

    if(stringr::str_detect(line,"xref: MSH:[A-Z]{1,2}[0-9]{1,}")){
      msh_id <- stringr::str_replace_all(line,"xref: MSH:","")
      msh_all <- c(msh_all, msh_id)
    }
    if(stringr::str_detect(line,"xref: UMLS:[A-Z]{1,2}[0-9]{1,}")){
      cui <- stringr::str_replace_all(
        stringr::str_match(line,"xref: UMLS:[A-Z]{1,2}[0-9]{1,}")[[1]],
        "xref: UMLS:","")
      cui_all <- c(cui_all, cui)
    }
    if(stringr::str_detect(
      line,
      "property_value: \"closeMatch\" http://linkedlifedata.com/resource/umls/id/")){
      cui_close <- stringr::str_replace_all(
        line,
        "property_value: \"closeMatch\" http://linkedlifedata.com/resource/umls/id/","")
      cui_close_all <- c(cui_close_all, cui_close)
    }
    if(stringr::str_detect(
      line,
      "property_value: exactMatch http://linkedlifedata.com/resource/umls/id/")){
      cui_exact <- stringr::str_replace_all(
        line,
        "property_value: exactMatch http://linkedlifedata.com/resource/umls/id/","")
      cui_exact_all <- c(cui_exact_all, cui_exact)
    }
    if(stringr::str_detect(
      line,
      "property_value: exactMatch NCIT:")){
      nci_t <- stringr::str_replace_all(
        line,
        "property_value: exactMatch NCIT:","")
      nci_all <- c(nci_all, nci_t)
    }
    i <- i + 1
  }

  efo_map <- as.data.frame(
    efo_map %>%
      dplyr::mutate(cui_close =
                      dplyr::if_else(cui_close == "",
                                     as.character(NA),cui_close)) %>%
      dplyr::mutate(nci_t = dplyr::if_else(nci_t == "",
                                           as.character(NA),nci_t)) %>%
      dplyr::mutate(msh = dplyr::if_else(msh == "",
                                         as.character(NA),msh)) %>%
      dplyr::mutate(cui = dplyr::if_else(cui == "",
                                         as.character(NA),cui)) %>%
      dplyr::mutate(cui_exact = dplyr::if_else(cui_exact == "",
                                         as.character(NA),cui_exact)) %>%
      tidyr::separate_rows(cui,sep = ",") %>%
      tidyr::separate_rows(cui_exact,sep = ",") %>%
      tidyr::separate_rows(cui_close,sep = ",") %>%
      tidyr::separate_rows(nci_t,sep = ",") %>%
      tidyr::separate_rows(msh,sep = ",") %>%
      dplyr::mutate(
        efo_name =
          stringr::str_replace(
            efo_name,
            " \\{http://www.co-ode.org/patterns#createdBy=","")) %>%
      dplyr::mutate(
        efo_name =
          stringr::str_replace(
            efo_name,'"http://www.ebi.ac.uk/ontology/webulous#OPPL_pattern"\\}',
            '')) %>%
      dplyr::distinct()
  )

  efo2name <- efo_map %>%
    dplyr::select(efo_id, efo_name) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplasm|tumor") &
          stringr::str_detect(tolower(efo_name), "colon|anal|anus|cecum|rectum|small intestine|rectal|colorectal"),
        "Colon/Rectum",
        as.character(NA)
      )
    )  %>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "bladder|urinary|ureter|urethra|urothelial") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplasm|tumor"),
        "Bladder/Urinary Tract",
        as.character(primary_site)
      )
    )%>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "pancreas|acinar cell|pancreatic") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplasm|tumor"),
        "Pancreas",
        as.character(primary_site)
      )
    )%>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "thymic|thymus") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplasm|tumor"),
        "Thymus",
        as.character(primary_site)
      )
    ) %>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name),"schwannoma") |
        (stringr::str_detect(tolower(efo_name), "neuroectodermal|peripheral nervous|peripheral nerve sheath") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplasm|tumor")),
        "Peripheral Nervous System",
        as.character(primary_site)
      )
    )%>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "thyroid") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplasm|tumor"),
        "Thyroid",
        as.character(primary_site)
      )
    )%>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "leukemia|myeloma") |
          (stringr::str_detect(tolower(efo_name), "t-cell|b-cell") &
             stringr::str_detect(tolower(efo_name), "cancer|tumor|neoplasm")) |
        (stringr::str_detect(tolower(efo_name), "myeloid") &
          stringr::str_detect(tolower(efo_name), "leukemia|neoplasm")),
        "Myeloid",
        as.character(primary_site)
      )
    ) %>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "gastric|stomach|esophagus|esophageal") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor"),
        "Esophagus/Stomach",
        as.character(primary_site)
      )
    ) %>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "eye|ocular |retinal ") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor"),
        "Eye",
        as.character(primary_site)
      )
    ) %>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "gallbladder|biliary tract|bile duct|hepatobiliary|cholangio|biliary intraepithelial") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor"),
        "Biliary Tract",
        as.character(primary_site)
      )
    ) %>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "chondroblastom|chordoma|giant cell") |
        (stringr::str_detect(tolower(efo_name), "bone") &
          stringr::str_detect(tolower(efo_name), "carcinoma|sarcom|cancer|neoplas|tumor")),
        "Bone",
        as.character(primary_site)
      )
    ) %>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
        !stringr::str_detect(tolower(efo_name), "bone") &
          stringr::str_detect(tolower(efo_name), "sarcom|soft tissue|leiomy|connective tissue"),
        "Soft Tissue",
        as.character(primary_site)
      )
    ) %>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "prostat") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor"),
        "Prostate",
        as.character(primary_site)
      )
    ) %>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "ovary|ovarian|fallopian tube") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|blastom|tumor"),
        "Ovary",
        as.character(primary_site)
      )
    ) %>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "uterine|endometri|female reproductive") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor"),
        "Uterus",
        as.character(primary_site)
      )
    ) %>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "cervix|cervical") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor"),
        "Cervix",
        as.character(primary_site)
      )
    ) %>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "head and neck|mouth|neck|glottis|larynx|parotid gland|salivary gland|oral squamous|tongue|nasal cavity|nasopharyngeal|oral cavity|oropharynx") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor"),
        "Head and Neck",
        as.character(primary_site)
      )
    ) %>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "liver|hepatocellular") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor"),
        "Liver",
        as.character(primary_site)
      )
    ) %>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name),"mesotheliom") |
        (stringr::str_detect(tolower(efo_name), "pleura") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|mesotheliom|neoplas|tumor")),
        "Pleura",
        as.character(primary_site)
      )
    ) %>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
          stringr::str_detect(tolower(efo_name), "lung|bronchi|bronchoalveolar") &
            stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor"),
          "Lung",
        as.character(primary_site)
      )
    ) %>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "neuroblastom|meningiom|glioblastom|medulloblastom|glioma|astrocytom") |
          (stringr::str_detect(tolower(efo_name), "central nervous|nervous system|brain") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor")),
        "CNS/Brain",
        as.character(primary_site)
      )
    ) %>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
          stringr::str_detect(tolower(efo_name), "breast") &
             stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor"),
        "Breast",
        as.character(primary_site)
      )
    ) %>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "kidney|renal|wilms") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor"),
        "Kidney",
        as.character(primary_site)
      )
    ) %>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name),"seminoma") |
        (stringr::str_detect(tolower(efo_name), "testis|testicular|embryonal|male reproductive|germ cell") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor")),
        "Testis",
        as.character(primary_site)
      )
    ) %>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name),"melanoma") |
          (stringr::str_detect(tolower(efo_name), "skin") &
             stringr::str_detect(tolower(efo_name), "squamous cell|carcinoma|cancer|neoplas|tumor")),
        "Skin",
        as.character(primary_site)
      )
    ) %>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
          (stringr::str_detect(tolower(efo_name), "vulva|vagina") &
             stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor")),
        "Vulva/Vagina",
        as.character(primary_site)
      )
    ) %>%
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
    ) %>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "lymphom") |
          (stringr::str_detect(tolower(efo_name), "cancer|neoplasm") &
          stringr::str_detect(tolower(efo_name), "lympho|hematopoietic"))
          ,
        "Lymphoid",
        as.character(primary_site)
      )
    ) %>%
    dplyr::mutate(
      primary_site = dplyr::if_else(
        stringr::str_detect(tolower(efo_name), "adrenal gland|adrenal") &
          stringr::str_detect(tolower(efo_name), "carcinoma|cancer|neoplas|tumor|blastom|lipom"),
        "Adrenal Gland",
        as.character(primary_site)
      )
    )
  cui_map <- efo_map %>%
    dplyr::filter(!is.na(cui)) %>%
    dplyr::select(efo_id, cui) %>%
    dplyr::left_join(dplyr::filter(umls_map$concept, main_term == T),
                     by = "cui") %>%
    dplyr::select(efo_id,cui,cui_name) %>%
    dplyr::distinct()

  cui_map_close <- efo_map %>%
    dplyr::filter(!is.na(cui_close)) %>%
    dplyr::select(efo_id, cui_close) %>%
    dplyr::left_join(dplyr::filter(umls_map$concept, main_term == T),
                     by = c("cui_close" = "cui")) %>%
    dplyr::select(efo_id,cui_close,cui_name) %>%
    dplyr::rename(cui = cui_close) %>%
    dplyr::distinct()

  nci_map <- efo_map %>%
    dplyr::filter(!is.na(nci_t)) %>%
    dplyr::select(efo_id, nci_t) %>%
    dplyr::left_join(umls_map$nci, by = "nci_t") %>%
    dplyr::left_join(dplyr::filter(umls_map$concept,
                                   main_term == T), by = "cui") %>%
    dplyr::filter(!is.na(cui)) %>%
    dplyr::select(efo_id,cui,cui_name) %>%
    dplyr::distinct()

  msh_map <- efo_map %>%
    dplyr::filter(!is.na(msh)) %>%
    dplyr::select(efo_id,msh) %>%
    dplyr::left_join(umls_map$msh, by = "msh") %>%
    dplyr::left_join(dplyr::filter(umls_map$concept, main_term == T),
                     by = "cui") %>%
    dplyr::filter(!is.na(cui)) %>%
    dplyr::select(efo_id,cui,cui_name) %>%
    dplyr::distinct()

  efo2xref <- dplyr::bind_rows(cui_map, cui_map_close, msh_map, nci_map) %>%
    dplyr::distinct() %>%
    dplyr::left_join(efo2name, by = "efo_id") %>%
    dplyr::filter(is.na(cui) | (!is.na(cui) & !is.na(cui_name)))

  return(list("efo2name" = efo2name,
              "efo2xref" = efo2xref))
}


map_disease_ontology <- function(
  skip_non_cui_mapped = T,
  release = "2021-07-30",
  umls_map = NULL,
  basedir = NULL){

  do_github_raw_url <-
    paste0("https://raw.githubusercontent.com/DiseaseOntology/",
           "HumanDiseaseOntology/main/src/ontology/releases/")

  release_dest <- stringr::str_replace_all(release,"-","")
  if(!file.exists(
    file.path(basedir,
              "data-raw",
              "do",
              paste0("TopNodes_Docancerslim.", release_dest,".obo")))){
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
  if(!file.exists(
    file.path(basedir, "data-raw","do",
              paste0("doid.",release_dest,".obo")))){
    download.file(
      url =
        paste0(do_github_raw_url, "doid.obo"),
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

  umls_map_lower <- umls_map$concept %>%
    dplyr::mutate(cui_name_lc = tolower(cui_name))

  do_ids <- do_index$id
  do_names <- do_index$name
  names(do_ids) <- NULL
  names(do_names) <- NULL
  do_map <- data.frame()
  do_map_top <- data.frame()
  while(i < length(do_ids)){
    do_id <- do_ids[i]
    if(stringr::str_detect(do_id,"DOID:")){
      do_name <- do_names[i]
      all_ancestors <-
        ontologyIndex::get_term_property(ontology = do_index,
                                         property="ancestors",
                                         term = do_id)
      all_xref <-
        ontologyIndex::get_term_property(ontology = do_index,
                                         property = "xref",
                                         term = do_id)
      all_subset <-
        ontologyIndex::get_term_property(ontology = do_index,
                                         property = "subset",
                                         term = do_id)
      cui <- NA
      do_cancer_slim <- FALSE
      do_rare_slim <- FALSE
      do_cancer_slim_top <- FALSE
      top_node_cancer_do_id <- NA
      top_node_cancer <- data.frame()
      cui_ids <- c()
      if(length(all_xref) > 0){
        j <- 1
        while(j <= length(all_xref)){
          if(stringr::str_detect(all_xref[j],"^UMLS_CUI:C[0-9]{1,}$")){
            cui <- stringr::str_replace(all_xref[j],"^UMLS_CUI:","")
            cui_ids <- c(cui_ids,cui)
          }
          j <- j + 1
        }
      }
      if(length(all_ancestors) > 0){
        all_ancestor_df <-
          data.frame("do_id" = all_ancestors,
                     stringsAsFactors = F)
        top_node_cancer <-
          as.data.frame(
            dplyr::inner_join(all_ancestor_df,
                              do_cancer_top_df,by="do_id") %>%
              dplyr::distinct())
      }
      if(length(all_subset) > 0){
        if("TopNodes_DOcancerslim" %in% all_subset){
          do_cancer_slim_top <- T
        }
        if("DO_cancer_slim" %in% all_subset){
          do_cancer_slim <- T
        }
      }


      if(length(cui_ids) == 0){
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
        while(u <= length(cui_ids)){
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

      if(nrow(top_node_cancer) > 0){
        k <- 1
        while(k <= nrow(top_node_cancer)){
          if(length(cui_ids) == 0){
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
            while(u <= length(cui_ids)){
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
        if(length(cui_ids) == 0){
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
          while(u <= length(cui_ids)){
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

  do_map_cui_name_matched <- do_map %>%
    dplyr::select(do_id,do_name) %>%
    dplyr::left_join(umls_map_lower, by = c("do_name" = "cui_name_lc")) %>%
    dplyr::filter(!is.na(cui)) %>%
    dplyr::select(do_id, cui) %>%
    dplyr::distinct()

  #do_map <- dplyr::left_join(do_map, do_map_cui_name_matched)

  if(skip_non_cui_mapped == T){
    do_map1 <- do_map %>% dplyr::filter(!is.na(cui)) %>% dplyr::distinct()
    do_map2 <- do_map %>% dplyr::filter(is.na(cui)) %>%
      dplyr::select(-c(cui)) %>%
      dplyr::inner_join(do_map_cui_name_matched,by = c("do_id"))

    do_map <- dplyr::bind_rows(do_map1, do_map2)
  }

  ## manual correction of erroneous or missing UMLS cross-references
  do_map <- do_map %>%
    ## sarcoma (non-existent UMLS_CUI = C0153519)
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:1115",
                                       "C1261473",
                                       as.character(cui))) %>%
    ## colorectal cancer (UMLS_CUI = C1261473, not part of tree)
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:9256",
                                       "C0009404", as.character(cui))) %>%
    ## hepatocellular carcinoma (missing UMLS_CUI)
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:684",
                                       "C2239176", as.character(cui))) %>%
    ## chronic leukemia (missing UMLS_CUI)
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:1036",
                                       "C1279296", as.character(cui))) %>%
    ## lymphoid leukemia (missing UMLS_CUI)
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:10747",
                                       "C0023448", as.character(cui))) %>%
    ## brain glioma (missing UMLS_CUI)
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:0060108",
                                       "C0349661", as.character(cui))) %>%
    ## ovarian serious carcinoma (missing UMLS_CUI)
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:0050933",
                                       "C0279663", as.character(cui))) %>%
    ## mucosal melanoma (missing UMLS_CUI)
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:0050929",
                                       "C3898222", as.character(cui))) %>%
    ## acral lentiginous melanoma (missing UMLS_CUI)
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:6367",
                                       "C0346037", as.character(cui))) %>%
    ## brain glioma (missing UMLS_CUI)
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:0080146",
                                       "C0279584",as.character(cui))) %>%
    ## correct urothelial carcinoma (transitional cell carcinoma)
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:2671",
                                       "C2145472", as.character(cui))) %>%
    ## colon cancer
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:219",
                                       "C0699790", as.character(cui))) %>%
    ## stomach cancer
    dplyr::mutate(cui = dplyr::if_else(do_id == "DOID:10534",
                                       "C0699791", as.character(cui))) %>%
    ## skip myelodysplastic syndrome
    dplyr::filter(do_id != "DOID:0050908") %>%

    dplyr::bind_rows(
      data.frame(do_id = "DOID:1788",
                 do_name = "peritoneal mesothelioma",
                 cui = "C1377610", do_cancer_slim = T,
                 stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(do_id = "DOID:1790",
                                do_name = "mesothelioma",
                                cui = "C0025500", do_cancer_slim = T,
                                stringsAsFactors = F)) %>%
    dplyr::bind_rows(
      data.frame(do_id = "DOID:1612", do_name = "breast cancer",
                 cui = "C0678222",  do_cancer_slim = T,
                 stringsAsFactors = F)) %>%
    dplyr::bind_rows(
      data.frame(do_id = "DOID:12603", do_name = "acute leukemia",
                 cui = "C0085669",  do_cancer_slim = T,
                 stringsAsFactors = F)) %>%
    dplyr::bind_rows(
      data.frame(do_id = "DOID:1793", do_name = "pancreatic cancer",
                 cui = "C0030297",  do_cancer_slim = T,
                 stringsAsFactors = F)) %>%
    dplyr::bind_rows(
      data.frame(do_id = "DOID:1793", do_name = "pancreatic cancer",
                 cui = "C0887833",  do_cancer_slim = T,
                 stringsAsFactors = F)) %>%
    dplyr::bind_rows(
      data.frame(do_id = "DOID:1793", do_name = "pancreatic cancer",
                 cui = "C0235974",  do_cancer_slim = T,
                 stringsAsFactors = F)) %>%
    dplyr::bind_rows(
      data.frame(do_id = "DOID:8923", do_name = "skin melanoma",
                 cui = "C0025202",  do_cancer_slim = T,
                 stringsAsFactors = F)) %>%
    dplyr::bind_rows(
      data.frame(do_id = "DOID:1324", do_name = "lung cancer",
                 cui = "C0684249",  do_cancer_slim = T,
                 stringsAsFactors = F)) %>%
    dplyr::bind_rows(
      data.frame(do_id = "DOID:0050908", do_name = "myelodysplastic syndrome",
                 cui = "C3463824",  do_cancer_slim = T,
                 stringsAsFactors = F)) %>%
    dplyr::bind_rows(
      data.frame(do_id = "DOID:9256", do_name = "colorectal cancer",
                 cui = "C0699790", do_cancer_slim = T,
                 stringsAsFactors = F)) %>%
    dplyr::bind_rows(
      data.frame(do_id = "DOID:0060081",
                 do_name = "triple-receptor negative breast cancer",
                 cui = "C3539878", do_cancer_slim = T,
                 stringsAsFactors = F)) %>%
    dplyr::bind_rows(
      data.frame(do_id = "DOID:0060079",
                 do_name = "her2-receptor positive breast cancer",
                 cui = "C1960398", do_cancer_slim = T,
                 stringsAsFactors = F)) %>%
    dplyr::bind_rows(
      data.frame(do_id = "DOID:0050745",
                 do_name = "diffuse large b-cell lymphoma",
                 cui = "C0079744", do_cancer_slim = T,
                 stringsAsFactors = F)) %>%
    dplyr::bind_rows(
      data.frame(do_id = "DOID:0050861",
                 do_name = "colorectal adenocarcinoma",
                 cui = "C0338106", do_cancer_slim = T,
                 stringsAsFactors = F)) %>%
    dplyr::bind_rows(
      data.frame(do_id = "DOID:0060105",
                 do_name = "brain medulloblastoma",
                 cui = "C1332188", do_cancer_slim = T,
                 stringsAsFactors = F)) %>%
    dplyr::bind_rows(
      data.frame(do_id = "DOID:0050744",
                 do_name = "anaplastic large cell lymphoma",
                 cui = "C0206180", do_cancer_slim = T,
                 stringsAsFactors = F)) %>%
    dplyr::bind_rows(
      data.frame(do_id = "DOID:0050746",
                 do_name = "mantle cell lymphoma",
                 cui = "C0334634", do_cancer_slim = T,
                 stringsAsFactors = F)) %>%
    dplyr::bind_rows(
      data.frame(do_id = "DOID:0050749",
                 do_name = "peripheral t-cell lymphoma",
                 cui = "C0079774", do_cancer_slim = T,
                 stringsAsFactors = F)) %>%
    dplyr::bind_rows(
      data.frame(do_id = "DOID:0050873",
                 do_name = "follicular lymphoma",
                 cui = "C0024301", do_cancer_slim = T,
                 stringsAsFactors = F)) %>%
    dplyr::bind_rows(
      data.frame(do_id = "DOID:0070004", do_name = "myeloma",
                 cui = "C0026764", do_cancer_slim = T,
                 stringsAsFactors = F)) %>%
    dplyr::bind_rows(
      data.frame(do_id = "DOID:3965",
                 do_name = "merkel cell carcinoma",
                 cui = "C0007129", do_cancer_slim = T,
                 stringsAsFactors = F))

  return(do_map)
}

map_umls <- function(update = T,
                     basedir = NULL){

  if(is.null(basedir)){
    rlogging::message("Please specifiy valid base directory")
  }
  if(!dir.exists(basedir)){
    rlogging::message("Base directory does not exist")
  }

  for(fn in c('MGCONSO','NAMES','MGREL_1','MGREL_2')){
    if(!file.exists(
      file.path(basedir, "data-raw", "umls", paste0(fn,".csv.gz"))) | update == T){
      download.file(
        paste0("ftp://ftp.ncbi.nlm.nih.gov/pub/medgen/csv/",fn,".csv.gz"),
        destfile = file.path(
          basedir, "data-raw", "umls", paste0(fn,".csv.gz"))
      )
    }
  }


  if(!file.exists(file.path(
    basedir, "data-raw", "umls", "Neoplasm_Core.txt"))){
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
      na.strings = c(""), stringsAsFactors = F) %>%
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
  umls_rel <- rbind(umls_rel_1, umls_rel_2) %>%
    dplyr::filter(REL == "CHD" & SUPPRESS == "N") %>%
    dplyr::select(CUI1,CUI2) %>%
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
      stringsAsFactors = F) %>%
    dplyr::select(CUI,name) %>%
    dplyr::rename(STR = name) %>%
    dplyr::mutate(main_term = TRUE)

  concept_summary_data <-
    read.csv(
      gzfile(
        file.path(basedir,"data-raw","umls","MGCONSO.csv.gz")),
      stringsAsFactors = F) %>%
    dplyr::select(CUI,SAB,STR) %>%
    dplyr::distinct() %>%
    dplyr::left_join(concept_names_main, by = c("CUI","STR")) %>%
    dplyr::rename(cui = CUI, source = SAB, cui_name = STR) %>%
    dplyr::mutate(main_term =
                    dplyr::if_else(is.na(main_term),
                                   FALSE,TRUE,TRUE))

  # saveRDS(concept_summary_data,
  #         file.path(
  #           basedir,"output","cui_umls_data.rds"
  #         )
  # )

  umls_nci <- dplyr::filter(umls, SAB == "NCI") %>%
    dplyr::select(CUI,SCUI) %>%
    dplyr::distinct() %>%
    dplyr::rename(cui = CUI, nci_t = SCUI)

  umls_msh <- dplyr::filter(umls, SAB == "MSH") %>%
    dplyr::select(CUI,SDUI) %>%
    dplyr::distinct() %>%
    dplyr::rename(cui = CUI, msh = SDUI)

  umls <- list("concept" = concept_summary_data,
               "relation" = umls_rel,
               "nciXref" = umls_nci,
               "mshXref" = umls_msh)
  return(umls)
}

onco_pheno_mapping <- function(umls_map){

  main_types_minor <-
    read.table(file="data-raw/oncotree/oncotree.main_types_minor.txt",
               header = F,stringsAsFactors = F,
               quote="",sep="\t") %>%
    dplyr::rename(main_type = V1) %>%
    dplyr::mutate(minor_type = T)

  nci_map <- umls_map$nci %>%
    dplyr::rename(cui_2 = cui)

  cui_name_map <- umls_map$concept %>%
    dplyr::filter(main_term == T) %>%
    dplyr::select(cui,cui_name) %>%
    dplyr::distinct()

  ot_tumor_types <-
    jsonlite::fromJSON("http://oncotree.mskcc.org/api/tumorTypes")
  cui <- data.frame(
    "cui" =
      as.character(ot_tumor_types$externalReferences$UMLS),
    stringsAsFactors = F)
  nci_t <- data.frame(
    "nci_t" =
      as.character(ot_tumor_types$externalReferences$NCI),
    stringsAsFactors = F)
  onco_tree <- data.frame("main_type" = ot_tumor_types$mainType,
                          "name" = ot_tumor_types$name,
                          "code" = ot_tumor_types$code,
                          "tissue" = ot_tumor_types$tissue,
                          "level" = ot_tumor_types$level,
                          "cui" = cui,
                          "nci_t" = nci_t, stringsAsFactors = F) %>%
    dplyr::mutate(cui = dplyr::if_else(cui == "NULL",
                                       as.character(NA),
                                       as.character(cui))) %>%
    dplyr::mutate(nci_t = dplyr::if_else(nci_t == "NULL",
                                         as.character(NA),
                                         as.character(nci_t))) %>%
    dplyr::filter(!is.na(main_type)) %>%
    dplyr::arrange(tissue,main_type,desc(level)) %>%
    dplyr::filter(code != "SRCCR") %>%
    dplyr::mutate(tissue =
                    dplyr::if_else(tissue == "Bowel",
                                   "Colon/Rectum",as.character(tissue))) %>%
    dplyr::mutate(main_type =
                    dplyr::if_else(main_type == "Non-Hodgkin Lymphoma",
                                   "Lymphoma Non_Hodgkin, NOS",
                                   as.character(main_type))) %>%
    dplyr::mutate(main_type =
                    dplyr::if_else(main_type == "Hodgkin Lymphoma",
                                   "Lymphoma Hodgkin, NOS",
                                   as.character(main_type))) %>%

    dplyr::mutate(name =
                    dplyr::if_else(name == "Bowel","Colon/Rectum",
                                   as.character(name))) %>%
    dplyr::mutate(main_type =
                    dplyr::if_else(main_type == "Bowel Cancer, NOS" |
                                     main_type == "Small Bowel Cancer",
                                   "Colorectal Cancer",
                                   as.character(main_type)))

  cui_name_map_lower <- cui_name_map %>%
    dplyr::mutate(cui_name_lc = tolower(cui_name))

  umls_map_lower <- umls_map$concept %>%
    dplyr::mutate(cui_name_lc = tolower(cui_name))

  onco_tree_cui <- onco_tree %>%
    dplyr::filter(!is.na(cui)) %>%
    dplyr::left_join(cui_name_map,by=c("cui"))

  onco_tree_cui1 <- onco_tree_cui %>%
    dplyr::filter(!is.na(cui) & !is.na(cui_name))
  onco_tree_cui2 <- onco_tree_cui %>%
    dplyr::filter(!is.na(cui) & is.na(cui_name)) %>%
    dplyr::select(-cui_name)

  onco_tree_cui_matched <- onco_tree %>%
    dplyr::filter(is.na(cui)) %>%
    dplyr::bind_rows(onco_tree_cui2) %>%
    dplyr::select(-c(cui,nci_t)) %>%
    dplyr::mutate(name_lc =
                    tolower(stringr::str_replace(name,"(, NOS)$|(, Other)$",""))) %>%
    dplyr::mutate(name_lc =
                    stringr::str_replace_all(
                      name_lc, "^mds with", "myelodysplastic syndrome with")) %>%
    dplyr::mutate(name_lc =
                    stringr::str_replace_all(
                      name_lc, "^aml with",
                      "acute myeloid leukemia with")) %>%
    dplyr::mutate(name_lc =
                    stringr::str_replace_all(name_lc,"-grade "," grade ")) %>%
    dplyr::inner_join(cui_name_map_lower, by = c("name_lc" = "cui_name_lc"))  %>%
    dplyr::left_join(nci_map, by = c("cui" = "cui_2")) %>%
    dplyr::select(-c(name_lc))


  found <- dplyr::bind_rows(onco_tree_cui1, onco_tree_cui_matched)
  remain <- onco_tree %>%
    dplyr::anti_join(found, by = c("name")) %>%
    dplyr::select(-c(cui,nci_t)) %>%
    dplyr::mutate(name_lc = tolower(stringr::str_replace(name,"(, NOS)$|(, Other)$",""))) %>%
    dplyr::mutate(name_lc = stringr::str_replace_all(name_lc,"^mds with","myelodysplastic syndrome with")) %>%
    dplyr::mutate(name_lc = stringr::str_replace_all(name_lc,"^aml with","acute myeloid leukemia with")) %>%
    dplyr::mutate(name_lc = stringr::str_replace_all(name_lc,"-grade "," grade ")) %>%
    dplyr::inner_join(umls_map_lower, by = c("name_lc" = "cui_name_lc"))  %>%
    dplyr::select(-name_lc) %>%
    dplyr::left_join(nci_map, by = c("cui" = "cui_2"))

  remain1 <- dplyr::filter(remain, source == "NCI") %>%
    dplyr::select(-source) %>%
    dplyr::distinct()

  remain2 <- remain %>% dplyr::anti_join(remain1, by = c("name")) %>%
    dplyr::filter(name != "Transient Abnormal Myelopoiesis" | (name == "Transient Abnormal Myelopoiesis" & source == "MSH")) %>%
    dplyr::filter(name != "Brenner Tumor, Malignant" | (name == "Brenner Tumor, Malignant" & source == "MSH")) %>%
    dplyr::filter(name != "Neuroendocrine Tumor, NOS" | (name == "Neuroendocrine Tumor, NOS" & source == "MSH")) %>%
    dplyr::select(-source)

  mapped_all <- as.data.frame(dplyr::bind_rows(found, remain1, remain2) %>%
                                dplyr::select(-main_term) %>%
                                dplyr::group_by(main_type, tissue, code, name, level, cui, cui_name) %>%
                                dplyr::summarise(nci_t = paste(nci_t, collapse="&"),
                                                 .groups = "drop"))

  remain <- dplyr::anti_join(onco_tree, mapped_all, by = c("name")) %>%
    dplyr::select(-c(cui,nci_t))

  #onco_tree_mapped <- onco_tree %>% dplyr::filter(!is.na(cui_name))
  level1_mapped <- remain %>%
    dplyr::filter(level == 1) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Adrenal Gland","C0750887",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Breast Invasive Carcinoma, NOS" | name == "Breast Invasive Cancer, NOS","C0853879",
                                       as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Ampulla of Vater","C0262401",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Biliary Tract","C0005426",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Bladder/Urinary Tract","C0042076",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Bone","C0005967",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Colon/Rectum","C0009402",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Breast","C1458155",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Cervix","C0302592",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "CNS/Brain","C0153633",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Esophagus/Stomach","C0152018",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Eye","C0848866",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Head and Neck","C0018671",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Kidney","C0022665",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Liver","C0023903",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Lung","C0024121",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Lymphoid","C0024299",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Ovary/Fallopian Tube","C0919267",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Myeloid","C0023418",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Pancreas","C0235974",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Penis","C0153601",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Peripheral Nervous System","C0031118",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Peritoneum","C0153467",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Pleura","C0153494",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Prostate","C0600139",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Skin","C0037286",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Soft Tissue","C0037579",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Testis","C0039590",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Thymus","C0751552",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Thyroid","C0007115",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Uterus","C0042138",as.character(cui))) %>%
    dplyr::mutate(cui = dplyr::if_else(name == "Vulva/Vagina","C0042237",as.character(cui))) %>%
    dplyr::filter(name != "Other") %>%
    dplyr::left_join(cui_name_map, by = c("cui"))

  mapped_all <- mapped_all %>% dplyr::bind_rows(level1_mapped) %>%
    dplyr::rename(ot_level = level, ot_code = code) %>%
    dplyr::left_join(main_types_minor,by=c("main_type")) %>%
    dplyr::distinct() %>%
    dplyr::mutate(main_type = dplyr::if_else(stringr::str_detect(main_type,"^(Prostate Cancer|Peripheral Nervous System|Penile Cancer|Thyroid Cancer|Pancreatic Cancer|Soft Tissue Cancer|Head and Neck Cancer|Hepatobiliary Cancer|Hodgkin Lymphoma|Leukemia|Melanoma|Mesothelioma|Glioma|Germ Cell Tumor|Endometrial Cancer|Colorectal Cancer|Cervical Cancer|Breast Cancer|Bone Cancer)$"),
                                             paste0(main_type,", NOS"),as.character(main_type))) %>%
    dplyr::mutate(main_type = dplyr::if_else(main_type == "Bladder Cancer","Bladder/Urinary Tract Cancer, NOS",as.character(main_type))) %>%
    dplyr::mutate(main_type = dplyr::if_else(main_type == "CNS Cancer","CNS/Brain Cancer, NOS",as.character(main_type))) %>%
    dplyr::mutate(main_type = dplyr::if_else(main_type == "Miscellaneous Neuroepithelial Tumor","CNS/Brain Cancer, NOS",as.character(main_type))) %>%
    dplyr::mutate(main_type = dplyr::if_else(main_type == "Esophagogastric Cancer","Esophageal/Stomach Cancer, NOS",as.character(main_type))) %>%
    dplyr::mutate(main_type = dplyr::if_else(main_type == "Ovarian Cancer","Ovarian/Fallopian Tube Cancer, NOS",as.character(main_type))) %>%
    dplyr::mutate(main_type = dplyr::if_else(main_type == "Thymic Tumor","Thymic Cancer, NOS",as.character(main_type))) %>%
    dplyr::mutate(main_type = dplyr::if_else(main_type == "Soft Tissue Sarcoma","Soft Tissue Cancer, NOS",as.character(main_type))) %>%
    dplyr::mutate(main_type = dplyr::if_else(main_type == "Vulvar Carcinoma" | main_type == "Vaginal Cancer","Vulvar/Vaginal Cancer, NOS",as.character(main_type))) %>%
    dplyr::mutate(main_type = dplyr::if_else(main_type == "Small Bowel Cancer","Bowel Cancer, NOS",as.character(main_type))) %>%
    dplyr::mutate(main_type = dplyr::if_else(main_type == "Peripheral Nervous System, NOS" | main_type == "Nerve Sheath Tumor","Peripheral Nervous System Cancer, NOS",as.character(main_type))) %>%
    dplyr::mutate(main_type = dplyr::if_else(main_type == "Renal Cell Carcinoma","Kidney Cancer, NOS",as.character(main_type))) %>%
    dplyr::mutate(main_type = dplyr::if_else(main_type == "Retinoblastoma","Eye Cancer, NOS",as.character(main_type))) %>%
    dplyr::mutate(main_type = dplyr::if_else(main_type == "Ampullary Cancer","Ampullary Carcinoma, NOS",as.character(main_type))) %>%
    dplyr::mutate(main_type = stringr::str_replace(main_type,",NOS",", NOS")) %>%
    dplyr::filter(main_type != "Miscellaneous Brain Tumor") %>%
    dplyr::filter(name != "Paraganglioma") %>%
    dplyr::filter(main_type != "Cancer of Unknown Primary" | (main_type == "Cancer of Unknown Primary" & ot_code == "CUP")) %>%
    dplyr::mutate(cui = dplyr::if_else(ot_code == "CUP","C0027667",as.character(cui)))

  onco_tree <- mapped_all %>%
    dplyr::arrange(main_type) %>%
    dplyr::filter(!is.na(cui_name)) %>%
    dplyr::select(-c(nci_t,ot_code,ot_level)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(primary_site = tissue, group = main_type) %>%
    dplyr::filter(cui != "C0346202" |
                    cui == "C0346202" &
                    group == "Cervical Cancer, NOS") %>%
    dplyr::filter(cui != "C1332912" |
                    cui == "C1332912" &
                    group == "Cervical Cancer, NOS") %>%
    dplyr::filter(cui != "C0024299" |
                    cui == "C0024299" &
                    group == "Lymphatic Cancer, NOS") %>%
    dplyr::distinct()

  syndromes_ding <-
    openxlsx::read.xlsx("data-raw/other/ding_pancancer_cell_2018.xlsx",
                        sheet = 4,startRow = 1) %>%
    tidyr::separate_rows(cui_syndrome,sep=";") %>%
    dplyr::select(cui_syndrome) %>%
    dplyr::filter(!is.na(cui_syndrome)) %>%
    dplyr::rename(cui = cui_syndrome) %>%
    dplyr::mutate(cui = stringr::str_trim(cui), name = "Hereditary Cancer Syndrome, NOS") %>%
    dplyr::distinct()

  susceptibility_ding <-
    openxlsx::read.xlsx("data-raw/other/ding_pancancer_cell_2018.xlsx",
                        sheet = 4,startRow = 1) %>%
    tidyr::separate_rows(cui_susceptibility,sep=";") %>%
    dplyr::select(cui_susceptibility) %>%
    dplyr::filter(!is.na(cui_susceptibility)) %>%
    dplyr::rename(cui = cui_susceptibility) %>%
    dplyr::mutate(cui = stringr::str_trim(cui),
                  name = "Hereditary Cancer Susceptibility, NOS") %>%
    dplyr::filter(cui != "C0238198" & cui != "C0027819" & cui != "C1333989") %>%
    dplyr::distinct()

  hereditary_cancers <- openxlsx::read.xlsx("data-raw/other/tumor_nos_umls.xlsx",sheet = 1) %>%
    dplyr::mutate(cui = stringr::str_trim(cui)) %>%
    dplyr::filter(name == "Hereditary_Cancer_Susceptibility_NOS") %>%
    dplyr::mutate(name = "Hereditary Cancer Susceptibility, NOS") %>%
    dplyr::bind_rows(syndromes_ding, susceptibility_ding) %>%
    dplyr::mutate(group = name) %>%
    dplyr::distinct()

  main_cancers_oncotree <- dplyr::bind_rows(onco_tree, hereditary_cancers) %>%
    dplyr::mutate(minor_type = dplyr::if_else(is.na(minor_type),FALSE,as.logical(minor_type))) %>%
    dplyr::select(cui, group, tissue, minor_type) %>%
    dplyr::filter(cui != "C0008479" | cui == "C0008479" & group == "Bone Cancer, NOS") %>%
    dplyr::filter(cui != "C0279672" | cui == "C0279672" & group == "Soft Tissue Cancer, NOS") %>%
    dplyr::filter(cui != "C0029463" | cui == "C0029463" & group == "Bone Cancer, NOS") %>%
    dplyr::filter(cui != "C0023467" | cui == "C0023467" & group == "Leukemia, NOS") %>%
    dplyr::filter(cui != "C1266144" | cui == "C1266144" & group == "Hereditary Cancer Susceptibility, NOS") %>%
    dplyr::filter(cui != "C0008497" | cui == "C0008497" & group == "Germ Cell Tumor, NOS") %>%
    dplyr::filter(cui != "C0025149" | cui == "C0025149" & group == "Hereditary Cancer Susceptibility, NOS") %>%
    dplyr::filter(cui != "C0031511" | cui == "C0031511" & group == "Hereditary Cancer Susceptibility, NOS") %>%
    dplyr::filter(cui != "C0035335" | cui == "C0035335" & group == "Eye Cancer, NOS") %>%
    dplyr::filter(cui != "C0334663" | cui == "C0334663" & group == "Soft Tissue Cancer, NOS") %>%
    dplyr::filter(cui != "C1266111" | cui == "C1266111" & group == "Soft Tissue Cancer, NOS") %>%
    dplyr::filter(cui != "C1332564" | cui == "C1332564" & group == "Bladder/Urinary Tract Cancer, NOS") %>%
    dplyr::filter(cui != "C0007134" | cui == "C0007134" & group == "Kidney Cancer, NOS") %>%
    dplyr::filter(cui != "C0334524" | cui == "C0334524" & tissue == "Vulva/Vagina") %>%
    dplyr::filter(cui != "C0346185" | cui == "C0346185" & tissue == "Ovary/Fallopian Tube") %>%
    dplyr::bind_rows(data.frame(cui = "C2608055", group = "Hereditary Cancer Susceptibility, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C0686584", tissue = "Myeloid", group = "Leukemia, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C0699791", tissue = "Esophagus/Stomach", group = "Esophageal/Stomach Cancer, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C0677055", tissue = "Vulva/Vagina", group = "Vulvar/Vaginal Cancer, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C1880169", tissue = "Soft Tissue", group = "Soft Tissue Cancer, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C0025500", tissue = "Pleura", group = "Pleural Cancer, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C0238122", tissue = "Ovary/Fallopian Tube", group = "Ovarian/Fallopian Tube Cancer, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C0023798", tissue = "Soft Tissue", group = "Soft Tissue Cancer, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C0279000", tissue = "Liver", group = "Liver Cancer, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C0348374", tissue = "CNS/Brain", group = "CNS/Brain Cancer, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C0346647", tissue = "Pancreas", group = "Pancreatic Cancer, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C0220624", tissue = "CNS/Brain", group = "CNS/Brain Cancer, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C0038874", tissue = "CNS/Brain", group = "CNS/Brain Cancer, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C0006118", tissue = "CNS/Brain", group = "CNS/Brain Cancer, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C0936223", tissue = "Prostate", group = "Prostate Cancer, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C0376358", tissue = "Prostate", group = "Prostate Cancer, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C0033578", tissue = "Prostate", group = "Prostate Cancer, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C0007102", tissue = "Colon/Rectum", group = "Colorectal Cancer, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C0009404", tissue = "Colon/Rectum", group = "Colorectal Cancer, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C0034885", tissue = "Colon/Rectum", group = "Colorectal Cancer, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C0009375", tissue = "Colon/Rectum", group = "Colorectal Cancer, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C1333976", tissue = "Liver", group = "Liver Cancer, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C0750952", tissue = "Biliary Tract", group = "Biliary Tract Cancer, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C0030297", tissue = "Pancreas", group = "Pancreatic Cancer, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C0015558", tissue = "Ovary/Fallopian Tube", group = "Ovarian/Fallopian Tube Cancer, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C0042995", tissue = "Vulva/Vagina", group = "Vulvar/Vaginal Cancer, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C0206765", tissue = "Soft Tissue", group = "Soft Tissue Cancer, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C0024623", tissue = "Esophagus/Stomach", group = "Esophageal/Stomach Cancer, NOS", stringsAsFactors = F)) %>%
    dplyr::bind_rows(data.frame(cui = "C0038356", tissue = "Esophagus/Stomach", group = "Esophageal/Stomach Cancer, NOS", stringsAsFactors = F)) %>%

    dplyr::mutate(minor_type = dplyr::if_else(is.na(minor_type),FALSE,as.logical(minor_type))) %>%
    dplyr::distinct() %>%
    dplyr::arrange(group)

  all_cancer_subtypes <- data.frame()
  i <- 1
  for(c in main_cancers_oncotree$cui){
    if(is.na(c)){
      next
    }
    gr <- unique(main_cancers_oncotree[!is.na(main_cancers_oncotree$cui) & main_cancers_oncotree$cui == c, ]$group)
    ti <- unique(main_cancers_oncotree[!is.na(main_cancers_oncotree$cui) & main_cancers_oncotree$cui == c, ]$tissue)
    minor_type <- unique(main_cancers_oncotree[!is.na(main_cancers_oncotree$cui) & main_cancers_oncotree$cui == c,]$minor_type)

    cancer_subtypes <- get_umls_children(c, umls_map = umls_map, umls_cui2name = cui_name_map) %>%
      dplyr::mutate(group = gr, tissue = ti, minor_type = minor_type)

    all_cancer_subtypes <- dplyr::bind_rows(all_cancer_subtypes, cancer_subtypes)
    cat(i, c,gr,nrow(all_cancer_subtypes),sep=" - ")
    cat("\n")
    i <- i + 1
  }

  all_cancer_subtypes_final <- all_cancer_subtypes %>%
    dplyr::filter(cui != "C0920349") %>% # Short Limb Dwarfism-Saddle Nose-Spinal Alterations-Metaphyseal Striation Syndrome
    dplyr::filter(cui != "C0009324") %>% # Ulcerative colitis
    dplyr::filter(cui != "C0566602") %>% # Primary sclerosing cholangitis
    dplyr::filter(cui != "C0010346") %>% # Crohn disease
    dplyr::filter(cui != "C0341332") %>% # Indeterminate colitis
    dplyr::filter(cui != "C0021390") %>% # Inflammatory bowel disease
    dplyr::filter(cui != "C0014527") %>% # Epidermolysis bullosa
    dplyr::filter(cui != "C0008029") %>% # Fibrous dysplasia of jaw
    dplyr::filter(cui != "C3854181") %>% # Nevus sebaceous
    dplyr::filter(cui != "C1839840") %>% # 46,XY sex reversal 8
    dplyr::filter(cui != "C0221026") %>% # X-linked agammaglobulinemia
    dplyr::filter(cui != "C1846545") %>% # Autoimmune lymphoproliferative syndrome type 2B
    dplyr::filter(cui != "C0238339") %>% # Hereditary pancreatitis
    dplyr::filter(group != "Peripheral Nervous System Cancer, NOS" |
                    (group == "Peripheral Nervous System Cancer, NOS" & stringr::str_detect(tolower(cui_name),"schwannom|neuroblastom|neurofibr|neuroendocrine|sheath|nerve|neuro|pheochromo|paraganglio"))) %>%
    dplyr::filter(!stringr::str_detect(tolower(cui_name),"diamond-blackfan|emochromatosis|hrombocytopenia| mononucleosis|performance status")) %>%
    dplyr::filter(!stringr::str_detect(tolower(cui_name),"sarcoidosis|myelofibrosis|neuromyelitis|myelodysplasia")) %>%
    dplyr::filter(!stringr::str_detect(tolower(cui_name),"familial pterygium|polycystic kidney| bullosa|nodular goiter")) %>%
    dplyr::filter(!stringr::str_detect(tolower(cui_name),"chalazion|neutropenia|nevus| lipidosis| myelitis| sclerosus")) %>%
    dplyr::filter(!stringr::str_detect(tolower(cui_name),"minimum residual disease|myasthenia gravis|alpha-1-antitrypsin")) %>%
    dplyr::filter(!stringr::str_detect(tolower(cui_name),"aplastic anemia|atrophic gastritis|leukoplakia|agammaglobulinemia|dyskeratosis")) %>%
    dplyr::filter(!stringr::str_detect(tolower(cui_name),"hypertension|erythroplakia|heel spur|grannuloma annulare| goiter")) %>%
    dplyr::filter(!stringr::str_detect(tolower(cui_name),"parapsoriasis | encephalitis | progression| heart disease| hyperplasia")) %>%
    dplyr::filter(!(stringr::str_detect(cui_name,"^Leio|Angioleiomyoma| Leio| leio") & stringr::str_detect(tissue,"Esophagus"))) %>%
    dplyr::filter(!(stringr::str_detect(cui_name,"(O|o)varian|(T|t)esticular|seminoma") & stringr::str_detect(tissue,"CNS/Brain"))) %>%
    dplyr::filter(!(stringr::str_detect(cui_name,"Wilms") & stringr::str_detect(tissue,"CNS|Bladder"))) %>%
    dplyr::filter(!(stringr::str_detect(cui_name,"(K|k)idney") & stringr::str_detect(tissue,"Bladder"))) %>%
    #dplyr::filter(!(stringr::str_detect(cui_name,"(E|e)sophag") & stringr::str_detect(tissue,"Head"))) %>%
    dplyr::filter(!(stringr::str_detect(cui_name,"^Lung| lung ") & stringr::str_detect(tissue,"Skin"))) %>%
    dplyr::filter(!(stringr::str_detect(tolower(cui_name), "lympho|lympha") & !stringr::str_detect(tolower(cui_name),"skin|papulosis") & (!is.na(tissue) & tissue == "Skin"))) %>%
    dplyr::filter(!(stringr::str_detect(tolower(cui_name), "ewing sarcom") & !stringr::str_detect(tolower(cui_name),"neuro") & (!is.na(tissue) & tissue == "CNS/Brain"))) %>%
    dplyr::filter(!(stringr::str_detect(tolower(cui_name), "sarcoma") & (!is.na(tissue) & tissue == "Myeloid" | tissue == "Lymphoid"))) %>%
    dplyr::filter(!(stringr::str_detect(tolower(cui_name), "renal cell") & (!is.na(tissue) & tissue == "Bladder/Urinary Tract"))) %>%
    
    dplyr::arrange(group) %>%
    dplyr::mutate(group = stringr::str_replace_all(group," |/|, ","_")) %>%
    dplyr::mutate(group = stringr::str_replace_all(group,"_and_","_And_")) %>%
    dplyr::mutate(tissue = dplyr::if_else(
      tissue == "Other","Other/Unknown",
      as.character(tissue))) %>%
    dplyr::mutate(primary_site = tissue) %>%
    dplyr::mutate(group =
                    dplyr::if_else(stringr::str_detect(tolower(cui_name),
                                                       "thyroid") &
                                     !is.na(primary_site) &
                                     primary_site == "Head and Neck",
                                   "Thyroid_Cancer_NOS",
                                   as.character(group))) %>%
    dplyr::mutate(
      tissue =
        dplyr::if_else(stringr::str_detect(tolower(cui_name),
                                           "thyroid") &
                         !is.na(primary_site) &
                         primary_site == "Head and Neck",
                       "Thyroid",
                       as.character(tissue))) %>%
    dplyr::mutate(
      primary_site =
        dplyr::if_else(stringr::str_detect(tolower(cui_name),
                                           "thyroid") &
                         primary_site == "Head and Neck" &
                         !is.na(primary_site),
                       "Thyroid",
                       as.character(primary_site))) %>%

    dplyr::mutate(
      group =
        dplyr::if_else(stringr::str_detect(tolower(cui_name),
                                           "barrett") &
                         primary_site == "Head and Neck" &
                         !is.na(primary_site),
                       "Esophageal_Stomach_Cancer_NOS",
                       as.character(group))) %>%
    dplyr::mutate(
      tissue =
        dplyr::if_else(stringr::str_detect(tolower(cui_name),
                                           "barrett") &
                         primary_site == "Head and Neck"&
                         !is.na(primary_site),
                       "Esophagus/Stomach",
                       as.character(tissue))) %>%
    dplyr::mutate(
      primary_site =
        dplyr::if_else(stringr::str_detect(tolower(cui_name),
                                           "barrett") &
                         primary_site == "Head and Neck" &
                         !is.na(primary_site),
                       "Esophagus/Stomach",
                       as.character(primary_site))) %>%
    dplyr::distinct()

  # remove overlapping entries in hereditary cancer syndrome group and cancer susceptibility group
  hereditary_condition_entries <-
    dplyr::filter(all_cancer_subtypes_final, group == "Hereditary_Cancer_Susceptibility_NOS")
  hereditary_syndrome_entries <-
    dplyr::filter(all_cancer_subtypes_final, group == "Hereditary_Cancer_Syndrome_NOS")
  hereditary_condition_entries <-
    dplyr::anti_join(hereditary_condition_entries,
                     dplyr::select(hereditary_syndrome_entries, cui))
  sporadic_cancer_entries <-
    dplyr::filter(all_cancer_subtypes_final, group != "Hereditary_Cancer_Susceptibility_NOS" &
                    group != "Hereditary_Cancer_Syndrome_NOS")

  final_tree <- dplyr::bind_rows(hereditary_condition_entries,
                                 hereditary_syndrome_entries,
                                 sporadic_cancer_entries) %>%
    dplyr::filter(!is.na(cui))

  final_tree_slim <- final_tree %>% dplyr::filter(minor_type == F)
  l <- list("full" = final_tree, "slim" = final_tree_slim)
  return(l)

}

get_umls_children <- function(c, umls_map = NULL, umls_cui2name = NULL){
  children_level0 <- data.frame("cui" = c, stringsAsFactors = F) %>%
    dplyr::inner_join(umls_cui2name, by = c("cui"))
  cui_children <- dplyr::filter(umls_map$relation, cui == c) %>%
    dplyr::select(cui_child)

  ## LEVEL 1
  children_level1 <- data.frame()
  if(nrow(cui_children) > 0){
    children_level1 <- data.frame("cui" = cui_children$cui_child,
                                  "level" = 1, stringsAsFactors = F) %>%
      dplyr::inner_join(umls_cui2name, by = c("cui")) %>%
      dplyr::select(-level)
  }

  ## LEVEL 2
  children_level2 <- data.frame()
  for(c1 in unique(children_level1$cui)){
    children <- dplyr::filter(umls_map$relation, cui == c1)
    if(nrow(children) > 0){
      children <- children %>%
        dplyr::select(cui_child) %>%
        dplyr::rename(cui = cui_child)
      children_level2 <- dplyr::bind_rows(children_level2, children) %>%
        dplyr::mutate(level = 2)
    }
  }
  if(nrow(children_level2) > 0){
    children_level2 <- children_level2 %>%
      dplyr::select(-level) %>%
      dplyr::distinct() %>%
      dplyr::inner_join(umls_cui2name, by = c("cui"))
  }

  ## LEVEL 3
  children_level3 <- data.frame()
  for(c2 in unique(children_level2$cui)){
    children <- dplyr::filter(umls_map$relation, cui == c2)
    if(nrow(children) > 0){
      children <- children %>%
        dplyr::select(cui_child) %>%
        dplyr::rename(cui = cui_child)
      children_level3 <- dplyr::bind_rows(children_level3, children) %>%
        dplyr::mutate(level = 3)
    }
  }
  if(nrow(children_level3) > 0){
    children_level3 <- children_level3 %>%
      dplyr::select(-level) %>%
      dplyr::distinct() %>%
      dplyr::inner_join(umls_cui2name, by = c("cui"))
  }


  ## LEVEL 4
  children_level4 <- data.frame()
  for(c3 in unique(children_level3$cui)){
    children <- dplyr::filter(umls_map$relation, cui == c3)
    if(nrow(children) > 0){
      children <- children %>%
        dplyr::select(cui_child) %>%
        dplyr::rename(cui = cui_child)
      children_level4 <- dplyr::bind_rows(children_level4, children) %>%
        dplyr::mutate(level = 4)
    }
  }
  if(nrow(children_level4) > 0){
    children_level4 <- children_level4 %>%
      dplyr::select(-level) %>%
      dplyr::distinct() %>%
      dplyr::inner_join(umls_cui2name, by = c("cui"))
  }

  ## LEVEL 5
  children_level5 <- data.frame()
  for(c4 in unique(children_level4$cui)){
    children <- dplyr::filter(umls_map$relation, cui == c4)
    if(nrow(children) > 0){
      children <- children %>%
        dplyr::select(cui_child) %>%
        dplyr::rename(cui = cui_child)
      children_level5 <- dplyr::bind_rows(children_level5, children) %>%
        dplyr::mutate(level = 5)
    }
  }
  if(nrow(children_level5) > 0){
    children_level5 <- children_level5 %>%
      dplyr::select(-level) %>%
      dplyr::distinct() %>%
      dplyr::inner_join(umls_cui2name, by = c("cui"))
  }


  ## LEVEL 5
  children_level6 <- data.frame()
  for(c5 in unique(children_level5$cui)){
    children <- dplyr::filter(umls_map$relation, cui == c5)
    if(nrow(children) > 0){
      children <- children %>%
        dplyr::select(cui_child) %>%
        dplyr::rename(cui = cui_child)
      children_level6 <- dplyr::bind_rows(children_level6, children) %>%
        dplyr::mutate(level = 6)
    }
  }
  if(nrow(children_level6) > 0){
    children_level6 <- children_level6 %>%
      dplyr::select(-level) %>%
      dplyr::distinct() %>%
      dplyr::inner_join(umls_cui2name, by = c("cui"))
  }

  all_children <- children_level0 %>%
    dplyr::bind_rows(children_level1) %>%
    dplyr::bind_rows(children_level2) %>%
    dplyr::bind_rows(children_level3) %>%
    dplyr::bind_rows(children_level4) %>%
    dplyr::bind_rows(children_level5) %>%
    dplyr::bind_rows(children_level6) %>%
    dplyr::distinct()

  return(all_children)
}


disease_ontology_release <- '2021-08-17v2'
efo_release <- 'v3.33.0'
oncotree_stable_release <- 'oncotree_2020_10_01'

## Get UMLS / DiseaseOntology / EFO mappings
umls_map <- map_umls(update = T,
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

efo2xref <- efo_map$efo2xref
efo2name <- efo_map$efo2name

## Use OncoTree as starting point for phenotype cross-map
onco_pheno_map <- onco_pheno_mapping(umls_map = umls_map)

## Cross-reference OncoTree with EFO and DO
for(m in c('slim','full')){
  onco_pheno_map[[m]] <- onco_pheno_map[[m]] %>%
    dplyr::left_join(dplyr::select(efo_map$efo2xref,cui,efo_id,efo_name),
                     by = "cui") %>%
    dplyr::left_join(do_map, by = "cui") %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      efo_id =
        dplyr::if_else(
          group == "Colorectal_Cancer_NOS" &
            stringr::str_detect(tolower(cui_name),
                                "metasta"),
          "EFO:1001480",as.character(efo_id))) %>%
    dplyr::mutate(
      efo_name =
        dplyr::if_else(
          group == "Colorectal_Cancer_NOS" &
            stringr::str_detect(tolower(cui_name),
                                "metasta"),
          "metastatic colorectal cancer",
          as.character(efo_name))) %>%
    dplyr::mutate(disease_ontology_release = disease_ontology_release,
                  efo_release = efo_release,
                  oncotree_release = oncotree_stable_release)
}

version_date <- format(Sys.Date(), format="%Y%m%d")

opm_slim <- onco_pheno_map$slim
opm_full <- onco_pheno_map$full

usethis::use_data(umls_map, overwrite = T)
usethis::use_data(do_map, overwrite = T)
usethis::use_data(efo2xref, overwrite = T)
usethis::use_data(efo2name, overwrite = T)
usethis::use_data(opm_slim, overwrite = T)
usethis::use_data(opm_full, overwrite = T)



# write.table(oncotree_umls_map$slim,file=paste0("phenotype_ontology_cancer.",version_date,".tsv"),sep="\t",col.names = T,row.names = F,quote = F)
# system(paste0('ln -sF phenotype_ontology_cancer.',version_date,'.tsv phenotype_ontology_cancer.tsv'))
# saveRDS(oncotree_umls_map$slim,file="phenotype_ontology_cancer.rds")

# write.table(onco_pheno_map$slim,
#             file=paste0("output/onco_pheno_map_",version_date,".tsv"),
#             sep="\t",col.names = T,row.names = F,quote = F)
# system(paste0('ln -sF output/onco_pheno_map_',version_date,'.tsv onco_pheno_map.tsv'))
# write.table(onco_pheno_map$full,
#             file=paste0("output/onco_pheno_map_full_",version_date,".tsv"),
#             sep="\t",col.names = T,row.names = F,quote = F)
# system(paste0('ln -sF output/onco_pheno_map_full_', version_date,'.tsv onco_pheno_map_full.tsv'))
# saveRDS(onco_pheno_map$slim,file="output/onco_pheno_map.rds")
# saveRDS(onco_pheno_map$full,file="output/onco_pheno_map_full.rds")
