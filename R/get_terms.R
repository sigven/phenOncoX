#' @title
#' Get oncology-relevant phenotype terms 
#'
#' @description
#' Downloads and returns a dataset with an expanded set of oncology-relevant 
#' phenotype terms (UMLS), using OncoTree as a starting point. The expansion
#' has been conducted by tracking the MeSH/UMLS tree structure with 
#' phenotype terms. Each record comes with cross-references to Disease 
#' Ontology (DO), Experimental Factor Ontology (EFO), and ICD10 wherever 
#' available. The dataset comes as a `list` object, with two elements:
#'
#' * `metadata` - a data frame with metadata regarding annotation 
#'                resources used
#' * `records` - a data frame with phenotype terms
#'
#' @param cache_dir Local directory for data download
#' @param force_download Logical indicating if local cache should be overwritten
#' (set to TRUE to re-download if file exists in cache)
#' @param site limit phenotype terms to those relevant for a given
#' primary tumor type/site. Possible values: 
#' 'Adrenal Gland', 'Ampulla of Vater', 'Biliary Tract', 'Bladder/Urinary Tract', 
#' 'Bone', 'Breast', 'Cervix', 'CNS/Brain', 'Colon/Rectum', 'Esophagus/Stomach', 
#' 'Eye', 'Head and Neck', 'Kidney', 'Liver', 'Lung', 'Lymphoid', 'Myeloid', 
#' 'Other/Unknown', 'Ovary/Fallopian Tube', 'Pancreas', 'Penis', 
#' 'Peripheral Nervous System', 'Peritoneum', 'Pleura', 'Prostate', 
#' 'Skin', 'Soft Tissue', 'Testis', 'Thymus', 'Thyroid', 'Uterus', 
#' 'Vulva/Vagina'
#' @param ignore_minor_type logical indicating if minor tumor types
#' should be excluded or not
#' 
#' @returns
#' \strong{metadata} - A data frame with 5 rows and 6 columns:
#' \itemize{
#'   \item \emph{source} - gene annotation source
#'   \item \emph{annotation_data} - type of annotations used
#'   \item \emph{url} - URL of annotation resource
#'   \item \emph{citation} - publication to cite for annotation source
#'   (citation; PMID)
#'   \item \emph{version} - version used
#'   \item \emph{abbreviation} - abbreviation used in column names of records
#'  }
#'
#' \strong{records} - A data frame the following columns:
#' \itemize{
#'   \item \emph{primary_site} - Primary tumor type/site
#'   \item \emph{ot_main_type} - Main tumor type (OncoTree)
#'   \item \emph{ot_name} - Phenotype name (OncoTree)
#'   \item \emph{ot_level} - Tree level (OncoTree)
#'   \item \emph{ot_code} - Phenotype code (OncoTree)
#'   \item \emph{ot_code_path} - Tree code path (OncoTree)
#'   \item \emph{cui} - Concept unique identifier (UMLS)
#'   \item \emph{cui_name} - Phenotype name (UMLS)
#'   \item \emph{efo_id} - EFO identifier
#'   \item \emph{efo_name} - EFO name
#'   \item \emph{do_id} - DO identifier
#'   \item \emph{do_name} - DO name
#'   \item \emph{do_cancer_slim} - DO identifier part of DO cancer slim 
#'   dictionary (TRUE/FALSE)
#'   \item \emph{minor_type} - logical indicating whether the term is part
#'   of a minor tumor type/site 
#'   \item \emph{icd10_code} - ICD10 identifier
#'   \item \emph{source} - term indicating the dataset source

#' }
#'
#' @examples
#' 
#' \dontrun{
#' library(phenOncoX)
#' oncology_terms <- get_terms(cache_dir = tempdir())
#' }
#' 
#' @export
#'

get_terms <- function(cache_dir = NA, 
                      force_download = FALSE,
                      site = NA,
                      ignore_minor_type = FALSE) {
    dat <- get_pox_data(
        cache_dir = cache_dir,
        force_download = force_download,
        db = "oncotree_expanded"
    )
    
    if (ignore_minor_type == TRUE) {
        dat$records <- dat$records |>
            dplyr::filter(
                .data$minor_type == FALSE
            )
    }
    
    if (NROW(dat$records) > 0 & !is.na(site)) {
        
        allowed_psites <- 
            unique(dat$records$primary_site)
        allowed_psites <- 
            allowed_psites[
                !is.na(unique(dat$records$primary_site))
            ]
        
        if (!(site %in% allowed_psites)) {
            lgr::lgr$error(paste0(
                "Argument site = '",
                site, "' is not allowed - possible values: '",
                paste(allowed_psites, collapse = "', '"),"'"
            ))
        }else {
            lgr::lgr$info(paste0(
                "Limiting phenotype terms to '", site, "'"
            ))
            dat$records <- dat$records |>
                dplyr::filter(!is.na(.data$primary_site)) |>
                dplyr::filter(
                    .data$primary_site == site
                )
        }
    }
    
    return(dat)
}
