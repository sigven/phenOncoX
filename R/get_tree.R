#' @title
#' Get core phenotype terms from OncoTree
#'
#' @description
#' Downloads and returns a dataset with curated phenotype terms (UMLS) from 
#' OncoTree. The dataset comes as a `list` object, with two elements:
#'
#' * `metadata` - a data frame with metadata regarding annotation 
#'                resources used
#' * `records` - a data frame with phenotype terms
#'
#' @param cache_dir Local directory for data download
#' @param force_download Logical indicating if local cache should be overwritten
#' (set to TRUE to re-download if file exists in cache)
#' @param site limit OncoTree terms to those relevant for a given
#' primary tumor type/site. Possible values: 
#' 'Adrenal Gland', 'Ampulla of Vater', 'Biliary Tract', 'Bladder/Urinary Tract', 
#' 'Bone', 'Breast', 'Cervix', 'CNS/Brain', 'Colon/Rectum', 'Esophagus/Stomach', 
#' 'Eye', 'Head and Neck', 'Kidney', 'Liver', 'Lung', 'Lymphoid', 'Myeloid', 
#' 'Other/Unknown', 'Ovary/Fallopian Tube', 'Pancreas', 'Penis', 
#' 'Peripheral Nervous System', 'Peritoneum', 'Pleura', 'Prostate', 
#' 'Skin', 'Soft Tissue', 'Testis', 'Thymus', 'Thyroid', 'Uterus', 
#' 'Vulva/Vagina'
#' @param max_tree_depth consider only terms up to a given depth in
#' the OncoTree (integer from 1-6, default: 6)
#' @param ignore_minor_type logical indicating if minor tumor types
#' should be excluded or not
#' 
#' @returns
#' \strong{metadata} - A data frame with 1 row and 6 columns:
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
#' \strong{records} - A data frame with the following columns:
#' \itemize{
#'   \item \emph{primary_site} - Primary tumor type/site
#'   \item \emph{ot_main_type} - Main tumor type (OncoTree)
#'   \item \emph{ot_name} - Phenotype name (OncoTree)
#'   \item \emph{ot_level} - Tree level (OncoTree)
#'   \item \emph{ot_code} - Phenotype code (OncoTree)
#'   \item \emph{ot_code_path} - Tree code path (OncoTree)
#'   \item \emph{cui} - Concept unique identifier (UMLS)
#'   \item \emph{cui_name} - Phenotype name (UMLS)
#'   \item \emph{minor_type} - logical indicating whether the term is part
#'   of a minor tumor type/site 
#'   \item \emph{source} - term indicating the dataset source

#' }
#'
#' @examples
#' 
#' \dontrun{
#' library(phenOncoX)
#' oncology_terms <- get_tree(cache_dir = tempdir())
#' }
#' 
#' @export
#'

get_tree <- function(cache_dir = NA, 
                      force_download = FALSE,
                      site = NA,
                      max_tree_depth = 6,
                      ignore_minor_type = FALSE) {
    
    lgr::lgr$appenders$console$set_layout(
        lgr::LayoutFormat$new(timestamp_fmt = "%Y-%m-%d %T")
    )
    
    dat <- get_pox_data(
        cache_dir = cache_dir,
        force_download = force_download,
        db = "oncotree_core"
    )
    
    if (!is.integer(max_tree_depth) & !is.numeric(max_tree_depth)) {
        lgr::lgr$error(paste0(
            "Argument max_tree_depth = '",
            max_tree_depth, "' is not an integer"
        ))
        return(dat)
    }else{
        if (max_tree_depth < 1 | max_tree_depth > 6) {
            lgr::lgr$error(paste0(
                "Argument max_tree_depth = '",
                max_tree_depth, "' must be within the range of 1 to 6"
            ))
            return(dat)
        }else{
            if(max_tree_depth < 6){
                lgr::lgr$info(paste0(
                    "Limiting OncoTree terms to those with a ",
                    "maximum tree depth of = '",
                    max_tree_depth, "'"
                ))
            }
        }
    }
    dat$records <- dat$records |>
        dplyr::filter(
            .data$ot_level <= max_tree_depth
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
                "Argument primary_site = '",
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
