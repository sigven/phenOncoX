#' @title
#' Get auxiliary phenotype ontology maps 
#'
#' @description
#' Downloads and returns multiple phenotype dictionaries used for
#' cross-referencing of terms, including UMLS, EFO, DO, and ICD10. 
#' The dataset comes as a `list` object, with two elements:
#'
#' * `metadata` - a data frame with metadata regarding annotation 
#'                resources used
#' * `records` - a list object with multiple lists/data.frames
#'
#' @param cache_dir Local directory for data download
#' @param force_download Logical indicating if local cache should be overwritten
#' (set to TRUE to re-download if file exists in cache)
#' 
#' @returns
#' \strong{metadata} - A data frame with 4 rows and 6 columns:
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
#' \strong{records} - A list with the following elements:
#' \itemize{
#'   \item \emph{umls} - List with UMLS data dictionaries
#'   \item \emph{do} - Data frame with DO identifiers/names
#'   \item \emph{efo} - List with EFO data dictionaries
#'   \item \emph{icd10} - Data frame with CUI to ICD10 cross-reference
#' }
#'
#' @examples
#' 
#' \dontrun{
#' library(geneOncoX)
#' oncology_terms <- get_aux_maps(cache_dir = tempdir())
#' }
#' 
#' @export
#'

get_aux_maps <- function(cache_dir = NA, 
                      force_download = FALSE,
                      primary_site = NA,
                      ignore_minor_type = FALSE) {
    dat <- get_pox_data(
        cache_dir = cache_dir,
        force_download = force_download,
        db = "auxiliary_maps"
    )
    
    return(dat)
}
