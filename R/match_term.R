#' @title
#' Match a query term with oncology-relevant ontology terms 
#'
#' @description
#' <Description of how the function works>
#'
#' @param query_term query phenotype term
#' @param match_type how to perform ontology-based matching
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
#' 
#' @export
#'

match_term <- function(query_term = NA,
                       match_type = "exact",
                       cache_dir = NA, 
                      force_download = FALSE,
                      site = NA) {
    
    ## TODO: check that match type is any of allowed values
    ## (to be determined)
    
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
    
    ## TODO: find matches between dat$records and query_term
    
    ## if(match_type == "exact"){
    ##
    ## }
    
    ## and so on
    
    #return(dat)
}
