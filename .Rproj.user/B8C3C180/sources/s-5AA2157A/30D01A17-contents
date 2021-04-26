source(file.path(here::here(),'R','ontology_db_utils.R'))

disease_ontology_release <- '2021-03-29'
efo_release <- 'v3.29.1'
oncotree_stable_release <- 'oncotree_2020_10_01'

## Get UMLS / DiseaseOntology / EFO mappings
umls_map <- map_umls(update = F)
saveRDS(umls_map,file="output/umls_map.rds")
do_map <- map_disease_ontology(skip_non_cui_mapped = T, 
                               release = disease_ontology_release, 
                               umls_map = umls_map)
saveRDS(do_map,file="output/do_map.rds")
efo_map <- map_efo(umls_map = umls_map, 
                   efo_release = efo_release, 
                   update = T)
saveRDS(efo_map$efo2xref, file="output/efo_map.rds")
saveRDS(efo_map$efo2name, file="output/efo2name.rds")

## Use OncoTree as starting point for phenotypes
onco_pheno_map <- map_oncotree_ontology(umls_map = umls_map)

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
# write.table(oncotree_umls_map$slim,file=paste0("phenotype_ontology_cancer.",version_date,".tsv"),sep="\t",col.names = T,row.names = F,quote = F)
# system(paste0('ln -sF phenotype_ontology_cancer.',version_date,'.tsv phenotype_ontology_cancer.tsv'))
# saveRDS(oncotree_umls_map$slim,file="phenotype_ontology_cancer.rds")

write.table(onco_pheno_map$slim, 
            file=paste0("output/onco_pheno_map_",version_date,".tsv"), 
            sep="\t",col.names = T,row.names = F,quote = F)
system(paste0('ln -sF output/onco_pheno_map_',version_date,'.tsv onco_pheno_map.tsv'))
write.table(onco_pheno_map$full, 
            file=paste0("output/onco_pheno_map_full_",version_date,".tsv"), 
            sep="\t",col.names = T,row.names = F,quote = F)
system(paste0('ln -sF output/onco_pheno_map_full', version_date,'.tsv onco_pheno_map_full.tsv'))
saveRDS(onco_pheno_map$slim,file="output/onco_pheno_map.rds")
saveRDS(onco_pheno_map$full,file="output/onco_pheno_map_full.rds")
