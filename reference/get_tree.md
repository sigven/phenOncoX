# Get core phenotype terms from OncoTree

Downloads and returns a dataset with curated phenotype terms (UMLS) from
OncoTree. The dataset comes as a `list` object, with two elements:

- `metadata` - a data frame with metadata regarding annotation resources
  used

- `records` - a data frame with phenotype terms

## Usage

``` r
get_tree(
  cache_dir = NA,
  force_download = FALSE,
  site = NA,
  max_tree_depth = 6,
  ignore_minor_type = FALSE
)
```

## Arguments

- cache_dir:

  Local directory for data download

- force_download:

  Logical indicating if local cache should be overwritten (set to TRUE
  to re-download if file exists in cache)

- site:

  limit OncoTree terms to those relevant for a given primary tumor
  type/site. Possible values: 'Adrenal Gland', 'Ampulla of Vater',
  'Biliary Tract', 'Bladder/Urinary Tract', 'Bone', 'Breast', 'Cervix',
  'CNS/Brain', 'Colon/Rectum', 'Esophagus/Stomach', 'Eye', 'Head and
  Neck', 'Kidney', 'Liver', 'Lung', 'Lymphoid', 'Myeloid',
  'Other/Unknown', 'Ovary/Fallopian Tube', 'Pancreas', 'Penis',
  'Peripheral Nervous System', 'Peritoneum', 'Pleura', 'Prostate',
  'Skin', 'Soft Tissue', 'Testis', 'Thymus', 'Thyroid', 'Uterus',
  'Vulva/Vagina'

- max_tree_depth:

  consider only terms up to a given depth in the OncoTree (integer from
  1-6, default: 6)

- ignore_minor_type:

  logical indicating if minor tumor types should be excluded or not

## Value

**metadata** - A data frame with 1 row and 6 columns:

- *source* - gene annotation source

- *annotation_data* - type of annotations used

- *url* - URL of annotation resource

- *citation* - publication to cite for annotation source (citation;
  PMID)

- *version* - version used

- *abbreviation* - abbreviation used in column names of records

**records** - A data frame with the following columns:

- *primary_site* - Primary tumor type/site

- *ot_main_type* - Main tumor type (OncoTree)

- *ot_name* - Phenotype name (OncoTree)

- *ot_level* - Tree level (OncoTree)

- *ot_code* - Phenotype code (OncoTree)

- *ot_code_path* - Tree code path (OncoTree)

- *cui* - Concept unique identifier (UMLS)

- *cui_name* - Phenotype name (UMLS)

- *minor_type* - logical indicating whether the term is part of a minor
  tumor type/site

- *source* - term indicating the dataset source

## Examples

``` r
if (FALSE) { # \dontrun{
library(phenOncoX)
oncology_terms <- get_tree(cache_dir = tempdir())
} # }
```
