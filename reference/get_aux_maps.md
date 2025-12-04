# Get auxiliary phenotype ontology maps

Downloads and returns multiple phenotype dictionaries used for
cross-referencing of terms, including UMLS, EFO, DO, and ICD10. The
dataset comes as a `list` object, with two elements:

- `metadata` - a data frame with metadata regarding annotation resources
  used

- `records` - a list object with multiple lists/data.frames

## Usage

``` r
get_aux_maps(cache_dir = NA, force_download = FALSE)
```

## Arguments

- cache_dir:

  Local directory for data download

- force_download:

  Logical indicating if local cache should be overwritten (set to TRUE
  to re-download if file exists in cache)

## Value

**metadata** - A data frame with 4 rows and 6 columns:

- *source* - gene annotation source

- *annotation_data* - type of annotations used

- *url* - URL of annotation resource

- *citation* - publication to cite for annotation source (citation;
  PMID)

- *version* - version used

- *abbreviation* - abbreviation used in column names of records

**records** - A list with the following elements:

- *umls* - List with UMLS data dictionaries

- *do* - Data frame with DO identifiers/names

- *efo* - List with EFO data dictionaries

- *icd10* - Data frame with CUI to ICD10 cross-reference

## Examples

``` r
if (FALSE) { # \dontrun{
library(phenOncoX)
oncology_terms <- get_aux_maps(cache_dir = tempdir())
} # }
```
