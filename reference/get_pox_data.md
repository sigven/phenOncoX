# Function that retrieves phenOncoX data from Google Drive/local cache

Function that retrieves phenOncoX data from Google Drive/local cache

## Usage

``` r
get_pox_data(cache_dir = NA, force_download = FALSE, db = "oncotree_core")
```

## Arguments

- cache_dir:

  Local directory for data download

- force_download:

  Logical indicating if local cache should be overwritten (set to TRUE
  to re-download if file exists in cache)

- db:

  type of dataset to be retrieved

## Value

pre-processed ontology terms records/metadata or auxiliary maps
