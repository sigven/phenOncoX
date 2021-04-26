### oncoPhenoMap - crossmapped phenotype ontologies for the oncology domain

#### Overview

This project attempts to provide a crossmapped set of phenotype ontology terms that are linked to cancer phenotypes. The mapping is semi-manually created, using [OncoTree](http://oncotree.mskcc.org/#/home) as the starting point for a list of [UMLS/MedGen](https://www.ncbi.nlm.nih.gov/medgen/) phenotype terms per cancer subtype/primary site. We make cross-mappings between phenotype terms from the [Experimental Factor Ontology (EFO)](https://github.com/EBISPOT/efo) as well as the [Disease Ontology (DO)](https://disease-ontology.org/).

Currently, the following versions are used to create the mapping:

 - OncoTree (2020_10_01)
 - Experimental Factor Ontology v29.1 (2021-04-20)
 - Disease Ontology (March 2021 release)

 The key output is found as tab-separated values (TSV) files in the output folder:

 1. `output/onco_pheno_map_full_20210426.tsv`
 2. `output/onco_pheno_map_20210426.tsv`


#### Contact

sigven AT ifi.uio.no
