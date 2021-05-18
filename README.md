### oncoPhenoMap - crossmapped phenotype ontology terms for the oncology domain

#### Overview

This R data package attempts to provide a crossmapped set of phenotype ontology terms attributed to cancer phenotypes. The mapping is semi-manually created, using [OncoTree](http://oncotree.mskcc.org/#/home) as the starting point for a list of [UMLS/MedGen](https://www.ncbi.nlm.nih.gov/medgen/) phenotype terms per cancer subtype/primary site. We make cross-mappings between phenotype terms from the [Experimental Factor Ontology (EFO)](https://github.com/EBISPOT/efo) as well as the [Disease Ontology (DO)](https://disease-ontology.org/).

Currently, the following versions are used to create the mapping:

 - OncoTree (2020_10_01)
 - Experimental Factor Ontology v30.0 (2021-05-17)
 - Disease Ontology (April 2021 release)


#### Usage

`devtools::install_github('sigven/oncoPhenoMap')`

- Phenotype mapping records are found in the exported _opm_slim_ object:

`head(oncoPhenoMap::opm_slim)`


#### Contact

sigven AT ifi.uio.no
