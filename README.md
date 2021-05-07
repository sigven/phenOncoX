### oncoPhenoMap - crossmapped phenotype ontologies for the oncology domain

#### Overview

This R package attempts to provide a crossmapped set of phenotype ontology terms attributed to cancer phenotypes. The mapping is semi-manually created, using [OncoTree](http://oncotree.mskcc.org/#/home) as the starting point for a list of [UMLS/MedGen](https://www.ncbi.nlm.nih.gov/medgen/) phenotype terms per cancer subtype/primary site. We make cross-mappings between phenotype terms from the [Experimental Factor Ontology (EFO)](https://github.com/EBISPOT/efo) as well as the [Disease Ontology (DO)](https://disease-ontology.org/).

Currently, the following versions are used to create the mapping:

 - OncoTree (2020_10_01)
 - Experimental Factor Ontology v29.1 (2021-04-20)
 - Disease Ontology (April 2021 release)


#### Usage

`devtools::install_github('sigven/oncoPhenoMap')`
`head(oncoPhenoMap::onco_pheno_slim)`


#### Contact

sigven AT ifi.uio.no
