### oncoPhenoMap - crossmapped phenotype ontology terms for the oncology domain

#### Overview

An ontological definition of disease enables each type of disease to be singularly classified in a formalized structure. By intention, the use of disease ontology terms should facilitate a cross-link of information between separate disease-related knowledge resources for a given domain. However, multiple disease ontology frameworks have been developed for human disease (i.e. [OncoTree](http://oncotree.mskcc.org/#/home), [Experimental Factor Ontology (EFO)](https://github.com/EBISPOT/efo), [Disease Ontology (DO)](https://github.com/DiseaseOntology/HumanDiseaseOntology), [UMLS](https://www.ncbi.nlm.nih.gov/medgen/), [ICD-10)](https://www.who.int/standards/classifications/classification-of-diseases), and they are used to different extents across knowledge resources in the oncology domain, most importantly:

-   gene-disease associations
-   drug-disease indications
-   disease-mutation associations (prognostic, drug sensitivity/resistance etc.)

In order to integrate such knowledge resources, there is henceforth a need to cross-link or map the entries across disease ontologies to the extent it is possible.

**oncoPhenoMap** is an R data package that attempts to address this challenge. In short, **oncoPhenoMap** provides a global cross-mapped set of phenotype ontology terms attributed to cancer phenotypes.

The mapping established within **oncoPhenoMap** is semi-manually curated, using [OncoTree](http://oncotree.mskcc.org/#/home) as the starting point for a list of [UMLS](https://www.ncbi.nlm.nih.gov/medgen/) phenotype terms per cancer subtype/primary site. Next, **oncoPhenoMap** appends a number of phenotypes attributed to heritable cancer conditions. Furthermore, each cancer subtype entry in *OncoTree* is expanded with additional subtypes that are found in the *UMLS* child-parent hierarchy of disease terms.

For each entry in the final list of phenotype terms, we make cross-mappings with phenotype terms from [EFO](https://github.com/EBISPOT/efo), [DO](https://disease-ontology.org/), and the [ICD10 classification](https://www.who.int/standards/classifications/classification-of-diseases).

Currently, the following versions are used to create the mapping:

-   OncoTree (2020_10_01)
-   Experimental Factor Ontology v34.0 (2021-09-16)
-   Disease Ontology (August 2021 release)

**IMPORTANT NOTE**: The mapping established by **oncoPhenoMap** attempts to be comprehensive, but we acknowledge that the presence of missing or erroneous cross-references might still occur.

#### Usage

`devtools::install_github('sigven/oncoPhenoMap')`

-   Key phenotype mapping records are found in three exported data frames:

    1.  *oncoPhenoMap::oncotree_basic* - all terms in OncoTree (UMLS-mapped), and curated/error-corrected as much as possible
    2.  *oncoPhenoMap::oncotree_expanded_full* - all terms in OncoTree (UMLS-mapped), with the additions of children in the UMLS hierarchy of disease terms, and terms related to heritable cancer conditions (complete set)
    3.  *oncoPhenoMap::oncotree_expanded_slim* - all terms in OncoTree (UMLS-mapped), with the additions of children in the UMLS hierarchy of disease terms, and terms related to heritable cancer conditions (slim set - ignoring minor/rare tumor types)

    All records in *1.*, *2.*, and *3.* contain cross-references to corresponding **EFO**, **DO** and **ICD10** terms or codes.

    An additional list structured coined *oncoPhenoMap::auxiliary_maps* provides mapping data across different disease ontologies

#### Contact

sigven AT ifi.uio.no
