# NUCLEON &mdash; Nuclide and Element Ontology

The Nuclide and Element Ontology (NUCLEON) is a machine-processable knowledge resource that represents 118 chemical elements and 4,105 nuclides, including 716 metastable isomers, with computable nuclear properties.

The ontology declares each element and nuclide as both a class and an individual: the class structure generates a taxonomy 
that aligns with vocabularies of chemistry and clinical medicine, and the individual entities support quantitative assertions, 
such as half-life. Each of 4,901 decay branches is an individual carrying its mode of decay and branching ratio; 4,669 branches 
have an explicit link to their decay product (daughter) nuclide. 

Content was assembled from authoritative resources including the International Atomic Energy Agency (IAEA) Live Chart of Nuclides; 
symbols and names were verified against standard chemical nomenclature. 
NUCLEON is linked to the Ontology of Chemical Elements, SNOMED Clinical Terms, 
and DICOM radiopharmaceutical isotope and PET radionuclide code sequences.

The ontology incorporates 107,617 triples across 4,227 classes and 9,156 named individuals. Element names are provided in seven languages. 

* [NUCLEON ontology in Turtle format](NUCLEON.ttl)
* [Example SQARQL queries](nucleon_queries.rq)

NUCLEON is also available via [NCBO BioPortal](https://bioportal.bioontology.org/ontologies/NUCLEON).
