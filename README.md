# NUCLEON &mdash; Nuclide and Element Ontology

The Nuclide and Element Ontology (NUCLEON) is a machine-processable knowledge resource that represents the 118 known chemical elements and 4,105 nuclides, including 716 metastable isomers, with computable nuclear properties.

Each element and nuclide is declared as both a class and an individual (OWL 2 punning). The class structure provides a taxonomy that aligns with vocabularies of chemistry and clinical medicine; the individuals carry the quantitative assertions, such as half-life. Each of 5,872 decay branches is an individual carrying its mode of decay and, where the source reports one, its branching ratio; 5,629 branches link explicitly to their decay product (daughter) nuclide.

Content was assembled from authoritative resources including the International Atomic Energy Agency (IAEA) Live Chart of Nuclides; symbols and names were verified against IUPAC nomenclature. NUCLEON is linked to the Ontology of Chemical Elements, SNOMED CT, and the DICOM radiopharmaceutical isotope (CID 18) and PET radionuclide (CID 4020) context groups, so a nuclide can be resolved directly from a DICOM Radionuclide Code Sequence. Element names are provided in seven languages.

NUCLEON is also available from [NCBO BioPortal](https://bioportal.bioontology.org/ontologies/NUCLEON).

## Contents

| Path | Description |
| --- | --- |
| [`NUCLEON.ttl`](NUCLEON.ttl) | The ontology, OWL 2 DL in Turtle |
| [`nucleon_queries.rq`](nucleon_queries.rq) | Example SPARQL queries (Q1&ndash;Q7) and validation queries (V1&ndash;V5) |
| [`scripts/validate.py`](scripts/validate.py) | Runs every query, reports content counts and structural checks |
| [`scripts/add_isomer_branches.py`](scripts/add_isomer_branches.py) | Derives decay branches for metastable isomers (used for v1.1) |
| [`validation/report.txt`](validation/report.txt) | Output of `validate.py` for the current release |
| [`source-data/`](source-data) | IAEA LiveChart export and the decay-mode reference table |

## Content (v1.1)

| Entity or property | Count |
| --- | ---: |
| Triples | 114,496 |
| Elements | 118 |
| Nuclides | 4,105 |
| &nbsp;&nbsp;metastable isomers | 716 |
| &nbsp;&nbsp;flagged stable | 249 |
| &nbsp;&nbsp;with a half-life | 3,637 |
| &nbsp;&nbsp;with at least one decay branch | 3,744 |
| Decay branches | 5,872 |
| &nbsp;&nbsp;with a decay product | 5,629 |
| Decay modes | 32 |
| `skos:exactMatch` assertions | 287 |
| Nuclides with a DICOM code | 87 |

Of the 243 branches without a decay product, 207 are fission branches (products vary by event), 5 bundle fission with electron capture/beta-plus decay, and 31 lead to a daughter not covered by the source data. Of the 361 nuclides without a branch, 234 are stable and 127 are radioactive with no decay mode reported by the source (50 ground states, 77 isomers). Fifteen primordial nuclides flagged stable (e.g. Bi-209, W-180) also carry the alpha or double-beta branch the source lists for them.

## Validation

```
pip install rdflib
python3 scripts/validate.py
```

The script runs every query in `nucleon_queries.rq` over asserted triples only (no reasoner is needed) and exits non-zero if any validation query or structural check fails. V1 confirms nucleon conservation across every asserted daughter link; because daughters are derived with that same arithmetic, it verifies internal consistency, not independent physical correctness. The ontology classifies as consistent with no unsatisfiable classes under HermiT.

## Using the ontology

- Physical values are asserted on the **individual** reading of each IRI. Protégé, the OWL API, rdflib, and SPARQL endpoints read both readings. **owlready2 reads a punned IRI only as a class** and silently returns no property values; query with rdflib or SPARQL instead.
- `:halfLife` is always in seconds (`xsd:double`); the source value and unit are kept in `:halfLifeValue` and `:halfLifeUnit`.
- Following the source's convention, delayed-particle channels (e.g. beta-delayed neutron emission) are sub-channels of their parent beta decay, so a nuclide's branching ratios can sum to more than 100.
- Uncertainties are not represented; values are central values. Limit or approximate branching ratios carry an `rdfs:comment`.
- Nuclides used in nuclear medicine can be selected by `:dicomContextGroup` (CID 18 or CID 4020).

## Changes in v1.1

- Added 971 decay branches for 636 metastable isomers (e.g. Kr-81m, In-113m, Tc-94m, Mn-52m), parsed from the source decay text; daughters derived as for ground states.
- Tc-99m: isomeric transition corrected to 99.9963%; 0.0037% beta-minus branch to Ru-99 added.
- Removed `:isMedicalIsotope`, whose criterion was undocumented; use `:dicomContextGroup`.
- Added license, validation script, and validation report.

## License

NUCLEON is released under the [Creative Commons Attribution 4.0 International license](LICENSE) (CC BY 4.0). SNOMED CT identifiers are included as cross-references only; use of SNOMED CT content is subject to the SNOMED International license.
