# Manuscript figures

| Figure | Source | Script | Outputs |
|---|---|---|---|
| 1. Classes and properties | `template/Fig1_template.pptx` (author's style) | `make_class_figure.py` | `class_figure.pptx` (editable), `class_figure.pdf`, `class_figure.png` |
| 2. Construction and validation | — | `make_methods_figure.py` | `methods_figure.svg`, `.pdf`, `.png` |
| 3. Chart of nuclides | `NUCLEON.ttl` | `make_chart_figure.py` | `chart_figure.svg`, `.pdf`, `.png`, `chart_figure_counts.txt` |

Figure 1 is built by cloning the header, panel, arrow and label shapes of the
template slide. The PDF and PNG here were exported with LibreOffice
(`soffice --headless --convert-to pdf class_figure.pptx`), which substitutes a
font for Aptos; export from PowerPoint for the final version.

Figure 3 reads the ontology: each ground state is one cell, shaded by the decay
mode of its branch with the largest branching ratio (beta-delayed particle
emission grouped with its beta decay). Cells containing a nuclide with a
`:dicomCode` (ground state or isomer) are outlined.
