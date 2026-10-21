# Source data

| File | Content | Source | Used for |
| --- | --- | --- | --- |
| `IAEA nuclide livechart.csv` | Ground states of all known nuclides | IAEA LiveChart of Nuclides data API, `fields=ground_states&nuclides=all` (https://nds.iaea.org/relnsd/vcharthtml/api_v0_guide.html). ENSDF evaluations; masses AME2020; abundances NUBASE2020. Extracted 18 October 2023 (last column). | Ground-state half-lives, widths, spin-parities, abundances, decay modes, branching ratios, Q-values |
| `nubase_3.mas20.txt` | NUBASE2020 table: ground states, isomers and levels | Kondev FG, Wang M, Huang WJ, Naimi S, Audi G. The NUBASE2020 evaluation of nuclear physics properties. *Chinese Physics C* 2021;45:030001. doi:10.1088/1674-1137/abddae. Atomic Mass Data Center distribution file, as redistributed in https://github.com/radioactivedecay/datasets | All metastable isomers; ground-state values the IAEA export lacks |
| `decay-modes.txt` | Decay-mode definitions and their changes in Z and A | NUBASE2020 (decay-mode notation) | `:changeInAtomicNumber` and `:changeInMassNumber` of the 32 decay modes |

Refresh the IAEA file with `python3 scripts/download_iaea.py` (keeps a dated copy of the previous file), then rerun the steps in `validation/README.md`.
