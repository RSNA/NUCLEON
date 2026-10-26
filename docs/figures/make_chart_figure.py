#!/usr/bin/env python3
"""Draw Figure 3: NUCLEON coverage on the chart of nuclides (Z against N).

Each ground state is one cell, shaded by its dominant decay mode (the branch
with the largest branching ratio, grouped as below). Cells of nuclides that
carry a DICOM code (ground state or isomer) are outlined. Twelve uranium
placeholders (U-243 to U-255) are not plotted; see EXCLUDE. Reads NUCLEON.ttl;
writes chart_figure.svg, .pdf and .png (600 dpi at 180 mm width) and
chart_figure_counts.txt next to this script.

Usage: python3 docs/figures/make_chart_figure.py [NUCLEON.ttl]
"""
import collections
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import rdflib

TTL = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "..", "..", "NUCLEON.ttl")
HERE = os.path.dirname(os.path.abspath(__file__))

# decay-mode individual -> display group
GROUP = {}
for m in ("BetaMinusDecay", "DoubleBetaMinusDecay") + tuple(
        f"BetaMinusDelayed{x}" for x in ("NeutronEmission", "TwoNeutronEmission", "ThreeNeutronEmission",
                                         "FourNeutronEmission", "FiveNeutronEmission", "SixNeutronEmission",
                                         "SevenNeutronEmission", "ProtonEmission", "AlphaEmission")):
    GROUP[m] = "β⁻"
for m in ("BetaPlusDecay", "PositronEmission", "ElectronCapture", "DoubleBetaPlusDecay", "DoubleElectronCapture",
          "BetaPlusDelayedProtonEmission", "BetaPlusDelayedTwoProtonEmission", "BetaPlusDelayedAlphaEmission"):
    GROUP[m] = "β⁺ / EC"
GROUP["AlphaDecay"] = "α"
GROUP["ClusterDecay"] = "α"          # heavy-cluster emitters are all α emitters; never dominant
for m in ("ProtonEmission", "TwoProtonEmission"):
    GROUP[m] = "p"
for m in ("NeutronEmission", "TwoNeutronEmission"):
    GROUP[m] = "n"
for m in ("SpontaneousFission", "BetaPlusDelayedFission", "BetaMinusDelayedFission"):
    GROUP[m] = "SF"
GROUP["IsomericTransition"] = "IT"   # isomers only; ground states never use it

# colours: validated as an adjacency palette (spatial neighbours on the chart)
# with the dataviz validator; every region is also labelled directly
STYLE = [  # group, fill, legend label
    ("stable", "#1a1a1a", "Stable"),
    ("β⁻", "#2a78d6", r"$\beta^{-}$ decay"),
    ("β⁺ / EC", "#e34948", r"$\beta^{+}$ decay / electron capture"),
    ("α", "#eda100", "α decay"),
    ("SF", "#008300", "Spontaneous fission"),
    ("p", "#1baf7a", "Proton emission"),
    ("n", "#e87ba4", "Neutron emission"),
    ("unknown", "#c8c8c8", "Radioactive, mode not known"),
]
COLOR = {k: c for k, c, _ in STYLE}

g = rdflib.Graph()
g.parse(TTL)
Q = """
PREFIX : <https://rsna.org/radlex/nucleon#>
SELECT ?n ?z ?nn ?m ?st ?mode ?br ?dicom WHERE {
  ?n :atomicNumber ?z ; :neutronNumber ?nn .
  ?n a :Nuclide .
  OPTIONAL { ?n :metastableState ?m }
  OPTIONAL { ?n :isStable ?st }
  OPTIONAL { ?n :dicomCode ?dicom }
  OPTIONAL { ?n :hasDecayBranch ?b . ?b :hasDecayMode ?mode . OPTIONAL { ?b :branchingRatio ?br } }
}"""
cells = {}
best = {}
dicom = {}
for r in g.query(Q):
    key = (int(r.z), int(r.nn))
    name = r.n.split("#")[1]
    if r.dicom is not None:
        dicom.setdefault(key, set()).add(name)
    if r.m is not None:
        continue                      # isomers share the ground state's cell
    cells.setdefault(key, {"name": name, "stable": bool(r.st and r.st.toPython())})
    if r.mode is not None:
        br = float(r.br) if r.br is not None else -1.0
        grp = GROUP[r.mode.split("#")[1]]
        if key not in best or br > best[key][0]:
            best[key] = (br, grp)

# The IAEA export carries twelve uranium entries beyond U-242 (U-243 to U-255,
# no U-252) as empty ENSDF placeholders: no half-life, no discovery, absent from
# NUBASE2020 or estimated there. They stay in the ontology but would draw a
# spurious row past the end of the uranium chain, so they are not plotted.
EXCLUDE = {f"U-{a}" for a in range(243, 256)}

cat = {}
for key, c in cells.items():
    if c["name"] in EXCLUDE:
        continue
    cat[key] = "stable" if c["stable"] else (best[key][1] if key in best else "unknown")
counts = collections.Counter(cat.values())
assert "IT" not in counts, "a ground state is dominated by isomeric transition"

# ---- draw --------------------------------------------------------------------------
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 7, "axes.linewidth": 0.6,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "svg.fonttype": "none", "mathtext.fontset": "custom",
                     "mathtext.rm": "Liberation Sans", "mathtext.it": "Liberation Sans:italic",
                     "mathtext.bf": "Liberation Sans:bold"})
W_MM = 180
fig, ax = plt.subplots(figsize=(W_MM / 25.4, W_MM / 25.4 * 0.70))
for (z, n), k in cat.items():
    ax.add_patch(Rectangle((n - 0.5, z - 0.5), 1, 1, facecolor=COLOR[k], edgecolor="white", linewidth=0.12))
for (z, n) in dicom:
    ax.add_patch(Rectangle((n - 0.5, z - 0.5), 1, 1, facecolor="none", edgecolor="#000000", linewidth=0.9))

# magic numbers
for mn in (2, 8, 20, 28, 50, 82, 126):
    ax.axvline(mn, color="#8c8c8c", lw=0.35, ls=(0, (2, 2)), zorder=0)
    if mn <= 118:
        ax.axhline(mn, color="#8c8c8c", lw=0.35, ls=(0, (2, 2)), zorder=0)

maxn = max(n for _, n in cat)
ax.set_xlim(-1, maxn + 2)
ax.set_ylim(-1, 120)
ax.set_aspect("equal")
ax.set_xlabel("Neutron number N")
ax.set_ylabel("Atomic number Z")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.set_xticks(range(0, maxn + 1, 20))
ax.set_yticks(range(0, 121, 20))
ax.tick_params(length=2.5, colors="#333333")

# direct labels for regions (secondary encoding: identity never colour alone)
def lab(x, y, s, **kw):
    ax.text(x, y, s, fontsize=7, color="#1a1a1a", **kw)
lab(118, 52, r"$\mathbf{\beta^{-}}$", ha="center")
lab(40, 58, r"$\mathbf{\beta^{+}}$ / EC", ha="center")
for txt, xy, xyt in (("α", (124, 90), (112, 99)), ("SF", (157, 104), (170, 96))):
    ax.annotate(txt, xy=xy, xytext=xyt, fontsize=7, fontweight="bold", ha="center", va="center",
                arrowprops=dict(arrowstyle="-", lw=0.5, color="#333333", shrinkA=2, shrinkB=1))

# selected DICOM-coded radionuclides
LABELS = {"F-18": (9, 9), "Ga-68": (31, 37), "Tc-99": (43, 56), "I-131": (53, 78),
          "Lu-177": (71, 106), "Ac-225": (89, 136)}
TEXT = {"F-18": "¹⁸F", "Ga-68": "⁶⁸Ga", "Tc-99": "⁹⁹ᵐTc", "I-131": "¹³¹I", "Lu-177": "¹⁷⁷Lu", "Ac-225": "²²⁵Ac"}
OFF = {"F-18": (-7, 17), "Ga-68": (-16, 9), "Tc-99": (-18, 10), "I-131": (14, -12),
       "Lu-177": (16, -10), "Ac-225": (14, -12)}
for k, (z, n) in LABELS.items():
    assert (z, n) in dicom, k
    dx, dy = OFF[k]
    ax.annotate(TEXT[k], xy=(n, z), xytext=(n + dx, z + dy), fontsize=6.5, ha="center", va="center",
                arrowprops=dict(arrowstyle="-", lw=0.5, color="#333333", shrinkA=1, shrinkB=2.5))

# legend, top left
from matplotlib.patches import Patch
handles = [Patch(facecolor=c, edgecolor="none", label=f"{l} ({counts[k]:,})") for k, c, l in STYLE if counts[k]]
handles.append(Patch(facecolor="white", edgecolor="black", linewidth=0.9,
                     label=f"Carries a DICOM code ({len(dicom)} cells)"))
leg = ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=6.5, handlelength=1.1,
                handleheight=1.1, borderaxespad=0.2, title=f"Ground states by dominant decay mode (n = {len(cat):,})",
                title_fontsize=6.8)
leg._legend_box.align = "left"

fig.tight_layout(pad=0.3)
p = os.path.join(HERE, "chart_figure")
fig.savefig(p + ".svg")
fig.savefig(p + ".pdf")
fig.savefig(p + ".png", dpi=600)
report = [f"ground states plotted: {len(cat)} of {len(cells)}",
          f"not plotted: {', '.join(sorted(c['name'] for c in cells.values() if c['name'] in EXCLUDE))}", f"cells with a DICOM-coded nuclide: {len(dicom)}",
          f"DICOM-coded nuclides: {sum(len(v) for v in dicom.values())}"]
report += [f"{l}: {counts[k]}" for k, _, l in STYLE]
open(p + "_counts.txt", "w").write("\n".join(report) + "\n")
print("\n".join(report))
