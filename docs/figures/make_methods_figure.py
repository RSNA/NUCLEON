#!/usr/bin/env python3
"""Draw the NUCLEON methods figure (pipeline from sources to release).

Writes methods_figure.svg, .pdf and .png (600 dpi at 180 mm width) next to
this script. Numbers are those of NUCLEON v1.6 and validation/README.md.
"""
import os

import cairosvg

W, H = 1300, 1040
FONT = "Liberation Sans, Arial, Helvetica, sans-serif"
MONO = "Liberation Mono, Courier New, monospace"
INK, MUTED, LINE, FILL, ACCENT, ACCENT_FILL = "#1a1a1a", "#555555", "#8c8c8c", "#f4f4f4", "#1f5fa8", "#e8f0fa"

out = []


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def text(x, y, s, size=17, weight="normal", color=INK, family=FONT, anchor="start", italic=False):
    st = ' font-style="italic"' if italic else ""
    out.append(f'<text x="{x}" y="{y}" font-family="{family}" font-size="{size}" font-weight="{weight}" '
               f'fill="{color}" text-anchor="{anchor}"{st}>{s}</text>')


def box(x, y, w, h, title, lines, accent=False, num=None):
    out.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="6" fill="{ACCENT_FILL if accent else FILL}" '
               f'stroke="{ACCENT if accent else LINE}" stroke-width="{2 if accent else 1.2}"/>')
    tx = x + 14
    if num is not None:
        out.append(f'<circle cx="{x + 24}" cy="{y + 25}" r="13" fill="{ACCENT if accent else INK}"/>')
        text(x + 24, y + 31, str(num), 16, "bold", "#ffffff", anchor="middle")
        tx = x + 46
    text(tx, y + 31, esc(title), 19, "bold")
    yy = y + 57
    for ln in lines:
        if isinstance(ln, tuple):      # (string, options)
            s, opt = ln
            text(x + 14, yy, s, **opt)
        else:
            text(x + 14, yy, ln, 17, color=INK)
        yy += 22


def arrow(x1, y1, x2, y2, color=LINE, width=2, dashed=False):
    d = ' stroke-dasharray="6 5"' if dashed else ""
    out.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{width}"{d} '
               f'marker-end="url(#ah{"A" if color == ACCENT else ""})"/>')


def header(x, w, label):
    text(x, 34, label, 21, "bold", MUTED)
    out.append(f'<line x1="{x}" y1="46" x2="{x + w}" y2="46" stroke="{MUTED}" stroke-width="1.5"/>')


# ---- columns -----------------------------------------------------------------------
SX, SW = 20, 300     # sources
BX, BW = 370, 460    # build
VX, VW = 880, 400    # validation
header(SX, SW, "Sources")
header(BX, BW, "Build (scripted, reproducible)")
header(VX, VW, "Validation")

GAP = 16


def hgt(n):          # box height for n body lines
    return 57 + 22 * (n - 1) + 15


steps = [
    ("Load ground states", ["Half-life in seconds (\u0127 ln2 / \u0393 for resonances);",
                            "spin-parity, abundance, decay modes."], False),
    ("Generate isomers", ["NUBASE state index \u2192 Tc-99m, Lu-177m3 \u2026;",
                          "excitation energy, half-life, decay modes."], False),
    ("Create decay branches, derive daughters",
     [("One :DecayBranch per decay mode, with its ratio.", {}),
      ("Z\u2032 = Z + \u0394Z,  A\u2032 = A + \u0394A  \u2192  daughter ground state", {"size": 18, "weight": "bold", "color": ACCENT}),
      ("No daughter for fission; cluster decay subtracts", {}),
      ("the cluster. 7,399 of 7,691 branches linked.", {})], True),
    ("Add emissions", ["Type, energy and yield per decay (\u2265 0.1%, \u2265 1 keV)",
                       "for DICOM-coded nuclides and their decay chains."], False),
    ("Map external codes", ["skos:exactMatch to OCE and SNOMED CT;",
                            ":dicomCode as in the Radionuclide Code Sequence."], False),
    ("Mirror values onto classes", ["Display copies on each punned class",
                                    "(:displayValue) for viewers such as BioPortal."], False),
]
rows, y = [], 64
for title, lines, acc in steps:
    rows.append((y, hgt(len(lines))))
    y += hgt(len(lines)) + GAP
BOTTOM = rows[-1][0] + rows[-1][1]
for k, ((y, h), (title, lines, acc)) in enumerate(zip(rows, steps), 1):
    box(BX, y, BW, h, title, lines, accent=acc, num=k)
for k in range(len(rows) - 1):
    y1 = rows[k][0] + rows[k][1]
    arrow(BX + 24, y1 + 1, BX + 24, rows[k + 1][0] - 2)

# sources, aligned with the step each feeds
sources = [
    (0, "IAEA LiveChart of Nuclides", ["ENSDF evaluations; export of", "18 October 2023"]),
    (1, "NUBASE2020", ["Kondev et al. 2021;", "2,019 metastable isomers"]),
    (2, "Decay-mode definitions", ["\u0394Z and \u0394A for 32 modes", "(NUBASE2020 notation)"]),
    (3, "ICRP Publication 107", ["Nuclear decay data for", "dosimetric calculations (2008)"]),
    (4, "IUPAC \u00b7 OCE \u00b7 SNOMED CT", ["DICOM PS3.16: CID 18,", "CID 4020 (release 2026c)"]),
]
for row, title, lines in sources:
    y, h = rows[row]
    sh = hgt(2)
    sy = y + (h - sh) / 2
    box(SX, sy, SW, sh, title, lines)
    arrow(SX + SW + 2, sy + sh / 2, BX - 3, sy + sh / 2, ACCENT if row == 2 else LINE)

# validation and release, bottom-aligned with the build column
vboxes = [
    ("Internal consistency", ["SPARQL checks V1\u2013V6: nucleon conservation,",
                              "A = Z + N, no self-decay, no dangling links,",
                              "element Z, display copies match. HermiT:",
                              "consistent, no unsatisfiable classes."]),
    ("Daughters vs independent data", ["ICRP-107: 1,510 links agree",
                                       "ENSDF: 1,688 links agree",
                                       ("0 disagreements", {"weight": "bold", "color": ACCENT})]),
    ("Values vs sources", ["IAEA export: 6 documented exceptions",
                           "Gamma lines vs ENSDF: 204 of 206 agree"]),
    ("Manual audit", ["Stratified random sample of 100 nuclides,",
                      "two independent reviewers"]),
]
REL = ("Release: NUCLEON.ttl (OWL 2 DL)",
       ["118 elements \u00b7 5,408 nuclides (2,019 isomers)",
        "7,691 decay branches \u00b7 3,990 emissions",
        ("BioPortal \u00b7 GitHub: scripts, queries, reports", {"color": MUTED})])
rel_h = hgt(3)
RY = BOTTOM - rel_h
vy = 64
vend = RY - 34
avail = vend - vy - GAP * (len(vboxes) - 1)
need = sum(hgt(len(l)) for _, l in vboxes)
pad = (avail - need) / len(vboxes)
for title, lines in vboxes:
    h = hgt(len(lines)) + pad
    box(VX, vy, VW, h, title, lines)
    vy += h + GAP
box(VX, RY, VW, rel_h, REL[0], REL[1])
arrow(VX + VW / 2, vend + 2, VX + VW / 2, RY - 3, LINE, 2.5)
# build -> validation bracket arrow
bx = BX + BW
out.append(f'<path d="M {bx + 6} {rows[0][0] + 10} H {bx + 20} V {BOTTOM - 10} H {bx + 6}" '
           f'fill="none" stroke="{LINE}" stroke-width="2"/>')
arrow(bx + 20, 300, VX - 3, 300, LINE, 2.5)

# ---- worked example ------------------------------------------------------------------
EY = BOTTOM + 52
text(SX, EY + 4, "Worked example: deriving the daughters of Bi-213", 21, "bold", MUTED)
out.append(f'<line x1="{SX}" y1="{EY + 16}" x2="{W - 20}" y2="{EY + 16}" stroke="{MUTED}" stroke-width="1.5"/>')


def nuc(cx, cy, name, zline, sub, accent=False):
    w, h = 190, 92
    out.append(f'<rect x="{cx - w / 2}" y="{cy - h / 2}" width="{w}" height="{h}" rx="8" '
               f'fill="{ACCENT_FILL if accent else "#ffffff"}" stroke="{ACCENT if accent else INK}" stroke-width="1.8"/>')
    text(cx, cy - 12, name, 22, "bold", anchor="middle")
    text(cx, cy + 12, zline, 17, anchor="middle")
    text(cx, cy + 34, sub, 15, color=MUTED, anchor="middle")


py = EY + 112
nuc(140, py, "Bi-213", "Z = 83, A = 213", "t½ 45.59 min")
nuc(700, py - 52, "Po-213", "Z = 84, A = 213", "derived daughter", accent=True)
nuc(700, py + 52, "Tl-209", "Z = 81, A = 209", "derived daughter", accent=True)
# branches
for dy, mode, pct, rule in ((-52, "β− decay", "97.86%", "ΔZ = +1, ΔA = 0"),
                            (52, "α decay", "2.14%", "ΔZ = −2, ΔA = −4")):
    y = py + dy
    out.append(f'<path d="M 235 {py} C 300 {py}, 300 {y}, 360 {y} H 598" fill="none" stroke="{ACCENT}" '
               f'stroke-width="2.2" marker-end="url(#ahA)"/>')
    text(380, y - 10, f"{mode}, {pct}", 18, "bold", ACCENT)
    text(380, y + 20, rule, 17, color=INK)
text(830, py - 30, "Each daughter is computed from the", 17)
text(830, py - 8, "parent's Z and A and the mode's ΔZ and ΔA,", 17)
text(830, py + 14, "never transcribed, so every link can be", 17)
text(830, py + 36, "re-derived and checked (query V1).", 17)

svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">'
       '<defs>'
       f'<marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
       f'<path d="M0,0 L10,5 L0,10 z" fill="{LINE}"/></marker>'
       f'<marker id="ahA" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
       f'<path d="M0,0 L10,5 L0,10 z" fill="{ACCENT}"/></marker>'
       '</defs>'
       f'<rect width="{W}" height="{H}" fill="#ffffff"/>' + "".join(out) + "</svg>")

here = os.path.dirname(os.path.abspath(__file__))
p = os.path.join(here, "methods_figure")
open(p + ".svg", "w", encoding="utf-8").write(svg)
cairosvg.svg2pdf(bytestring=svg.encode(), write_to=p + ".pdf")
px_w = round(180 / 25.4 * 600)   # 180 mm at 600 dpi
cairosvg.svg2png(bytestring=svg.encode(), write_to=p + ".png", output_width=px_w, output_height=round(px_w * H / W))
print("wrote", p + ".{svg,pdf,png}", px_w, "px wide")
