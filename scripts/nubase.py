"""Reader for the NUBASE2020 table (source-data/nubase_3.mas20.txt).

Kondev FG, Wang M, Huang WJ, Naimi S, Audi G. The NUBASE2020 evaluation of
nuclear physics properties. Chinese Physics C 2021;45:030001.
doi:10.1088/1674-1137/abddae

Column layout from the file header (1-based): A 1-3, Z 5-7, state index i 8,
A El 12-16, s 17, mass excess 19-31, excitation energy 43-54, its uncertainty
55-65, half-life 70-78, unit 79-80, half-life uncertainty 82-88, spin and
parity 89-102, ..., decay modes 120-.
"""
import math
import re

YEAR_S = 365.2422 * 86400  # NUBASE's year
UNIT_S = {"ys": 1e-24, "zs": 1e-21, "as": 1e-18, "fs": 1e-15, "ps": 1e-12, "ns": 1e-9, "us": 1e-6,
          "ms": 1e-3, "s": 1.0, "m": 60.0, "h": 3600.0, "d": 86400.0, "y": YEAR_S,
          "ky": 1e3 * YEAR_S, "My": 1e6 * YEAR_S, "Gy": 1e9 * YEAR_S, "Ty": 1e12 * YEAR_S,
          "Py": 1e15 * YEAR_S, "Ey": 1e18 * YEAR_S, "Zy": 1e21 * YEAR_S, "Yy": 1e24 * YEAR_S}
# NUBASE unit -> NUCLEON :halfLifeUnit (NUCLEON uses "a" for years and "µs")
UNIT_NUCLEON = {"us": "µs", "y": "a"}
ISOMER_STATES = {1, 2, 3, 4, 5, 6}   # 8, 9 are isobaric analog states


def suffix(i):
    return "" if i == 0 else "m" if i == 1 else f"m{i}"


def parse(path, zsym):
    """Return {nucleon_name: record} for ground states and isomers with a half-life."""
    out = {}
    for line in open(path, encoding="latin-1"):
        if line.startswith("#") or len(line) < 80:
            continue
        try:
            a, z, i = int(line[0:3]), int(line[4:7]), int(line[7])
        except ValueError:
            continue
        if z not in zsym or (i != 0 and i not in ISOMER_STATES):
            continue
        t_raw = line[69:78].strip()
        unit = line[78:80].strip()
        exc_raw = line[42:54].strip()
        rec = {
            "Z": z, "A": a, "i": i, "symbol": zsym[z],
            "name": f"{zsym[z]}-{a}{suffix(i)}",
            "nubase_label": line[11:17].rstrip(),
            "exc_keV": None, "exc_systematic": "#" in exc_raw,
            "stable": "stbl" in t_raw,
            "hl_value": None, "hl_unit": None, "hl_s": None, "hl_op": None, "hl_systematic": "#" in t_raw,
            "spin_parity": line[88:102].strip().replace("*", "").replace("#", "").strip() or None,
            "decay": line[119:].strip(),
        }
        if exc_raw:
            m = re.match(r"[-\d.]+", exc_raw.replace("#", ""))
            if m:
                rec["exc_keV"] = m.group(0)
        tv = t_raw.replace("stbl", "").replace("#", "").strip()
        if tv:
            op = None
            for o in (">=", "<=", ">", "<", "~"):
                if tv.startswith(o):
                    op, tv = o, tv[len(o):].strip()
                    break
            try:
                float(tv)
            except ValueError:
                tv = ""
            if tv and unit in UNIT_S:
                rec.update(hl_value=tv, hl_unit=unit, hl_s=float(tv) * UNIT_S[unit], hl_op=op)
        if i != 0 and rec["hl_s"] is None and not rec["stable"]:
            continue          # a level with no half-life is not a confirmed isomer
        if "p-unst" in t_raw:
            continue
        out[rec["name"]] = rec
    return out


# NUBASE decay-mode codes -> NUCLEON :DecayMode local names
MODE = {
    "A": "AlphaDecay", "B-": "BetaMinusDecay", "B+": "BetaPlusDecay", "EC": "ElectronCapture",
    "e+": "PositronEmission", "IT": "IsomericTransition", "SF": "SpontaneousFission",
    "p": "ProtonEmission", "2p": "TwoProtonEmission", "n": "NeutronEmission", "2n": "TwoNeutronEmission",
    "2B-": "DoubleBetaMinusDecay", "2B+": "DoubleBetaPlusDecay", "2EC": "DoubleElectronCapture",
    "B-n": "BetaMinusDelayedNeutronEmission", "B-2n": "BetaMinusDelayedTwoNeutronEmission",
    "B-3n": "BetaMinusDelayedThreeNeutronEmission", "B-4n": "BetaMinusDelayedFourNeutronEmission",
    "B-5n": "BetaMinusDelayedFiveNeutronEmission", "B-6n": "BetaMinusDelayedSixNeutronEmission",
    "B-7n": "BetaMinusDelayedSevenNeutronEmission", "B-p": "BetaMinusDelayedProtonEmission",
    "B-d": "BetaMinusDelayedDeuteronEmission", "B-t": "BetaMinusDelayedTritonEmission",
    "B-A": "BetaMinusDelayedAlphaEmission", "B-SF": "BetaMinusDelayedFission",
    "B+p": "BetaPlusDelayedProtonEmission", "B+2p": "BetaPlusDelayedTwoProtonEmission",
    "B+3p": "BetaPlusDelayedThreeProtonEmission", "B+A": "BetaPlusDelayedAlphaEmission",
    "B+SF": "BetaPlusDelayedFission", "ECp": "BetaPlusDelayedProtonEmission",
    "ECA": "BetaPlusDelayedAlphaEmission", "ECSF": "BetaPlusDelayedFission",
}
CLUSTER = re.compile(r"^(\d+)([A-Z][a-z]?)$")   # e.g. 14C, 24Ne


def decay_modes(s):
    """Parse a NUBASE decay string, e.g. 'IT~100;B-=0.0037 6' or 'A=91 3;B+ ?'.

    Returns [(mode_or_cluster, op, value_str_or_None)], unknown codes as ('?CODE', ...).
    Cluster decays are returned as ('CD:<A><El>', op, value)."""
    out = []
    s = re.sub(r"\[[^\]]*\]", "", s)          # drop bracketed notes such as [gs=0,m=100]
    for part in s.split(";"):
        part = part.strip()
        if not part or part.startswith("IS"):
            continue
        m = re.match(r"([0-9A-Za-z+\-]+?)\s*(=|~|<=|>=|<|>|\s\?|\?)\s*([0-9.]+(?:[eE][-+]?\d+)?)?", part + " ")
        if not m:
            out.append(("?" + part, None, None))
            continue
        code, op, val = m.group(1), m.group(2).strip(), m.group(3)
        if op == "?":
            val = None
        if code in MODE:
            out.append((MODE[code], op, val))
        elif CLUSTER.match(code):
            out.append(("CD:" + code, op, val))
        else:
            out.append(("?" + code, op, val))
    return out


def log_ratio(a, b):
    return abs(math.log10(a / b))
