#!/usr/bin/env python3
"""Download nuclide data from the IAEA LiveChart of Nuclides data API.

Ground states (one request):
    python3 scripts/download_iaea.py

  Writes source-data/IAEA nuclide livechart.csv, the file NUCLEON's ground-state
  values are loaded from, after keeping a dated copy of the previous file. Then
  compare and reload:
    python3 scripts/compare_to_source.py
    python3 scripts/reload_from_source.py NUCLEON.ttl --write

Excited states, for comparison with NUBASE2020 isomers (one request per nuclide,
about 3,400 requests; resumable):
    python3 scripts/download_iaea.py --levels

  Writes source-data/iaea_levels_with_half_life.csv: every excited level that the
  API reports with a half-life, i.e. the candidate isomers.

API: https://nds.iaea.org/relnsd/vcharthtml/api_v0_guide.html
Data: ENSDF (evaluated nuclear structure data); masses from AME2020;
abundances from NUBASE2020. Each file ends each row with its extraction date,
which should be cited, e.g. "IAEA LiveChart of Nuclides, data extracted
18 October 2023".

Standard library only (plus pandas for --levels). Run from the repository root.
"""
import argparse
import csv
import datetime
import io
import os
import shutil
import sys
import time
import urllib.error
import urllib.request

BASE = "https://nds.iaea.org/relnsd/v1/data?"
# The API sometimes refuses requests without a browser-like User-Agent (HTTP 403).
HEADERS = {"User-Agent": "Mozilla/5.0 (NUCLEON ontology build; https://github.com/RSNA/NUCLEON)"}
GS_FILE = os.path.join("source-data", "IAEA nuclide livechart.csv")
LEVELS_FILE = os.path.join("source-data", "iaea_levels_with_half_life.csv")


def get(query, retries=4):
    req = urllib.request.Request(BASE + query, headers=HEADERS)
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read().decode("utf-8")
        except (urllib.error.URLError, TimeoutError) as e:
            if attempt == retries - 1:
                raise
            wait = 5 * (attempt + 1)
            print(f"  {e}; retrying in {wait} s", file=sys.stderr)
            time.sleep(wait)


def download_ground_states():
    text = get("fields=ground_states&nuclides=all")
    rows = list(csv.reader(io.StringIO(text)))
    if len(rows) < 3000 or "half_life_sec" not in rows[0]:
        sys.exit(f"Unexpected response ({len(rows)} rows); file not replaced.\n{text[:500]}")
    if os.path.exists(GS_FILE):
        stamp = datetime.date.today().isoformat()
        backup = GS_FILE.replace(".csv", f" (before {stamp}).csv")
        shutil.copyfile(GS_FILE, backup)
        print(f"previous file kept as {backup}")
    with open(GS_FILE, "w", encoding="utf-8", newline="") as f:
        f.write(text)
    header = rows[0]
    date_col = next((i for i, h in enumerate(header) if "extraction" in h.lower()), None)
    dates = sorted({r[date_col] for r in rows[1:] if date_col is not None and len(r) > date_col})
    print(f"wrote {GS_FILE}: {len(rows) - 1} ground states; extraction date(s): {', '.join(dates) or 'n/a'}")


def download_levels(pause):
    import pandas as pd
    gs = pd.read_csv(GS_FILE)
    nuclides = [f"{int(z + n)}{s.lower()}" for z, n, s in zip(gs.z, gs.n, gs.symbol) if s != "n"]
    done = set()
    if os.path.exists(LEVELS_FILE):
        done = set(pd.read_csv(LEVELS_FILE, usecols=["query"])["query"])
    first = not os.path.exists(LEVELS_FILE)
    with open(LEVELS_FILE, "a", encoding="utf-8", newline="") as out:
        w = None
        for k, nuc in enumerate(nuclides, 1):
            if nuc in done:
                continue
            text = get(f"fields=levels&nuclides={nuc}")
            df = pd.read_csv(io.StringIO(text)) if text.strip() else pd.DataFrame()
            if not df.empty and "half_life" in df.columns:
                energy = pd.to_numeric(df.get("energy"), errors="coerce").fillna(0)
                has_hl = df["half_life"].notna() & df["half_life"].astype(str).str.strip().ne("")
                keep = df[(energy > 0) & has_hl].copy()
                keep.insert(0, "query", nuc)
                if not keep.empty:
                    if w is None:
                        w = csv.DictWriter(out, fieldnames=list(keep.columns))
                        if first:
                            w.writeheader()
                    w.writerows(keep.to_dict("records"))
            if k % 100 == 0:
                print(f"  {k}/{len(nuclides)}")
                out.flush()
            time.sleep(pause)
    print(f"wrote {LEVELS_FILE}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--levels", action="store_true", help="also download excited levels with a half-life")
    p.add_argument("--pause", type=float, default=0.5, help="seconds between level requests (default 0.5)")
    a = p.parse_args()
    if not os.path.isdir("source-data"):
        sys.exit("Run from the repository root.")
    download_ground_states()
    if a.levels:
        download_levels(a.pause)


if __name__ == "__main__":
    main()
