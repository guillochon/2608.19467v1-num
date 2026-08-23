#!/usr/bin/env python3
import json
from fractions import Fraction
from pathlib import Path

root = Path("/home/james/src/graph-likelihood")
rows = json.loads((root / "data/task1_n15_40.json").read_text())
rows += json.loads((root / "data/task1_n41_50.json").read_text())
for n in (15, 16, 17, 18, 19, 23, 25, 28, 33, 40, 50):
    r = next(x for x in rows if x["n"] == n)
    Lc = Fraction(r["L_c5"])
    Lk = Fraction(r["L_kbal"])
    ratio = Lc / Lk
    print(
        f"n={n} parts={r['best_parts']} beat={r['beat']} "
        f"ratio={float(ratio):.10g} exact={ratio.numerator}/{ratio.denominator}"
    )
    print(f"  L_c5  = {r['L_c5']}")
    print(f"  L_kbal= {r['L_kbal']}")
