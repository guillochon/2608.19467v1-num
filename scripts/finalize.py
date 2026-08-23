#!/usr/bin/env python3
import hashlib, json, subprocess
from pathlib import Path
from fractions import Fraction
root = Path("/home/james/src/graph-likelihood")
(root / "code").mkdir(exist_ok=True)
for src in ["likelihood.c", "Makefile"]:
    data = (root / src).read_bytes()
    (root / "code" / src).write_bytes(data)
rows = json.loads((root / "data/task1_n15_40.json").read_text())
rows += json.loads((root / "data/task1_n41_50.json").read_text())
tsv = root / "data/task1_witness_table.tsv"
with tsv.open("w") as f:
    f.write("n\tbest_parts\tbeat\tratio_float\tL_c5\tL_kbal\tratio_exact\n")
    for r in sorted(rows, key=lambda x: x["n"]):
        Lc = Fraction(r["L_c5"])
        Lk = Fraction(r["L_kbal"])
        ratio = Lc / Lk
        parts = ",".join(map(str, r["best_parts"]))
        f.write(
            f"{r['n']}\t{parts}\t{int(r['beat'])}\t{float(ratio):.16g}\t"
            f"{r['L_c5']}\t{r['L_kbal']}\t{ratio.numerator}/{ratio.denominator}\n"
        )
print("wrote", tsv)
print("--- sha256 ---")
for p in [
    "likelihood.c",
    "code/independent.py",
    "Makefile",
    "RESULTS.md",
    "LETTER.md",
]:
    h = hashlib.sha256((root / p).read_bytes()).hexdigest()
    print(h, p)
print("--- git ---")
subprocess.check_call(["git", "add", "-A"], cwd=root)
subprocess.check_call(
    [
        "git",
        "-c",
        "user.email=james@local",
        "-c",
        "user.name=james",
        "commit",
        "-m",
        "Exact C5-family table n=15-50, tower n=80, n=16/18 searches\n\n"
        "Independent Python evaluator plus C kbal/graph6. Witness table, "
        "growth-constant fit, C7/C9/Petersen blow-ups, complete-multipartite n=16.",
    ],
    cwd=root,
)
print(subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip())
print(subprocess.check_output(["git", "status", "--short"], cwd=root, text=True))
