#!/bin/bash
set -euo pipefail
cd /home/james/src/graph-likelihood
mkdir -p data logs
echo "=== independent table 15-40 ==="
/usr/bin/time -f 'elapsed %e sec maxrss %M KB' python3 code/independent.py table --nmin 15 --nmax 40 --out data/task1_n15_40.json
echo "=== C confirmation of each winner + K_bal ==="
python3 - <<'PY'
import json, subprocess, time, os
from fractions import Fraction
rows = json.load(open("data/task1_n15_40.json"))
out = []
for row in rows:
    n = row["n"]
    parts = row["best_parts"]
    t0 = time.perf_counter()
    r1 = subprocess.run(["./likelihood","c5parts",*map(str,parts)], capture_output=True, text=True, check=True)
    r2 = subprocess.run(["./likelihood","kbal",str(n)], capture_output=True, text=True, check=True)
    dt = time.perf_counter() - t0
    Lc = r1.stdout.strip().split("likelihood=")[1]
    Lk = r2.stdout.strip().split("likelihood=")[1]
    ok_c5 = Lc == row["L_c5"]
    ok_k = Lk == row["L_kbal"]
    rec = {"n": n, "parts": parts, "C_L_c5": Lc, "C_L_kbal": Lk, "match_c5": ok_c5, "match_kbal": ok_k, "seconds": dt}
    out.append(rec)
    flag = "OK" if ok_c5 and ok_k else "MISMATCH"
    print(f"n={n:3d} {flag} C5={ok_c5} Kbal={ok_k} t={dt:.2f}s")
json.dump(out, open("data/task1_C_confirm.json","w"), indent=2)
print("wrote data/task1_C_confirm.json")
if not all(r["match_c5"] and r["match_kbal"] for r in out):
    raise SystemExit(1)
PY
echo "=== thresholds 16 18 ==="
python3 code/independent.py threshold 16 18 | tee data/task4_thresholds.txt
