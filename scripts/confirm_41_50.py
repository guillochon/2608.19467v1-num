#!/usr/bin/env python3
import json, subprocess, time
from pathlib import Path
root = Path("/home/james/src/graph-likelihood")
rows = json.loads((root / "data/task1_n41_50.json").read_text())
out = []
for row in rows:
    n = row["n"]
    parts = row["best_parts"]
    t0 = time.perf_counter()
    r1 = subprocess.run(
        [str(root / "likelihood"), "c5parts", *map(str, parts)],
        capture_output=True, text=True, check=True, cwd=root,
    )
    r2 = subprocess.run(
        [str(root / "likelihood"), "kbal", str(n)],
        capture_output=True, text=True, check=True, cwd=root,
    )
    Lc = r1.stdout.strip().split("likelihood=")[1]
    Lk = r2.stdout.strip().split("likelihood=")[1]
    ok = Lc == row["L_c5"] and Lk == row["L_kbal"]
    dt = time.perf_counter() - t0
    print(f"n={n} {'OK' if ok else 'FAIL'} t={dt:.2f}s", flush=True)
    out.append({"n": n, "ok": ok, "seconds": dt})
(root / "data/task1_C_confirm_41_50.json").write_text(json.dumps(out, indent=2) + "\n")
print("all_ok", all(x["ok"] for x in out))
