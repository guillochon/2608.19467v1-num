#!/usr/bin/env python3
"""All complete-multipartite graphs on n=16 and n=18 via the subset recurrence."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from independent import (  # noqa: E402
    L_kbal,
    L_subset,
    complete_multipartite_adj_bits,
    complete_multipartite_aut_order,
)


def partitions(n: int):
    out = []

    def rec(remain, maxp, acc):
        if remain == 0:
            if len(acc) >= 2:
                out.append(tuple(sorted(acc, reverse=True)))
            return
        for p in range(1, min(maxp, remain) + 1):
            rec(remain - p, p, acc + [p])

    rec(n, n, [])
    # unique
    return sorted(set(out))


def scan(n: int) -> dict:
    t0 = time.perf_counter()
    Lk = L_kbal(n)
    parts_list = partitions(n)
    print(f"n={n} partitions={len(parts_list)}", flush=True)
    best = None
    rows = []
    for i, m in enumerate(parts_list):
        adj = complete_multipartite_adj_bits(m)
        aut = complete_multipartite_aut_order(m)
        L = L_subset(n, adj, aut)
        ratio = L / Lk
        rec = {
            "parts": list(m),
            "aut": aut,
            "L": f"{L.numerator}/{L.denominator}",
            "ratio_float": float(ratio),
            "beat": L < Lk,
        }
        rows.append(rec)
        if best is None or L < best[0]:
            best = (L, rec)
        if (i + 1) % 30 == 0:
            print(
                f"  {i+1}/{len(parts_list)} best={best[1]['parts']} ratio={best[1]['ratio_float']:.6g}",
                flush=True,
            )
    elapsed = time.perf_counter() - t0
    rows.sort(key=lambda r: r["ratio_float"])
    return {
        "n": n,
        "n_partitions": len(parts_list),
        "best": best[1],
        "any_beat": any(r["beat"] for r in rows),
        "seconds": elapsed,
        "top8": rows[:8],
    }


def main():
    results = []
    for n in (16, 18):
        r = scan(n)
        results.append(r)
        print(
            f"DONE n={n} best={r['best']['parts']} ratio={r['best']['ratio_float']:.6g} "
            f"beat={r['any_beat']} t={r['seconds']:.1f}s",
            flush=True,
        )
    path = Path("/home/james/src/graph-likelihood/data/task4b_multipartite.json")
    path.write_text(json.dumps(results, indent=2) + "\n")
    print("wrote", path)


if __name__ == "__main__":
    main()
