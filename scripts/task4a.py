#!/usr/bin/env python3
"""Task 4(a): blow-ups of C7, C9, Petersen at n=16 and n=18.

Enumerate positive part vectors up to Aut(H) (dihedral for cycles;
S_Petersen computed by backtracking). Compare L to L(K_bal).
"""
from __future__ import annotations

import itertools
import json
import sys
import time
from fractions import Fraction
from math import factorial
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from independent import (  # noqa: E402
    L_blowup,
    L_kbal,
    cycle_adj,
    d5_canonical,  # not used
    graph_automorphisms,
    petersen_adj,
)

ROOT = Path(__file__).resolve().parents[1]


def cyclic_images(m, p, reflections=True):
    m = tuple(m)
    out = []
    for i in range(p):
        rot = tuple(m[(i + j) % p] for j in range(p))
        out.append(rot)
        if reflections:
            out.append(tuple(rot[0] if j == 0 else rot[p - j] for j in range(p)))
            # standard reflection: reverse
            out.append(tuple(reversed(rot)))
    return out


def dihedral_canonical(m):
    p = len(m)
    imgs = []
    for i in range(p):
        rot = tuple(m[(i + j) % p] for j in range(p))
        imgs.append(rot)
        imgs.append(tuple(rot[(p - j) % p] for j in range(p)))
    return min(imgs)


def compositions_canonical(n, p, canonical_fn):
    reps = []
    for comb in itertools.combinations(range(n - 1), p - 1):
        cuts = (0,) + tuple(x + 1 for x in comb) + (n,)
        m = tuple(cuts[i + 1] - cuts[i] for i in range(p))
        if canonical_fn(m) == m:
            reps.append(m)
    return reps


def petersen_canonical(m, auts):
    m = tuple(m)
    return min(tuple(m[g[i]] for i in range(10)) for g in auts)


def scan(name, adj, n, canonical_fn):
    p = len(adj)
    t0 = time.perf_counter()
    reps = compositions_canonical(n, p, canonical_fn)
    Lk = L_kbal(n)
    best = None
    aut_H = graph_automorphisms(adj)
    records = []
    for i, m in enumerate(reps):
        L = L_blowup(adj, m, aut_H)
        rec = {"parts": list(m), "L": f"{L.numerator}/{L.denominator}"}
        records.append(rec)
        if best is None or L < best[0]:
            best = (L, rec)
        if (i + 1) % 50 == 0:
            print(f"  {name} n={n} {i+1}/{len(reps)} best_ratio={float(best[0]/Lk):.6g}", flush=True)
    elapsed = time.perf_counter() - t0
    Lbest = best[0]
    ratio = Lbest / Lk
    out = {
        "base": name,
        "n": n,
        "n_orbits": len(reps),
        "best_parts": best[1]["parts"],
        "L": best[1]["L"],
        "L_kbal": f"{Lk.numerator}/{Lk.denominator}",
        "beat": Lbest < Lk,
        "ratio_exact": f"{ratio.numerator}/{ratio.denominator}",
        "ratio_float": float(ratio),
        "seconds": elapsed,
        "n_aut_H": len(aut_H),
    }
    print(
        f"{name} n={n} orbits={len(reps)} best={out['best_parts']} "
        f"ratio={out['ratio_float']:.6g} beat={out['beat']} t={elapsed:.1f}s",
        flush=True,
    )
    return out


def main():
    outdir = ROOT / "data"
    outdir.mkdir(exist_ok=True)
    results = []
    for n in (16, 18):
        results.append(scan("C7", cycle_adj(7), n, dihedral_canonical))
        results.append(scan("C9", cycle_adj(9), n, dihedral_canonical))
        padj = petersen_adj()
        auts = graph_automorphisms(padj)
        results.append(
            scan("Petersen", padj, n, lambda m, a=auts: petersen_canonical(m, a))
        )
    path = outdir / "task4a_blowups.json"
    path.write_text(json.dumps(results, indent=2) + "\n")
    print("wrote", path)


if __name__ == "__main__":
    main()
