#!/usr/bin/env python3
"""Blow-ups of twin-free bases of order b at n=16 (paper §5, through b=8).

Coverage: every graph with a twin pair is an independent-set blow-up of a
smaller twin-free base (true twins are complements of false-twin blow-ups,
and L(G)=L(complement)). Over-generation of Aut(H)-images of part vectors
is harmless. We skip a vector when prod(m_i!) * |Aut(H)| < T, which is a
sound upper bound on |Aut| of the blow-up of a twin-free H.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from fractions import Fraction
from itertools import combinations
from math import factorial
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from independent import (  # noqa: E402
    L_blowup,
    L_kbal,
    blowup_aut_order,
    graph_automorphisms,
)

ROOT = Path(__file__).resolve().parents[1]
GENG = Path("/home/james/src/nauty2_8_9/geng")


def decode_g6(s: str) -> list[list[int]]:
    s = s.strip()
    n = ord(s[0]) - 63
    adj = [[] for _ in range(n)]
    bits: list[int] = []
    for ch in s[1:]:
        v = ord(ch) - 63
        for b in range(5, -1, -1):
            bits.append((v >> b) & 1)
    t = 0
    for j in range(1, n):
        for i in range(j):
            if t < len(bits) and bits[t]:
                adj[i].append(j)
                adj[j].append(i)
            t += 1
    return adj


def has_twins(adj: list[list[int]]) -> bool:
    n = len(adj)
    open_n = [set(adj[i]) for i in range(n)]
    for u in range(n):
        for v in range(u + 1, n):
            if open_n[u] == open_n[v]:
                return True
            nu = (open_n[u] | {u}) - {v}
            nv = (open_n[v] | {v}) - {u}
            if nu == nv:
                return True
    return False


def compositions(n: int, p: int) -> list[tuple[int, ...]]:
    out = []
    for comb in combinations(range(n - 1), p - 1):
        cuts = (0,) + tuple(x + 1 for x in comb) + (n,)
        out.append(tuple(cuts[i + 1] - cuts[i] for i in range(p)))
    return out


def prod_fact(m: tuple[int, ...]) -> int:
    p = 1
    for x in m:
        p *= factorial(x)
    return p


def canon_m(m: tuple[int, ...], auts: list[tuple[int, ...]]) -> tuple[int, ...]:
    p = len(m)
    return min(tuple(m[g[j]] for j in range(p)) for g in auts)


def scan_base(b: int, n: int, T: int, Lk: Fraction) -> dict:
    t0 = time.perf_counter()
    proc = subprocess.run(
        [str(GENG), "-q", str(b)],
        capture_output=True,
        text=True,
        check=True,
    )
    lines = [ln.strip() for ln in proc.stdout.splitlines() if ln.strip()]
    comps = compositions(n, b)
    pf = {m: prod_fact(m) for m in comps}
    n_graphs = 0
    n_twinfree = 0
    n_vecs = 0
    n_eval = 0
    n_beat = 0
    best = None
    records = []
    for gi, g6 in enumerate(lines):
        n_graphs += 1
        adj = decode_g6(g6)
        if has_twins(adj):
            continue
        n_twinfree += 1
        auts = graph_automorphisms(adj)
        A = len(auts)
        seen: set[tuple[int, ...]] = set()
        for m in comps:
            if pf[m] * A < T:
                continue
            cm = canon_m(m, auts)
            if cm in seen:
                continue
            seen.add(cm)
            autG = blowup_aut_order(cm, auts)
            n_vecs += 1
            if autG < T:
                continue
            n_eval += 1
            L = L_blowup(adj, cm, auts)
            ratio = L / Lk
            rec = {
                "base": b,
                "g6": g6,
                "m": list(cm),
                "aut": autG,
                "L": f"{L.numerator}/{L.denominator}",
                "ratio_float": float(ratio),
                "beat": L < Lk,
            }
            if rec["beat"]:
                n_beat += 1
                records.append(rec)
            if best is None or L < best[0]:
                best = (L, rec)
        if (gi + 1) % 200 == 0:
            print(
                f"  b={b} {gi+1}/{len(lines)} tf={n_twinfree} eval={n_eval} "
                f"best={best[1]['ratio_float'] if best else None:.6g}",
                flush=True,
            )
    elapsed = time.perf_counter() - t0
    out = {
        "base": b,
        "n": n,
        "T": T,
        "n_graphs": n_graphs,
        "n_twinfree": n_twinfree,
        "n_part_orbits_above_bound": n_vecs,
        "n_eval": n_eval,
        "n_beat": n_beat,
        "best": None if not best else best[1],
        "beats": records[:20],
        "seconds": elapsed,
    }
    print(
        f"DONE b={b} graphs={n_graphs} twinfree={n_twinfree} eval={n_eval} "
        f"best_ratio={out['best']['ratio_float'] if out['best'] else None} "
        f"beat={n_beat} t={elapsed:.1f}s",
        flush=True,
    )
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=16)
    p.add_argument("--T", type=int, default=1556)
    p.add_argument("--bmin", type=int, default=2)
    p.add_argument("--bmax", type=int, default=7)
    p.add_argument("--out", default="")
    args = p.parse_args()
    Lk = L_kbal(args.n)
    print(f"L(K_bal)={Lk.numerator}/{Lk.denominator} T={args.T}", flush=True)
    results = []
    for b in range(args.bmin, args.bmax + 1):
        results.append(scan_base(b, args.n, args.T, Lk))
    dest = Path(args.out) if args.out else ROOT / "data" / "task4_twin_blowups.json"
    dest.write_text(json.dumps({"n": args.n, "T": args.T, "bases": results}, indent=2) + "\n")
    print("wrote", dest, flush=True)


if __name__ == "__main__":
    main()
