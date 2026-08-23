#!/usr/bin/env python3
"""Exact Aut-threshold T_n = 1 / (L(K_bal) P_n) from bound (1).

Bound (1) of the paper:
    L(G) >= 1 / (|Aut(G)| * P_n)
with P_n = prod_{i=1}^n C(i-1, floor((i-1)/2)).

Hence L(G) < L(K_bal) forces |Aut(G)| > 1 / (L(K_bal) P_n).
The integer threshold is the least integer strictly greater than that
rational (equivalently ceil(T) if T is non-integral, else T+1).
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_ind_path = next(
    p
    for p in [
        ROOT / "code" / "independent.py",
        ROOT / "code" / "independent.py",
        ROOT / "code" / "independent.py",
    ]
    if p.exists()
)
_spec = importlib.util.spec_from_file_location("independent", _ind_path)
_ind = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_ind)


def _attr(*names):
    for n in names:
        if hasattr(_ind, n):
            return getattr(_ind, n)
    raise AttributeError(names)


L_kbal = _attr("L_kbal", "L_kbal")
central_binom_product = _attr("central_binom_product", "central_binom_product")
binom = _ind.binom


def integer_threshold(T: Fraction) -> int:
    if T.denominator == 1:
        return T.numerator + 1
    return T.numerator // T.denominator + 1


def Pn_factors(n: int) -> list[tuple[int, int, int]]:
    """List (i, floor((i-1)/2), C(i-1, floor((i-1)/2))) for i=1..n."""
    rows = []
    for i in range(1, n + 1):
        k = (i - 1) // 2
        rows.append((i, k, binom(i - 1, k)))
    return rows


def main() -> None:
    out: dict = {}
    paper = {13: 183, 14: 2552, 16: 1556, 18: 462}
    for n in (13, 14, 16, 18):
        L = L_kbal(n)
        Pn = central_binom_product(n)
        T = 1 / (L * Pn)
        thresh = integer_threshold(T)
        factors = Pn_factors(n)
        row = {
            "n": n,
            "L_kbal": f"{L.numerator}/{L.denominator}",
            "P_n": str(Pn),
            "P_n_factors": [
                {"i": i, "k": k, "C(i-1,k)": c} for i, k, c in factors
            ],
            "T_exact": f"{T.numerator}/{T.denominator}",
            "T_float": float(T),
            "integer_threshold": thresh,
            "paper": paper[n],
            "matches_paper": thresh == paper[n],
        }
        out[str(n)] = row
        print(
            f"n={n} P_n={Pn}  T={T.numerator}/{T.denominator} "
            f"≈ {float(T):.10f}  |Aut| >= {thresh}  paper={paper[n]}  "
            f"{'OK' if row['matches_paper'] else 'MISMATCH'}"
        )
    path = Path("/home/james/src/graph-likelihood/data/task4_T_exact.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", path)


if __name__ == "__main__":
    main()
