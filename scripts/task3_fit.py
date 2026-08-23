#!/usr/bin/env python3
"""Task 3: e(n) from exact L, least-squares fit, min-cost insertion constant."""
from __future__ import annotations

import json
import math
import re
import sys
from decimal import Decimal, getcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
getcontext().prec = 80
LN2 = Decimal(2).ln()
ENTROPY = Decimal(1) / (Decimal(4) * LN2)  # 1/(4 ln 2)
BIP = Decimal("0.5") - Decimal(1) / (Decimal(8) * LN2)


def log2_frac(num: int, den: int) -> Decimal:
    return (Decimal(num).ln() - Decimal(den).ln()) / LN2


def parse_c_tower(path: Path) -> dict:
    text = path.read_text()
    out = {}
    for m, n, L in re.findall(
        r"m=(\d+)\s+C5 quotient recurrence: n=(\d+).*?\nn=\d+ likelihood=(\d+)/(\d+)",
        text,
        flags=re.S,
    ):
        # fallback simpler
        pass
    for m, Lnum, Lden in re.findall(r"m=(\d+) .*?\nn=\d+ likelihood=(\d+)/(\d+)", text, flags=re.S):
        out[int(m)] = (int(Lnum), int(Lden))
    if not out:
        for n, a, b in re.findall(r"n=(\d+) likelihood=(\d+)/(\d+)", text):
            n = int(n)
            if n % 5 == 0:
                out[n // 5] = (int(a), int(b))
    return out


def e_from_L(n: int, num: int, den: int) -> Decimal:
    # e = -log2(L) / n^2 = log2(den/num) / n^2
    return log2_frac(den, num) / Decimal(n * n)


def lstsq(rows, use_logn_n=True):
    """Fit e = c + a * log2(n)/n + b/n  (or without a)."""
    # Design matrix
    ys = []
    xs = []
    for n, e in rows:
        y = float(e)
        row = [1.0, math.log2(n) / n, 1.0 / n]
        ys.append(y)
        xs.append(row)
    # normal equations
    k = 3
    A = [[0.0] * k for _ in range(k)]
    b = [0.0] * k
    for x, y in zip(xs, ys):
        for i in range(k):
            b[i] += x[i] * y
            for j in range(k):
                A[i][j] += x[i] * x[j]
    # Gaussian elimination
    M = [A[i][:] + [b[i]] for i in range(k)]
    for i in range(k):
        piv = max(range(i, k), key=lambda r: abs(M[r][i]))
        M[i], M[piv] = M[piv], M[i]
        den = M[i][i]
        for j in range(i, k + 1):
            M[i][j] /= den
        for r in range(k):
            if r == i:
                continue
            f = M[r][i]
            for j in range(i, k + 1):
                M[r][j] -= f * M[i][j]
    coef = [M[i][k] for i in range(k)]
    resid = []
    for x, y in zip(xs, ys):
        pred = sum(ci * xi for ci, xi in zip(coef, x))
        resid.append(y - pred)
    rms = math.sqrt(sum(r * r for r in resid) / len(resid))
    return coef, rms, resid


def min_cost_e(m: int) -> float:
    """Min insertion cost W/n^2 for balanced C5 blow-up, float log2."""
    n = 5 * m
    dim = m + 1
    N = dim ** 5

    def idx(a, b, c, d, e):
        return (((a * dim + b) * dim + c) * dim + d) * dim + e

    inf = 1e300
    dp = [inf] * N
    dp[idx(0, 0, 0, 0, 0)] = 0.0
    # iterate by total placed
    coords = list(range(dim))
    for s in range(n):
        # states with sum s
        for a in coords:
            for b in coords:
                for c in coords:
                    for d in coords:
                        e = s - a - b - c - d
                        if e < 0 or e > m:
                            continue
                        cur = dp[idx(a, b, c, d, e)]
                        if cur >= inf:
                            continue
                        counts = (a, b, c, d, e)
                        for i in range(5):
                            if counts[i] >= m:
                                continue
                            deg = counts[(i - 1) % 5] + counts[(i + 1) % 5]
                            # log2 C(s, deg)
                            if deg < 0 or deg > s:
                                add = inf
                            elif deg == 0 or deg == s:
                                add = 0.0
                            else:
                                add = (
                                    math.lgamma(s + 1)
                                    - math.lgamma(deg + 1)
                                    - math.lgamma(s - deg + 1)
                                ) / math.log(2.0)
                            nxt = list(counts)
                            nxt[i] += 1
                            j = idx(*nxt)
                            val = cur + add
                            if val < dp[j]:
                                dp[j] = val
    W = dp[idx(m, m, m, m, m)]
    # e ~ (W + log2(n! |Aut|))/n^2  wait no: A >= Nseq * 2^{-W}, L = A/(n! Aut)
    # conservative: e_term = W/n^2 is the leading piece of -log2 of one term
    # The true -log2 L = log2(n! Aut) - log2 A <= log2(n! Aut) + W
    # and >= W + log2(n! Aut) - log2(#terms) >= W + O(n log n)
    return W / (n * n), W


def main():
    py_map = {}
    for fname in ("task2_py_tower_m1_12.json", "task2_py_tower_m13_16.json"):
        p = ROOT / "data" / fname
        if p.exists():
            rows_js = json.loads(p.read_text())
            if isinstance(rows_js, dict) and "rows" in rows_js:
                rows_js = rows_js["rows"]
            for row in rows_js:
                py_map[row["m"]] = row

    # parse both C tower logs
    c_map = {}
    for fname in ("task2_C_tower_m1_12.txt", "task2_C_tower_m13_16.txt"):
        p = ROOT / "data" / fname
        if p.exists():
            c_map.update(parse_c_tower(p))

    rows = []
    print(f"{'m':>3} {'n':>4} {'e_py':>14} {'e_C':>14} {'matchL':>6}")
    for m in range(1, 17):
        n = 5 * m
        e_py = e_C = None
        match = ""
        if m in py_map:
            num, den = map(int, py_map[m]["L"].split("/"))
            e_py = e_from_L(n, num, den)
        if m in c_map:
            e_C = e_from_L(n, *c_map[m])
            if m in py_map:
                pnum, pden = map(int, py_map[m]["L"].split("/"))
                match = "OK" if (pnum, pden) == c_map[m] else "NO"
        e_use = e_py or e_C
        rows.append((n, m, e_use))
        print(
            f"{m:3d} {n:4d} {str(e_py)[:14] if e_py else '-':>14} "
            f"{str(e_C)[:14] if e_C else '-':>14} {match:>6}"
        )

    print("\nentropy 1/(4 ln 2) =", format(ENTROPY, ".12f"))
    print("bipartite 1/2-1/(8 ln 2) =", format(BIP, ".12f"))

    # fits on largest terms, odd/even m separately; store residuals
    fits = {}
    for label, pred in (
        ("odd m>=5", lambda m: m >= 5 and m % 2 == 1),
        ("even m>=6", lambda m: m >= 6 and m % 2 == 0),
        ("all m>=8", lambda m: m >= 8),
        ("all m>=10", lambda m: m >= 10),
        ("odd all m>=1", lambda m: m % 2 == 1),
        ("even all m>=2", lambda m: m % 2 == 0),
    ):
        data = [(n, float(e)) for n, m, e in rows if e is not None and pred(m)]
        if len(data) < 3:
            continue
        coef, rms, resid = lstsq(data)
        c, a, b = coef
        print(
            f"fit {label}: c={c:.10f} a={a:.6f} b={b:.6f} rms={rms:.3e} "
            f"c-entropy={c - float(ENTROPY):+.6f} N={len(data)}"
        )
        fits[label] = {
            "c": c,
            "a": a,
            "b": b,
            "rms": rms,
            "N": len(data),
            "c_minus_entropy": c - float(ENTROPY),
            "points": [
                {
                    "n": n,
                    "e": y,
                    "predicted": c + a * math.log2(n) / n + b / n,
                    "residual": r,
                }
                for (n, y), r in zip(data, resid)
            ],
        }

    # min-cost insertion
    print("\nmin-cost insertion W/n^2 (float):")
    mc = {}
    for m in (4, 6, 8, 10, 12):
        wn2, W = min_cost_e(m)
        mc[m] = wn2
        print(f"  m={m} n={5*m} W/n^2={wn2:.10f}  vs e_exact={float(rows[m-1][2]):.10f}")

    out = {
        "entropy": format(ENTROPY, ".20f"),
        "bipartite": format(BIP, ".20f"),
        "e_rows": [
            {"m": m, "n": n, "e": format(e, ".20f") if e is not None else None}
            for n, m, e in rows
        ],
        "min_cost_W_over_n2": mc,
        "fits": fits,
        "parity_e_rows": {
            "odd_m": [
                {"m": m, "n": n, "e": format(e, ".20f")}
                for n, m, e in rows
                if e is not None and m % 2 == 1
            ],
            "even_m": [
                {"m": m, "n": n, "e": format(e, ".20f")}
                for n, m, e in rows
                if e is not None and m % 2 == 0
            ],
        },
        "honest_framing": (
            "Fitted c lies in about 0.363-0.366 for the windows used, "
            "all slightly above 1/(4 ln 2) ≈ 0.36067. This is a numerical "
            "indication, not a proof that lim e(n) exceeds the entropy constant."
        ),
    }
    dest = ROOT / "data" / "task3_fit.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest)


if __name__ == "__main__":
    sys.exit(main())
