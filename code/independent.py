#!/usr/bin/env python3
"""Independent exact graph-likelihood evaluator for arXiv:2608.19467.

Written from the paper's recurrences, not from likelihood.c.

    A(K1) = 1
    A(G)  = sum_v A(G-v) / C(|V|-1, deg v)                           (Prop. 8)
    A(H[m]) = sum_{i: m_i>0} m_i A(H[m-e_i]) / C(n-1, sum_{j~i} m_j) (Prop. 11)
    L(G)  = A(G) / (n! |Aut(G)|)

    |Aut(C5[m])| = (prod m_i!) * |{g in D5 : m o g = m}|  (all m_i >= 1)

    b(n,0)=b(n,n)=1
    b(n,a) = b(n-1,a-1)/C(n-1,a-1) + b(n-1,a)/C(n-1,a)       (recurrence 6)
    L(K_{a,n-a}) = b(n,a) / (n! * (2 if 2a==n else 1))

Integer scaling: L_r = lcm_d C(r-1,d) (L_1=1), P[k]=prod_{r<=k} L_r,
Fhat = A * P[n] is an integer.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from decimal import Decimal, getcontext
from fractions import Fraction
from math import factorial, gcd, lcm
from typing import Dict, List, Optional, Sequence, Tuple

Parts = Tuple[int, ...]


def binom(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    k = min(k, n - k)
    r = 1
    for i in range(1, k + 1):
        r = r * (n - k + i) // i
    return r


def lcm_binomials(r: int) -> int:
    if r <= 1:
        return 1
    acc = 1
    nm1 = r - 1
    for d in range(nm1 + 1):
        acc = lcm(acc, binom(nm1, d))
    return acc


def reduce_frac(num: int, den: int) -> Tuple[int, int]:
    if den < 0:
        num, den = -num, -den
    g = gcd(num, den)
    return num // g, den // g


def d5_images(m: Sequence[int]) -> List[Parts]:
    m = tuple(m)
    out = []
    for i in range(5):
        for s in (1, -1):
            out.append(tuple(m[(i + s * j) % 5] for j in range(5)))
    return out


def d5_canonical(m: Sequence[int]) -> Parts:
    return min(d5_images(m))


def c5_aut_order(m: Sequence[int]) -> int:
    fact = 1
    for x in m:
        if x < 1:
            raise ValueError("c5_aut_order requires positive parts")
        fact *= factorial(x)
    stab = sum(1 for img in d5_images(m) if img == tuple(m))
    return fact * stab


def positive_c5_orbits(n: int) -> List[Parts]:
    reps = []
    for a in range(1, n - 3):
        for b in range(1, n - 2 - a):
            for c in range(1, n - 1 - a - b):
                for d in range(1, n - a - b - c):
                    e = n - a - b - c - d
                    if e < 1:
                        continue
                    m = (a, b, c, d, e)
                    if d5_canonical(m) == m:
                        reps.append(m)
    return reps


def most_balanced_c5(n: int) -> Parts:
    q, r = divmod(n, 5)
    raw = [q + 1] * r + [q] * (5 - r)
    return d5_canonical(raw)


def c5_deg(m: Sequence[int], i: int) -> int:
    return m[(i + 4) % 5] + m[(i + 1) % 5]


def c5_fhat(m: Sequence[int]) -> Tuple[int, int]:
    """Return (Fhat, P[n]) for a C5 blow-up (zeros allowed). A = Fhat / P[n]."""
    target = tuple(int(x) for x in m)
    n = sum(target)
    if n == 0:
        return 1, 1
    prev: Dict[Parts, int] = {(0, 0, 0, 0, 0): 1}
    P = 1
    for k in range(1, n + 1):
        Lk = lcm_binomials(k)
        P *= Lk
        cur: Dict[Parts, int] = {}

        def rec(i: int, left: int, acc: List[int]) -> None:
            if i == 4:
                acc[4] = left
                if left > target[4]:
                    return
                st = tuple(acc)
                total = 0
                for j in range(5):
                    if st[j] == 0:
                        continue
                    pred_l = list(st)
                    pred_l[j] -= 1
                    pred = tuple(pred_l)
                    ckd = binom(k - 1, c5_deg(st, j))
                    total += st[j] * prev[pred] * (Lk // ckd)
                cur[st] = total
                return
            hi = min(target[i], left)
            for v in range(hi + 1):
                acc[i] = v
                rec(i + 1, left - v, acc)

        rec(0, k, [0, 0, 0, 0, 0])
        prev = cur
    return prev[target], P


def L_c5(m: Sequence[int]) -> Fraction:
    m = tuple(m)
    n = sum(m)
    fhat, P = c5_fhat(m)
    aut = c5_aut_order(m)
    num, den = reduce_frac(fhat, P * factorial(n) * aut)
    return Fraction(num, den)


def b_row(n: int) -> List[Fraction]:
    cur = [Fraction(1)]
    for t in range(1, n + 1):
        nxt = [Fraction(0)] * (t + 1)
        nxt[0] = Fraction(1)
        nxt[t] = Fraction(1)
        for a in range(1, t):
            s = cur[a - 1] / binom(t - 1, a - 1)
            s += cur[a] / binom(t - 1, a)
            nxt[a] = s
        cur = nxt
    return cur


def L_bipartite(a: int, b: int) -> Fraction:
    n = a + b
    extra = 2 if a == b else 1
    return b_row(n)[a] / (factorial(n) * extra)


def L_kbal(n: int) -> Fraction:
    return L_bipartite(n // 2, n - n // 2)


def L_subset(n: int, adj_bits: Sequence[int], aut: int) -> Fraction:
    count = 1 << n
    psc = [0] * (n + 1)
    psc[1] = 1
    for r in range(2, n + 1):
        psc[r] = psc[r - 1] * lcm_binomials(r)
    tab = [0] * count
    for mask in range(1, count):
        k = mask.bit_count()
        if k == 1:
            tab[mask] = 1
            continue
        Lk = psc[k] // psc[k - 1]
        total = 0
        t = mask
        while t:
            v = (t & -t).bit_length() - 1
            sub = mask ^ (1 << v)
            d = (adj_bits[v] & sub).bit_count()
            total += tab[sub] * (Lk // binom(k - 1, d))
            t &= t - 1
        tab[mask] = total
    num, den = reduce_frac(tab[count - 1], psc[n] * factorial(n) * aut)
    return Fraction(num, den)


def k_bipartite_adj(a: int, b: int) -> List[int]:
    n = a + b
    adj = [0] * n
    for i in range(a):
        for j in range(a, n):
            adj[i] |= 1 << j
            adj[j] |= 1 << i
    return adj


def c5_family_scan(n: int) -> dict:
    t0 = time.perf_counter()
    orbits = positive_c5_orbits(n)
    prev: Dict[Parts, int] = {(0, 0, 0, 0, 0): 1}
    P = 1
    for k in range(1, n + 1):
        Lk = lcm_binomials(k)
        P *= Lk
        cur: Dict[Parts, int] = {}
        for a in range(k + 1):
            for b in range(k - a + 1):
                for c in range(k - a - b + 1):
                    for d in range(k - a - b - c + 1):
                        e = k - a - b - c - d
                        m = (a, b, c, d, e)
                        total = 0
                        for j, mj in enumerate(m):
                            if mj == 0:
                                continue
                            pred_l = list(m)
                            pred_l[j] -= 1
                            pred = tuple(pred_l)
                            ckd = binom(k - 1, c5_deg(m, j))
                            total += mj * prev[pred] * (Lk // ckd)
                        cur[m] = total
        prev = cur
    nfac = factorial(n)
    den_base = P * nfac
    best_L: Optional[Fraction] = None
    best_parts: Optional[List[int]] = None
    for m in orbits:
        fhat = prev[m]
        aut = c5_aut_order(m)
        num, den = reduce_frac(fhat, den_base * aut)
        L = Fraction(num, den)
        if best_L is None or L < best_L or (L == best_L and list(m) < best_parts):  # type: ignore[operator]
            best_L = L
            best_parts = list(m)
    elapsed = time.perf_counter() - t0
    assert best_L is not None and best_parts is not None
    return {
        "n": n,
        "n_orbits": len(orbits),
        "best_parts": best_parts,
        "best_L": f"{best_L.numerator}/{best_L.denominator}",
        "seconds": elapsed,
    }


def graph_automorphisms(adj: Sequence[Sequence[int]]) -> List[Tuple[int, ...]]:
    n = len(adj)
    neigh = [frozenset(adj[i]) for i in range(n)]
    deg = [len(neigh[i]) for i in range(n)]
    autos: List[Tuple[int, ...]] = []

    def rec(img: List[int], used: List[bool]) -> None:
        k = len(img)
        if k == n:
            autos.append(tuple(img))
            return
        for v in range(n):
            if used[v] or deg[v] != deg[k]:
                continue
            ok = True
            for u in range(k):
                if (u in neigh[k]) != (img[u] in neigh[v]):
                    ok = False
                    break
            if not ok:
                continue
            used[v] = True
            img.append(v)
            rec(img, used)
            img.pop()
            used[v] = False

    rec([], [False] * n)
    return autos


def blowup_aut_order(m: Sequence[int], aut_H: Sequence[Sequence[int]]) -> int:
    fact = 1
    for x in m:
        fact *= factorial(x)
    mt = tuple(m)
    p = len(m)
    stab = sum(1 for g in aut_H if all(mt[g[j]] == mt[j] for j in range(p)))
    return fact * stab


def L_blowup(adj: Sequence[Sequence[int]], m: Sequence[int], aut_H=None) -> Fraction:
    m = tuple(int(x) for x in m)
    p = len(m)
    n = sum(m)
    if aut_H is None:
        aut_H = graph_automorphisms(adj)
    prev: Dict[Parts, int] = {(0,) * p: 1}
    P = 1
    neigh = [tuple(adj[i]) for i in range(p)]
    for k in range(1, n + 1):
        Lk = lcm_binomials(k)
        P *= Lk
        cur: Dict[Parts, int] = {}

        def rec(i: int, left: int, acc: List[int]) -> None:
            if i == p - 1:
                acc[p - 1] = left
                if left > m[p - 1]:
                    return
                st = tuple(acc)
                total = 0
                for j in range(p):
                    if st[j] == 0:
                        continue
                    pred_l = list(st)
                    pred_l[j] -= 1
                    pred = tuple(pred_l)
                    ddeg = sum(st[u] for u in neigh[j])
                    ckd = binom(k - 1, ddeg)
                    total += st[j] * prev[pred] * (Lk // ckd)
                cur[st] = total
                return
            hi = min(m[i], left)
            for v in range(hi + 1):
                acc[i] = v
                rec(i + 1, left - v, acc)

        rec(0, k, [0] * p)
        prev = cur
    fhat = prev[m]
    aut = blowup_aut_order(m, aut_H)
    num, den = reduce_frac(fhat, P * factorial(n) * aut)
    return Fraction(num, den)


def cycle_adj(p: int) -> List[List[int]]:
    return [[(i - 1) % p, (i + 1) % p] for i in range(p)]


def petersen_adj() -> List[List[int]]:
    adj: List[List[int]] = [[] for _ in range(10)]
    for i in range(5):
        adj[i] += [(i + 1) % 5, (i + 4) % 5, i + 5]
        adj[i + 5] += [i, (i + 2) % 5 + 5, (i + 3) % 5 + 5]
    return [sorted(set(row)) for row in adj]


def central_binom_product(n: int) -> int:
    acc = 1
    for i in range(1, n + 1):
        acc *= binom(i - 1, (i - 1) // 2)
    return acc


def aut_threshold(n: int, L_star: Optional[Fraction] = None) -> Fraction:
    if L_star is None:
        L_star = L_kbal(n)
    return 1 / (L_star * central_binom_product(n))


def log2_fraction(num: int, den: int, prec: int = 80) -> str:
    getcontext().prec = prec
    val = (Decimal(num).ln() - Decimal(den).ln()) / Decimal(2).ln()
    return format(val, ".30f")


# Exact checkpoint strings from the task brief / Theorem 5.
CP = {
    "c5_1": Fraction(1, 270),
    "c5_2": Fraction(1217597, 145179288576000000),
    "c5_3": Fraction(63977511069907, 43503039787261205233506826321920000000000),
    "c5_4": Fraction(
        18253737912596545805697739,
        172081388803736230731527560120521881526959111746593503576064000000000000000,
    ),
    "c5_5": Fraction(
        89836876981426998730761750566013583995604129,
        7472886125258401659424666180865536229365989642933541083625988665391638703650134317041815225253232640000000000000000000000,
    ),
    "k78": Fraction(7628328998218493, 1044072954894268925604163831726080000000000),
    "n16": Fraction(6360639670045255013, 2149955028718278571604094162290343936000000000000),
    "n18": Fraction(
        2508406335896865033959107,
        161372813408313544688091786898512049839618706636800000000000000,
    ),
}


def cmd_check(_args: argparse.Namespace) -> int:
    failed = 0

    def report(name: str, got: Fraction, expect: Fraction) -> None:
        nonlocal failed
        ok = got == expect
        if not ok:
            failed += 1
        print(f"{'PASS' if ok else 'FAIL'} {name}")
        if not ok:
            print(f"  got      {got}")
            print(f"  expected {expect}")

    t0 = time.perf_counter()
    report("L(C5[1])", L_c5((1, 1, 1, 1, 1)), CP["c5_1"])
    report("L(C5[2])", L_c5((2, 2, 2, 2, 2)), CP["c5_2"])
    report("L(C5[3]) Thm5", L_c5((3, 3, 3, 3, 3)), CP["c5_3"])
    report("L(C5[4])", L_c5((4, 4, 4, 4, 4)), CP["c5_4"])
    report("L(C5[5])", L_c5((5, 5, 5, 5, 5)), CP["c5_5"])
    report("L(C5[3,3,3,3,4])", L_c5((3, 3, 3, 3, 4)), CP["n16"])
    report("L(C5[3,3,4,4,4])", L_c5((3, 3, 4, 4, 4)), CP["n18"])
    report("L(K7,8) rec6", L_bipartite(7, 8), CP["k78"])
    report(
        "L(K7,8) subset",
        L_subset(15, k_bipartite_adj(7, 8), factorial(7) * factorial(8)),
        CP["k78"],
    )
    c5_bits = [0] * 5
    for i in range(5):
        c5_bits[i] |= 1 << ((i + 1) % 5)
        c5_bits[i] |= 1 << ((i + 4) % 5)
    report("L(C5) subset", L_subset(5, c5_bits, 10), CP["c5_1"])
    report("L(C5[3]) Prop11 general", L_blowup(cycle_adj(5), (3, 3, 3, 3, 3)), CP["c5_3"])
    print(f"independent checkpoints in {time.perf_counter() - t0:.3f}s")
    return 1 if failed else 0


def cmd_c5(args: argparse.Namespace) -> int:
    m = tuple(args.parts)
    t0 = time.perf_counter()
    L = L_c5(m)
    dt = time.perf_counter() - t0
    print(f"n={sum(m)} likelihood={L.numerator}/{L.denominator}")
    print(f"# aut={c5_aut_order(m)} seconds={dt:.6f}", file=sys.stderr)
    return 0


def cmd_kbal(args: argparse.Namespace) -> int:
    n = args.n
    t0 = time.perf_counter()
    L = L_kbal(n)
    dt = time.perf_counter() - t0
    print(f"n={n} likelihood={L.numerator}/{L.denominator}")
    print(f"# seconds={dt:.6f}", file=sys.stderr)
    return 0


def cmd_table(args: argparse.Namespace) -> int:
    out = []
    for n in range(args.nmin, args.nmax + 1):
        scan = c5_family_scan(n)
        Lk = L_kbal(n)
        best = Fraction(scan["best_L"])
        ratio = best / Lk
        beat = best < Lk
        bal = list(most_balanced_c5(n))
        row = {
            "n": n,
            "best_parts": scan["best_parts"],
            "balanced_parts": bal,
            "n_orbits": scan["n_orbits"],
            "L_c5": scan["best_L"],
            "L_kbal": f"{Lk.numerator}/{Lk.denominator}",
            "beat": beat,
            "ratio_exact": f"{ratio.numerator}/{ratio.denominator}",
            "seconds_c5_scan": scan["seconds"],
            "best_is_balanced": scan["best_parts"] == bal,
        }
        out.append(row)
        flag = "WITNESS" if beat else "GAP"
        print(
            f"n={n:3d} {flag:7s} parts={scan['best_parts']} orbits={scan['n_orbits']:5d} "
            f"ratio={float(ratio):.8g} balanced={row['best_is_balanced']} "
            f"t={scan['seconds']:.2f}s"
        )
        sys.stdout.flush()
    if args.out:
        with open(args.out, "w") as f:
            json.dump(out, f, indent=2)
            f.write("\n")
    return 0


def cmd_tower(args: argparse.Namespace) -> int:
    rows = []
    getcontext().prec = 80
    for m in range(args.mmin, args.mmax + 1):
        t0 = time.perf_counter()
        L = L_c5((m, m, m, m, m))
        dt = time.perf_counter() - t0
        n = 5 * m
        elog = log2_fraction(L.denominator, L.numerator)
        e = Decimal(elog) / Decimal(n * n)
        row = {
            "m": m,
            "n": n,
            "L": f"{L.numerator}/{L.denominator}",
            "seconds": dt,
            "-log2_L": elog,
            "e(n)": format(e, ".20f"),
        }
        rows.append(row)
        print(f"m={m:2d} n={n:3d} e={format(e, '.12f')} t={dt:.2f}s")
        sys.stdout.flush()
    if args.out:
        with open(args.out, "w") as f:
            json.dump(rows, f, indent=2)
            f.write("\n")
    return 0


def cmd_threshold(args: argparse.Namespace) -> int:
    for n in args.ns:
        L = L_kbal(n)
        T = aut_threshold(n, L)
        if T.denominator == 1:
            thresh = T.numerator + 1
        else:
            thresh = T.numerator // T.denominator + 1
        print(f"n={n} L_kbal={L.numerator}/{L.denominator}")
        print(f"  T=1/(L P_n)={T.numerator}/{T.denominator}")
        print(f"  necessary |Aut| >= {thresh} to possibly beat K_bal")
        print(f"  T approx {float(T):.10f}")
    return 0


def cmd_blowup(args: argparse.Namespace) -> int:
    if args.base == "c5":
        adj = cycle_adj(5)
    elif args.base == "c7":
        adj = cycle_adj(7)
    elif args.base == "c9":
        adj = cycle_adj(9)
    elif args.base == "petersen":
        adj = petersen_adj()
    else:
        print("unknown base", file=sys.stderr)
        return 2
    m = tuple(args.parts)
    if len(m) != len(adj):
        print(f"expected {len(adj)} parts, got {len(m)}", file=sys.stderr)
        return 2
    t0 = time.perf_counter()
    L = L_blowup(adj, m)
    dt = time.perf_counter() - t0
    print(f"n={sum(m)} likelihood={L.numerator}/{L.denominator}")
    print(f"# seconds={dt:.6f}", file=sys.stderr)
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("check")
    s.set_defaults(func=cmd_check)

    s = sub.add_parser("c5")
    s.add_argument("parts", nargs=5, type=int)
    s.set_defaults(func=cmd_c5)

    s = sub.add_parser("kbal")
    s.add_argument("n", type=int)
    s.set_defaults(func=cmd_kbal)

    s = sub.add_parser("table")
    s.add_argument("--nmin", type=int, default=15)
    s.add_argument("--nmax", type=int, default=40)
    s.add_argument("--out", default="")
    s.set_defaults(func=cmd_table)

    s = sub.add_parser("tower")
    s.add_argument("--mmin", type=int, default=1)
    s.add_argument("--mmax", type=int, default=16)
    s.add_argument("--out", default="")
    s.set_defaults(func=cmd_tower)

    s = sub.add_parser("threshold")
    s.add_argument("ns", nargs="+", type=int)
    s.set_defaults(func=cmd_threshold)

    s = sub.add_parser("blowup")
    s.add_argument("base", choices=["c5", "c7", "c9", "petersen"])
    s.add_argument("parts", nargs="+", type=int)
    s.set_defaults(func=cmd_blowup)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
