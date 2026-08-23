#!/usr/bin/env python3
"""Evaluate a Holt–Royle / McKay vertex-transitive graph6 census.

Reads one-graph-per-line graph6 (possibly with nauty headers), canonicalises
with labelg, obtains |Aut| from countg -a, evaluates exact L via the subset
recurrence, and uses L(G)=L(complement) to evaluate each complementary pair
once.

Coverage (vertex-transitive case): the Holt–Royle census (Zenodo 4010122,
arXiv:1811.09015) lists every vertex-transitive graph of order n, including
disconnected graphs (k equal copies of a connected VT graph of order n/k).
Every such G has Aut(G) transitive, so this is the complete VT class.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from decimal import Decimal
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
L_subset = _attr("L_subset", "L_subset")
aut_threshold = _attr("aut_threshold", "aut_threshold")


def decode_g6(s: str) -> tuple[int, list[int]]:
    s = s.strip()
    if not s or s[0] in ">#":
        raise ValueError(s)
    n = ord(s[0]) - 63
    adj = [0] * n
    bits: list[int] = []
    for ch in s[1:]:
        v = ord(ch) - 63
        for b in range(5, -1, -1):
            bits.append((v >> b) & 1)
    t = 0
    for j in range(1, n):
        for i in range(j):
            if t < len(bits) and bits[t]:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
            t += 1
    return n, adj


def encode_g6(n: int, adj: list[int]) -> str:
    bits: list[int] = []
    for j in range(1, n):
        for i in range(j):
            bits.append(1 if (adj[i] >> j) & 1 else 0)
    while len(bits) % 6:
        bits.append(0)
    out = [chr(63 + n)]
    for t in range(0, len(bits), 6):
        c = 0
        for b in bits[t : t + 6]:
            c = (c << 1) | b
        out.append(chr(63 + c))
    return "".join(out)


def complement_adj(n: int, adj: list[int]) -> list[int]:
    full = (1 << n) - 1
    return [((~adj[i]) & full) ^ (1 << i) for i in range(n)]


def run_lines(cmd: list[str], lines: list[str]) -> list[str]:
    proc = subprocess.run(
        cmd,
        input="".join(x if x.endswith("\n") else x + "\n" for x in lines),
        capture_output=True,
        text=True,
        check=True,
    )
    return [
        ln.strip()
        for ln in proc.stdout.splitlines()
        if ln.strip() and not ln.startswith(">")
    ]


def batch_canon(g6_lines: list[str], labelg: str) -> list[str]:
    return run_lines([labelg], g6_lines)


def batch_aut(g6_lines: list[str], g6aut: str) -> list[int]:
    """Parse `g6aut` output (`graph6 aut_order`) in input order."""
    proc = subprocess.run(
        [g6aut],
        input="".join(x if x.endswith("\n") else x + "\n" for x in g6_lines),
        capture_output=True,
        text=True,
        check=True,
    )
    auts: list[int] = []
    for line in proc.stdout.splitlines():
        line = line.strip()
        if not line or line.startswith(">"):
            continue
        parts = line.split()
        if len(parts) < 2:
            raise RuntimeError(f"cannot parse g6aut line: {line!r}")
        try:
            auts.append(int(Decimal(parts[-1])))
        except Exception as exc:
            raise RuntimeError(f"cannot parse g6aut line: {line!r}") from exc
    if len(auts) != len(g6_lines):
        raise RuntimeError(
            f"g6aut returned {len(auts)} auts for {len(g6_lines)} graphs"
        )
    return auts


def load_g6_dir(path: Path) -> list[str]:
    raw: list[str] = []
    files = sorted(path.rglob("*")) if path.is_dir() else [path]
    for f in files:
        if not f.is_file():
            continue
        if f.suffix in {".tar", ".gz", ".tgz"}:
            continue
        text = f.read_text(errors="replace")
        for ln in text.splitlines():
            ln = ln.strip()
            if ln and ln[0] not in ">#":
                raw.append(ln)
    return raw


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("source", help="directory of graph6 files, or a single .g6")
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--T", type=int, required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--labelg", default="")
    p.add_argument("--g6aut", default="")
    p.add_argument("--countg", default="", help="unused; kept for compatibility")
    p.add_argument("--provenance", default="")
    args = p.parse_args()

    labelg = args.labelg or shutil.which("labelg")
    g6aut = args.g6aut or shutil.which("g6aut")
    if not g6aut:
        cand = ROOT / "code" / "g6aut"
        if cand.is_file():
            g6aut = str(cand)
    if not labelg or not g6aut:
        sys.exit("need labelg and g6aut")

    src = Path(args.source)
    raw = load_g6_dir(src)
    print(f"raw lines: {len(raw)} from {src}", flush=True)
    t0 = time.perf_counter()

    # Drop graphs of the wrong order.
    filtered = []
    for g6 in raw:
        try:
            n, _ = decode_g6(g6)
        except ValueError:
            continue
        if n == args.n:
            filtered.append(g6)
    print(f"order-n graphs: {len(filtered)}", flush=True)

    canon = batch_canon(filtered, labelg)
    unique = sorted(set(canon))
    print(f"canonical unique: {len(unique)}", flush=True)

    auts = batch_aut(unique, g6aut)
    aut_map = dict(zip(unique, auts))

    # Pair G with complement; evaluate the representative with fewer edges
    # (ties: lexicographically smaller canonical graph6).
    Lk = L_kbal(args.n)
    T_frac = aut_threshold(args.n, Lk)
    print(f"L(K_bal)={Lk.numerator}/{Lk.denominator}", flush=True)
    print(f"T={T_frac.numerator}/{T_frac.denominator} ≈ {float(T_frac):.6f}", flush=True)

    seen: set[str] = set()
    pairs: list[tuple[str, str]] = []
    for g6 in unique:
        if g6 in seen:
            continue
        n, adj = decode_g6(g6)
        cg6 = encode_g6(n, complement_adj(n, adj))
        ccanon = batch_canon([cg6], labelg)[0]
        seen.add(g6)
        seen.add(ccanon)
        e = sum(x.bit_count() for x in adj) // 2
        maxe = n * (n - 1) // 2
        if ccanon == g6:
            pairs.append((g6, g6))  # self-complementary
        elif e < maxe - e or (e == maxe - e and g6 <= ccanon):
            pairs.append((g6, ccanon))
        else:
            pairs.append((ccanon, g6))
    print(f"complement-pairs (incl self-comp): {len(pairs)}", flush=True)

    records = []
    best = None
    for i, (g6, cg6) in enumerate(pairs):
        n, adj = decode_g6(g6)
        aut = aut_map.get(g6)
        if aut is None:
            aut = batch_aut([g6], g6aut)[0]
            aut_map[g6] = aut
        L = L_subset(n, adj, aut)
        ratio = L / Lk
        rec = {
            "id": i,
            "g6": g6,
            "g6_complement": None if cg6 == g6 else cg6,
            "self_complementary": cg6 == g6,
            "aut": aut,
            "edges": sum(x.bit_count() for x in adj) // 2,
            "L": f"{L.numerator}/{L.denominator}",
            "ratio": f"{ratio.numerator}/{ratio.denominator}",
            "ratio_float": float(ratio),
            "beat": L < Lk,
            "aut_ge_T": aut >= args.T,
        }
        records.append(rec)
        if best is None or L < best[0]:
            best = (L, rec)
        if (i + 1) % 25 == 0:
            print(
                f"  {i+1}/{len(pairs)} best_ratio={best[1]['ratio_float']:.6g}",
                flush=True,
            )

    records.sort(key=lambda r: r["ratio_float"])
    elapsed = time.perf_counter() - t0
    n_below_T = sum(1 for r in records if not r["aut_ge_T"])
    n_beat = sum(1 for r in records if r["beat"])
    out = {
        "n": args.n,
        "T": args.T,
        "provenance": args.provenance,
        "raw_lines": len(raw),
        "canonical_unique_labelled": len(unique),
        "complement_pairs_evaluated": len(records),
        "self_complementary": sum(1 for r in records if r["self_complementary"]),
        "aut_ge_T": sum(1 for r in records if r["aut_ge_T"]),
        "aut_lt_T": n_below_T,
        "any_beat": n_beat > 0,
        "n_beat": n_beat,
        "best": None if not records else records[0],
        "seconds": elapsed,
        "graphs": records,
    }
    dest = Path(args.out)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=2) + "\n")
    tsv = dest.with_suffix(".tsv")
    with tsv.open("w") as f:
        f.write("id\tg6\taut\tedges\tL\tratio\tratio_float\tbeat\taut_ge_T\n")
        for r in records:
            f.write(
                f"{r['id']}\t{r['g6']}\t{r['aut']}\t{r['edges']}\t{r['L']}\t"
                f"{r['ratio']}\t{r['ratio_float']:.12g}\t{int(r['beat'])}\t"
                f"{int(r['aut_ge_T'])}\n"
            )
    print(
        f"DONE unique={len(unique)} pairs={len(records)} above_T={out['aut_ge_T']} "
        f"best={out['best']['g6'] if out['best'] else None} "
        f"ratio={out['best']['ratio_float'] if out['best'] else None} "
        f"beat={out['any_beat']} t={elapsed:.1f}s wrote {dest} and {tsv}",
        flush=True,
    )


if __name__ == "__main__":
    main()
