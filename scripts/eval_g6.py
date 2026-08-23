#!/usr/bin/env python3
"""Evaluate graph6 graphs of order n<=20 with exact subset likelihood.

If a nauty `countg` binary is available it is used for |Aut| and canonical
graph6; otherwise |Aut| may be passed as a second field.

Compares to L(K_bal). Writes a JSON table of unique (canonical) graphs
with |Aut| >= T.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from independent import L_kbal, L_subset, aut_threshold  # noqa: E402


def decode_g6(s: str) -> tuple[int, list[int]]:
    s = s.strip()
    if not s or s[0] == ">":
        raise ValueError(s)
    n = ord(s[0]) - 63
    adj = [0] * n
    bits = []
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


def complement_g6_bits(n: int, adj: list[int]) -> list[int]:
    full = (1 << n) - 1
    return [((~adj[i]) & full) ^ (1 << i) for i in range(n)]


def nauty_canon_aut(g6_lines: list[str], countg: str) -> list[tuple[str, int]]:
    """Return (canonical_g6, aut_order) for each input line, via countg -a."""
    proc = subprocess.run(
        [countg, "-a"],
        input="".join(x if x.endswith("\n") else x + "\n" for x in g6_lines),
        capture_output=True,
        text=True,
        check=True,
    )
    # countg -a prints: graph6  aut_group_size  (format depends on version)
    out = []
    for line in proc.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) == 1:
            # fallback: just g6
            continue
        # typical: "Og~~~~wW@~?  384" or "Og..  Aut=384"
        g6 = parts[0]
        aut_s = parts[-1]
        aut_s = aut_s.replace("Aut=", "").replace(";", "")
        out.append((g6, int(aut_s)))
    return out


def labelg_canon(g6_lines: list[str], labelg: str) -> list[str]:
    proc = subprocess.run(
        [labelg],
        input="".join(x if x.endswith("\n") else x + "\n" for x in g6_lines),
        capture_output=True,
        text=True,
        check=True,
    )
    return [ln.strip() for ln in proc.stdout.splitlines() if ln.strip() and not ln.startswith(">")]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("g6file")
    p.add_argument("--n", type=int, default=16)
    p.add_argument("--T", type=int, default=1556)
    p.add_argument("--out", default="")
    p.add_argument("--countg", default="")
    p.add_argument("--labelg", default="")
    args = p.parse_args()

    raw = [
        ln.strip()
        for ln in Path(args.g6file).read_text().splitlines()
        if ln.strip() and not ln.startswith("#")
    ]
    print(f"raw graphs: {len(raw)}", flush=True)

    labelg = args.labelg or shutil.which("labelg")
    countg = args.countg or shutil.which("countg")
    t0 = time.perf_counter()

    # Canonicalise to unique iso types.
    if labelg:
        print(f"canonicalising with {labelg}", flush=True)
        canon_lines = labelg_canon(raw, labelg)
        unique = sorted(set(canon_lines))
        # keep the representative with fewer 1-bits vs complement
        kept = []
        seen = set()
        for g6 in unique:
            n, adj = decode_g6(g6)
            e = sum(x.bit_count() for x in adj) // 2
            maxe = n * (n - 1) // 2
            if e > maxe - e:
                # prefer complement representative: re-encode later via labelg
                cadj = complement_g6_bits(n, adj)
                # use original if already processed complement
                key = g6
            else:
                key = g6
            if key in seen:
                continue
            seen.add(key)
            kept.append(g6)
        unique = kept
        print(f"unique labelled: {len(unique)}", flush=True)
    else:
        unique = sorted(set(raw))
        print("WARNING: no labelg; uniqueness is labelling-dependent", flush=True)

    Lk = L_kbal(args.n)
    T_frac = aut_threshold(args.n, Lk)
    print(f"L(K_bal)={Lk.numerator}/{Lk.denominator}", flush=True)
    print(f"T={T_frac.numerator}/{T_frac.denominator} ≈ {float(T_frac):.6f}", flush=True)

    records = []
    best = None
    for i, g6 in enumerate(unique):
        n, adj = decode_g6(g6)
        if n != args.n:
            continue
        # |Aut| via countg on a single graph if available
        aut = None
        if countg:
            try:
                r = subprocess.run(
                    [countg, "-a"],
                    input=g6 + "\n",
                    capture_output=True,
                    text=True,
                    check=True,
                )
                toks = r.stdout.split()
                for tok in reversed(toks):
                    tok = tok.replace("Aut=", "").replace(";", "")
                    if tok.isdigit():
                        aut = int(tok)
                        break
            except subprocess.CalledProcessError:
                aut = None
        if aut is None:
            print(f"skip (no aut) {g6}", flush=True)
            continue
        if aut < args.T:
            continue
        L = L_subset(n, adj, aut)
        ratio = L / Lk
        rec = {
            "g6": g6,
            "aut": aut,
            "edges": sum(x.bit_count() for x in adj) // 2,
            "L": f"{L.numerator}/{L.denominator}",
            "ratio": f"{ratio.numerator}/{ratio.denominator}",
            "ratio_float": float(ratio),
            "beat": L < Lk,
        }
        records.append(rec)
        if best is None or L < best[0]:
            best = (L, rec)
        if (i + 1) % 20 == 0:
            print(
                f"  {i+1}/{len(unique)} evaluated_above_T={len(records)} "
                f"best_ratio={best[1]['ratio_float']:.6g}",
                flush=True,
            )

    records.sort(key=lambda r: r["ratio_float"])
    elapsed = time.perf_counter() - t0
    out = {
        "n": args.n,
        "T": args.T,
        "raw": len(raw),
        "unique": len(unique),
        "evaluated_aut_ge_T": len(records),
        "best": None if not records else records[0],
        "any_beat": any(r["beat"] for r in records),
        "seconds": elapsed,
        "graphs": records,
    }
    dest = Path(args.out) if args.out else Path(args.g6file).with_suffix(".eval.json")
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print(
        f"DONE unique={len(unique)} above_T={len(records)} "
        f"best={out['best']['g6'] if out['best'] else None} "
        f"ratio={out['best']['ratio_float'] if out['best'] else None} "
        f"beat={out['any_beat']} t={elapsed:.1f}s wrote {dest}",
        flush=True,
    )


if __name__ == "__main__":
    main()
