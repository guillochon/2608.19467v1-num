# Exact graph-likelihood computations for arXiv:2608.19467

Computations for Severini–Weisstein, *The minimum of the graph likelihood*.
All reported likelihoods are exact rationals. Family minima are **upper bounds**
on \(\min_G L(G)\) and **witnesses against Conjecture 1** when they beat
\(K_{\lfloor n/2\rfloor,\lceil n/2\rceil}\). Nothing in this file is a claim of
a global minimum for \(n\ge 15\).

## Implementations (independent)

| Tag | Path | Method |
|---|---|---|
| **C** | `likelihood.c` | C11/GMP. C5: full quotient grid of `mpq`. Bipartite: subset recurrence (n≤20) and recurrence (6) (`kbal`). Explicit graphs: subset + graph6. |
| **Py** | `code/independent.py` | Python ints. C5: integer-scaled layer DP (Prop. 11) with \(L_r=\mathrm{lcm}_d\binom{r-1}{d}\). Bipartite: recurrence (6). Explicit graphs: integer-scaled subset. Blow-ups of a fixed base: Prop. 11. |

Every headline number below was checked by both, except: (i) C5 family *scans*
(all D5-orbits at a given n) are done in Py, with C re-evaluating the reported
winner and \(K_{\mathrm{bal}}\); (ii) C7/C9/Petersen scans are Prop. 11 in Py
only, with C5/K_bal checkpoints tying the same code to C.

Sources hashed after the computational run (see `logs/hardware.txt` for CPU/RAM):

```
likelihood.c            sha256  (see git)
code/independent.py     sha256  (see git)
```

Hardware: WSL2 Linux, 12th Gen Intel Core i5-12600K, 16 logical CPUs, 16 GB RAM,
ext4 on `/home/james/src/graph-likelihood`. Git: `1b4b3d6` (Task 0) plus later
commits on `master`.

---

## Task 0 — Checkpoints

Both implementations reproduce every supplied checkpoint.

| Object | Exact \(L\) | C | Py |
|---|---|---|---|
| \(C_5[1^5]\) | \(1/270\) | yes | yes |
| \(C_5[2^5]\) | \(1217597/145179288576000000\) | yes | yes |
| \(C_5[3^5]\) (Thm 5) | \(63977511069907/43503039787261205233506826321920000000000\) | yes | yes |
| \(K_{7,8}\) (Thm 5) | \(7628328998218493/1044072954894268925604163831726080000000000\) | yes | yes |
| \(C_5[4^5]\) | \(18253737912596545805697739/172081388803736230731527560120521881526959111746593503576064000000000000000\) | yes | yes |
| \(C_5[5^5]\) | \(89836876981426998730761750566013583995604129/7472886125258401659424666180865536229365989642933541083625988665391638703650134317041815225253232640000000000000000000000\) | yes | yes |
| \(C_5[3,3,3,3,4]\) | \(6360639670045255013/2149955028718278571604094162290343936000000000000\) | yes | yes |
| \(C_5[3,3,4,4,4]\) | \(2508406335896865033959107/161372813408313544688091786898512049839618706636800000000000000\) | yes | yes |

Additionally \(L(K_{8,8})\) and \(L(K_{9,9})\) from recurrence (6) match the
n=16 and n=18 subset recurrences in C (2^16 and 2^18 states).

`make check` (C, Thm 5 pair) and `python3 code/independent.py check` both pass.

---

## Result A — C5 family vs \(K_{\mathrm{bal}}\), n=15..50

For each n, Py evaluates **every** D5-orbit of positive 5-tuples summing to n
(no part-size cap) by one shared layer DP; C then recomputes \(L\) of the
reported minimiser and of \(K_{\lfloor n/2\rfloor,\lceil n/2\rceil}\). All
C confirmations matched.

**Headline.** Witnesses against Conjecture 1 exist in the C5 family at
n=15, 17, and every 19≤n≤50. There is **no** C5-family witness at n=16 or n=18.

The minimiser is always a 5-tuple with parts in \(\{\lfloor n/5\rfloor,\lceil n/5\rceil\}\).
When n=5k+3 the two smaller parts are placed **non-adjacent** on the 5-cycle
(alternating), not clumped: e.g. n=23 is (4,5,4,5,5) rather than (4,4,5,5,5).
That is still the most balanced composition, up to D5, once adjacency on C5 is
taken into account.

| n | best parts | beat \(K_{\mathrm{bal}}\)? | \(L(C_5[m])/L(K_{\mathrm{bal}})\) (float) | C+Py |
|---|---|---|---|---|
| 15 | (3,3,3,3,3) | WITNESS | 0.2012839596 | yes |
| 16 | (3,3,3,3,4) | GAP | 41.69091076 | yes |
| 17 | (3,3,3,4,4) | WITNESS | 0.3213635993 | yes |
| 18 | (3,3,4,4,4) | GAP | 20.35624533 | yes |
| 19 | (3,4,4,4,4) | WITNESS | 0.07215458633 | yes |
| 20 | (4,4,4,4,4) | WITNESS | 0.088188939 | yes |
| 21 | (4,4,4,4,5) | WITNESS | 0.018451569 | yes |
| 22 | (4,4,4,5,5) | WITNESS | 0.26318628 | yes |
| 23 | (4,5,4,5,5) | WITNESS | 0.0040481475 | yes |
| 24 | (4,5,5,5,5) | WITNESS | 0.06027056 | yes |
| 25 | (5,5,5,5,5) | WITNESS | 1.7271431e-06 | yes |
| 30 | (6,6,6,6,6) | WITNESS | 4.8606138e-07 | yes |
| 35 | (7,7,7,7,7) | WITNESS | 3.9886479e-14 | yes |
| 40 | (8,8,8,8,8) | WITNESS | 4.7068795e-15 | yes |
| 45 | (9,9,9,9,9) | WITNESS | 1.1800852e-24 | yes |
| 50 | (10^5) | WITNESS | 5.8567945e-26 | yes |

Exact fractions: `data/task1_n15_40.json`, `data/task1_n41_50.json`.
C confirmations: `data/task1_C_confirm.json`, `data/task1_C_confirm_41_50.json`.

Exact ratios at the two gap orders (Py, confirmed by C on both sides):

\[
\frac{L(C_5[3,3,3,3,4])}{L(K_{8,8})}
=\frac{6360639670045255013}{152566579964369860}
\approx 41.69091076,
\]

\[
\frac{L(C_5[3,3,4,4,4])}{L(K_{9,9})}
=\frac{2508406335896865033959107}{123225393272519524676648}
\approx 20.35624533.
\]

Exact Thm 5 ratio at n=15:

\[
\frac{L(C_5[3^5])}{L(K_{7,8})}
=\frac{1535460265677768}{7628328998218493}
\approx 0.2012839596.
\]

---

## Result C — Balanced C5 tower, m=1..16 (n=5m)

C (`c5blowup m`, mpq grid) and Py (layer DP) agree on L for m=1..12 (exact
fraction match). C computed m=13..16; Py independently computed m=13..16
and the exponents match to all printed digits.

\[
e(n)=-\frac{\log_2 L(C_5[m^5])}{n^2}
\]
with \(\log_2\) of the exact fraction (Decimal, 80 digits), not float division
of numerator by denominator.

| m | n | e(n) | C+Py |
|---|---|---|---|
| 1 | 5 | 0.323072623882 | yes |
| 2 | 10 | 0.367950079860 | yes |
| 3 | 15 | 0.396158266253 | yes |
| 4 | 20 | 0.406723439501 | yes |
| 5 | 25 | **0.408836513709** (peak) | yes |
| 6 | 30 | 0.408260351042 | yes |
| 8 | 40 | 0.405495235918 | yes |
| 10 | 50 | 0.402540743364 | yes |
| 12 | 60 | 0.399772502395 | yes |
| 14 | 70 | 0.397274746805 | C; Py e matches |
| 16 | 80 | **0.395088438283** | C; Py e matches |

Exact L: C logs in `data/task2_C_tower_m1_12.txt`, `data/task2_C_tower_m13_16.txt`;
Py JSON in `data/task2_py_tower_m1_12.json`, `data/task2_py_tower_m13_16.json`.

**Wall times (n=80 tower), same machine as `data/hardware.txt`:**
12th Gen Intel Core i5-12600K (8 cores / 16 threads), 16 GB RAM, WSL2 Linux.

| engine | range | wall | max RSS |
|---|---|---|---|
| C (`c5blowup`, mpq grid) | m=1–12 (n=5..60) | **17.50 s** | 65 MB |
| C | m=13–16 (n=65..80) | **170.79 s** | 361 MB |
| C | m=1–16 total | **188.3 s** | 361 MB |
| Py (layer DP) | m=1–12 | **7.76 s** (sum of per-m) | — |
| Py | m=16 (n=80) | **14.44 s** | — |
| Py | m=13–16 | **43.2 s** | — |

The C n=80 run (1.42 million quotient states of `mpq`) is the piece that sits
beyond the paper's \(O(n\,2^n)\) subset-recurrence reach.

---

## Result D — Growth constant

Entropy constant \(1/(4\ln 2)\approx 0.360673760222\).
Bipartite constant \(1/2-1/(8\ln 2)\approx 0.319663119889\).

Least squares of \(e(n)=c+a(\log_2 n)/n+b/n\) on the tower:

| sample | c | a | b | rms | c − entropy |
|---|---|---|---|---|---|
| odd m≥5 | 0.366003 | 0.754 | −2.431 | 5.1e-5 | +0.00533 |
| even m≥6 | 0.365199 | 0.780 | −2.535 | 4.5e-5 | +0.00453 |
| all m≥8 | 0.363442 | 0.850 | −2.845 | 1.0e-5 | +0.00277 |
| all m≥10 | 0.362818 | 0.878 | −2.967 | 2.9e-6 | +0.00214 |

All fitted c are **strictly above** the entropy constant, but the gap shrinks as
the window moves to larger m. The honest reading of these windows is
**c ≈ 0.363–0.366**, above \(1/(4\ln 2)\) in the fits **but not provably**.
Per-point residuals and the odd/even splits are in `data/task3_fit.json`
(`fits.*.points` and `parity_e_rows`). We do **not** claim a proof that
\(\lim e(n)>1/(4\ln 2)\).

The theoretically cleaner leading term is the min-cost insertion path W
(sum of \(\log_2\binom{k}{d}\) along a module sequence). Float DP gives
W/n² = 0.2867, 0.3073, 0.3189, 0.3264, 0.3317 at m=4,6,8,10,12, still rising
and still below e(n) by ~ \(\log_2(n!)/n^2\). If W/n² → c, then e(n)→c
because \(\log_2(n!|\mathrm{Aut}|)/n^2\to 0\). Larger-m min-cost DP would
tighten the constant; we did not push it far enough for a theorem-grade c.

**Open question for the authors.** If the C5-family constant is strictly
larger than the entropy constant, then balanced C5 blow-ups are eventually
rarer than a typical process output, yet they remain far from the global
minimum 1/2. Some other family must then overtake C5 blow-ups at a finite
order. The numerics are compatible with this but do not prove it.

---

## Result E — Gap orders n=16 and n=18

### Aut threshold (Lemma 17 / bound (1)) — exact rationals

Bound (1): \(L(G)\ge 1/(|\mathrm{Aut}(G)|\,P_n)\) with
\[
P_n=\prod_{i=1}^n\binom{i-1}{\lfloor(i-1)/2\rfloor}.
\]
Hence \(L(G)<L(K_{\mathrm{bal}})\) forces
\(|\mathrm{Aut}(G)| > T:=1/(L(K_{\mathrm{bal}})P_n)\). The integer threshold is
the least integer strictly greater than that rational.

For n=16 the central-binomial factors of \(P_{16}\) are
\[
1,1,2,3,6,10,20,35,70,126,252,462,924,1716,3432,6435,
\]
so \(P_{16}=9061429740221589431500800000\). Recurrence (6) gives
\[
L(K_{8,8})=\frac{7628328998218493}{107497751435913928580204708114517196800000000000},
\]
and therefore
\[
T_{16}=\frac{11863221866496000000}{7628328998218493}\approx 1555.1534116143.
\]
Integer threshold **1556**, matching the paper. For n=18,
\[
T_{18}=\frac{1016440849521377280000000}{2200453451294991512083}\approx 461.9233589891,
\]
integer threshold **462**. Full factor tables for n=13,14,16,18:
`data/task4_T_exact.json` (all `matches_paper: true`).
C subset recurrence agrees with (6) on \(L(K_{8,8})\) and \(L(K_{9,9})\).

### Paper §5 coverage (lifted, not reinvented)

The paper does **not** enumerate graphs. It enumerates conjugacy classes of
subgroups \(\Gamma\le S_n\) with \(|\Gamma|\ge T\), then takes unions of
\(\Gamma\)-orbitals on unordered pairs. Coverage: every subgroup of \(S_n\)
of order \(\ge T\) occurs in a chain of maximal subgroups descending from
\(S_n\), so appears (up to conjugacy in \(S_n\)) in that walk. Over-generation
is harmless; under-coverage would silently break the proof.

Forced-twin test (paper §5): if some \(u\neq v\) have \(\{u,w\}\) and \(\{v,w\}\)
in the same pair-orbit for every \(w\notin\{u,v\}\), then every \(\Gamma\)-invariant
graph makes \(u,v\) twins. At n=16 the paper reports 116597 classes above
T=1556, of which 116553 force twins (not expanded) and 44 twin-free classes
give 128 isomorphism classes, none below \(K_{8,8}\). Twin graphs are covered
only as blow-ups through base order 8, so unrestricted n=16 remains open there.

A pruning “product of constituent orders \(< T\)” is valid. A claim “covered
because the transitive constituents were enumerated” is **not** automatically
valid: \(\mathrm{Aut}(G)\) is a subdirect product of its constituents, so
\(|\mathrm{Aut}|\) can lie far below the product. The coverage lemma we use
is exactly the paper’s: for every 2-closed \(Q\le S_n\) with \(|Q|\ge T\),
some enumerated \(P\le Q\), hence every \(Q\)-invariant graph appears among
the unions of \(P\)-orbitals.

### Vertex-transitive graphs (unconditional)

If \(\mathrm{Aut}(G)\) is transitive then \(G\) is vertex-transitive.
We evaluate the **complete** Holt–Royle census of vertex-transitive graphs
(Zenodo [10.5281/zenodo.4010122](https://doi.org/10.5281/zenodo.4010122),
arXiv:1811.09015), which is the standard distribution of the McKay–Royle
VT census: every VT graph of order n, **including disconnected** graphs
(equal-order disjoint unions of a connected VT graph) and both a graph and
its complement (stored by valency). \(L(G)=L(\overline{G})\) halves the work.

| n | census lines (all valencies) | unique labelled | complement-pairs evaluated | \|Aut\|≥T | min \(L/L(K_{\mathrm{bal}})\) | beat \(K_{\mathrm{bal}}\)? |
|---|---|---|---|---|---|---|
| 16 | 286 | 286 | 143 | 19 | **1** (attained by \(K_{8,8}\cong\overline{2K_8}\)) | no |
| 18 | 380 | 380 | 190 | 31 | **1** (attained by \(K_{9,9}\cong\overline{2K_9}\)) | no |

At n=16 the 286 lines split as 272 connected + 14 disconnected, matching the
classical connected-VT count. Tables: `data/task4_vt16.tsv`, `data/task4_vt18.tsv`
(id, graph6, \|Aut\|, exact L, exact ratio). JSON with provenance:
`data/task4_vt16.json`, `data/task4_vt18.json`. Evaluator: Py subset recurrence
with nauty `labelg` canon and `code/g6aut` for \|Aut\| (16! printed in
scientific notation; parsed as an integer via `Decimal`). Wall: 60.4 s (n=16),
373.8 s (n=18).

**Theorem (VT case).** Among all vertex-transitive graphs on 16 vertices,
including disconnected ones, the unique minimum of \(L\) (up to complement)
is \(K_{8,8}\). The same holds at n=18 with unique minimum \(K_{9,9}\).
This does **not** use bound (1) except as a diagnostic: every census graph
was evaluated, including those with \|Aut\| < T.

Selected n=16 VT ratios (full table in the TSV):

| graph6 (sparser of G, complement) | \|Aut\| | edges | \(L/L(K_{8,8})\) |
|---|---|---|---|
| `O~~~~{??G@_F?N?N_Fw@~` (\(2K_8\)) | 3251404800 | 56 | 1 |
| `O{eCN{??JpnFXkXkbu?|` | 12288 | 56 | ≈ 16.565 |
| `O{eCN|_WGpbFXbXbbu?|` | 4096 | 56 | ≈ 19.226 |
| `OsaCB@_EWrKrXeFwB{B{?` | 4096 | 48 | ≈ 193.49 |

### Intransitive case (Step 3) — stalled

The remaining graphs with \|Aut\| ≥ 1556 are intransitive. Replicating the
paper’s S_{16} maximal-subgroup descent (116597 classes) needs GAP
(`MaximalSubgroupClassReps`, conjugacy in \(S_{16}\)). This machine has no
GAP install and no passwordless sudo, so the walk was not run. Scripts ready
for when GAP is available: `scripts/s16_walk.g` (logs every class, k, forced-twin
flag; expands only twin-free orbitals) and `scripts/vt16.g` (transitive-group
orbital cross-check). Without that walk we **do not** claim a theorem for
unrestricted n=16. n=18 is strictly harder (T=462 admits far more groups).

### Blow-ups of C7, C9, Petersen (Task 4a)

All **positive** part vectors, one representative per Aut(H)-orbit.
None beat \(K_{\mathrm{bal}}\). Ratios are enormous.

| base | n | orbits | best parts | L / L(K_bal) |
|---|---|---|---|---|
| C7 | 16 | 375 | (1,1,3,3,2,3,3) | ≈ 1.55e5 |
| C9 | 16 | 375 | (1,1,1,1,2,3,2,3,2) | ≈ 6.00e8 |
| Petersen | 16 | 76 | (1,1,1,1,2,1,2,2,1,4) | ≈ 3.16e4 |
| C7 | 18 | 912 | (2,2,3,3,2,3,3) | ≈ 3.77e5 |
| C9 | 18 | 1387 | (2^9) | ≈ 1.48e9 |
| Petersen | 18 | 291 | (1,1,2^8) | ≈ 1.37e5 |

Exact fractions: `data/task4a_blowups.json`. Py only (same Prop. 11 engine
that matched C on C5 checkpoints).

### Complete multipartite graphs (Task 4b, partial)

Every complete multipartite graph on n vertices is evaluated by the subset
recurrence with Aut order \(\prod_j (n_j!)\prod_s (m_s!)\). L is complement-invariant,
so this is also every disjoint union of cliques.

At **n=16** there are 230 partitions into ≥2 parts. The unique minimum in this
class is **K_{8,8}** (ratio 1). The runner-up is \(K_{7,9}\) at ratio ≈ 914.
**No complete-multipartite counterexample at n=16.**

At **n=18** there are 384 partitions. The unique minimum is **K_{9,9}**.
**No complete-multipartite counterexample at n=18.**

Exact records: `data/task4b_multipartite.json`.

### What is *not* proved

A proof of Conjecture 1 at n=16 still needs the intransitive half of the §5
walk (or an equivalent 2-closed-group census). The **vertex-transitive**
subclass is settled (unique min \(K_{8,8}\), including disconnected). Twin-free
graphs above T=1556 were already settled by the paper (128 iso classes).
Twin-forcing classes with quotient order ≥ 9 remain the obstruction, as in
the paper. n=18 is open on the same intransitive side. C has a `graph6`
mode for n≤20.

---

## Data and code

- `data/task1_n15_40.json`, `data/task1_n41_50.json` — family table
- `data/task2_*.txt`, `data/task2_*.json` — tower
- `data/task3_fit.json` — e(n), fits, residuals, odd/even splits
- `data/task4_T_exact.json` — exact T, P_n, central-binomial factors
- `data/task4a_blowups.json`, `data/task4b_multipartite.json`
- `data/task4_vt16.{json,tsv}`, `data/task4_vt18.{json,tsv}` — VT census
- `code/independent.py`, `likelihood.c`, `code/g6aut.c`, `scripts/`

Wall-clock: Task 0 <5s; table n=15–40 ~70s Py + 3s C; n=41–50 ~4 min Py + 13s C;
tower m=1–16 **188.3 s C** / ~51 s Py (see Result C); C7/C9/Petersen ~4 min;
n=16 complete multipartite ~85s; VT n=16 60 s; VT n=18 374 s.
