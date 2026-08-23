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
the window moves to larger m. An honest interval from this range of n is
**c ∈ (0.360, 0.367)**; we do **not** claim a proof that
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

### Aut threshold (Lemma 17 / bound (1))

\(P_n=\prod_{i=1}^n\binom{i-1}{\lfloor(i-1)/2\rfloor}\).
A graph with \(L(G)<L(K_{\mathrm{bal}})\) must have
\(|\mathrm{Aut}(G)| > 1/(L(K_{\mathrm{bal}})P_n)\).

| n | T = 1/(L P_n) | necessary \|Aut\| | matches paper |
|---|---|---|---|
| 16 | 11863221866496000000/7628328998218493 ≈ 1555.1534 | **≥ 1556** | yes |
| 18 | 1016440849521377280000000/2200453451294991512083 ≈ 461.9234 | **≥ 462** | yes |

(Py recurrence (6) for L; C subset agrees on L(K_8,8) and L(K_9,9).)

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

A proof that Conjecture 1 holds at n=16 (resp. 18) requires enumerating **all**
graphs with |Aut| ≥ 1556 (resp. 462), or an equivalent group-orbit census.
The paper already checked the twin-free n=16 graphs above the threshold (128
isomorphism classes) and left twin graphs with large quotients open. We have
not completed that census. C has a `graph6` mode for n≤20 to continue the
search.

---

## Data and code

- `data/task1_n15_40.json`, `data/task1_n41_50.json` — family table
- `data/task2_*.txt`, `data/task2_*.json` — tower
- `data/task3_fit.json` — e(n) and fits
- `data/task4a_blowups.json`, `data/task4b_multipartite.json`
- `code/independent.py`, `likelihood.c`, `scripts/`

Wall-clock (order of magnitude): Task 0 <5s; table n=15–40 ~70s Py + 3s C;
n=41–50 ~4 min Py + 13s C; tower m=1–16 ~3 min C, ~1 min Py; C7/C9/Petersen
~4 min; n=16 complete multipartite ~85s.
