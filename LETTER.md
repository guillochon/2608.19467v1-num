# Note to Severini–Weisstein on exact C5-family witnesses and the n=16,18 gaps

We independently reimplemented the deletion recurrence (Prop. 8), the
blow-up recurrence (Prop. 11), and the bipartite recurrence (6) in C11/GMP
and, separately, in integer-scaled Python. Both codes reproduce Theorem 5
and the additional exact checkpoints in your message.

**C5 family, n=15–50.** Enumerating all positive 5-tuples up to D5, the
family minimum is always a most-balanced tuple (for n=5k+3 the two smaller
parts sit non-adjacent on the cycle). It beats \(K_{\lfloor n/2\rfloor,\lceil n/2\rceil}\)
at n=15, 17, and every 19≤n≤50, and loses at n=16 and n=18. Exact ratios:

| n | parts | \(L(C_5[m])/L(K_{\mathrm{bal}})\) |
|---|---|---|
| 15 | (3,3,3,3,3) | 1535460265677768 / 7628328998218493 ≈ 0.20128396 |
| 16 | (3,3,3,3,4) | 6360639670045255013 / 152566579964369860 ≈ 41.69091 |
| 17 | (3,3,3,4,4) | 9900038976142070369457 / 30806348318129881169162 ≈ 0.3213636 |
| 18 | (3,3,4,4,4) | 2508406335896865033959107 / 123225393272519524676648 ≈ 20.35625 |

These are upper bounds on \(\min L\), not global minima. The n=16 and n=18
gaps in the C5 family are therefore real, so the failure set of Conjecture 1
is not monotone at this scale.

**Tower.** Balanced blow-ups \(C_5[m^5]\) through n=80 give
\(e(n)=-\log_2 L/n^2\) peaking at n=25 (\(e=0.40883651\)) and falling to
\(e(80)=0.39508844\). Fits \(e=c+a(\log_2 n)/n+b/n\) on the largest terms
produce c ≈ 0.363–0.366, all slightly above the entropy constant
\(1/(4\ln 2)\approx 0.36067\), with the gap shrinking as n grows. We do not
claim a proof that the C5 constant exceeds entropy; if it does, some other
family must overtake C5 blow-ups at finite order.

**n=16 and n=18 outside C5.** The Lemma 17 threshold is |Aut|≥1556 (n=16)
and |Aut|≥462 (n=18); we recomputed T as the exact rationals
\(11863221866496000000/7628328998218493\) and
\(1016440849521377280000000/2200453451294991512083\). Exhaustive Aut(H)-orbits
of blow-ups of C7, C9, and Petersen all lose to \(K_{\mathrm{bal}}\) by
factors \(>10^4\). Every complete multipartite graph on 16 (resp. 18) vertices
has unique minimum \(K_{8,8}\) (resp. \(K_{9,9}\)).

The complete Holt–Royle / McKay–Royle census of vertex-transitive graphs,
including disconnected ones (286 graphs on 16 vertices, 380 on 18), was
evaluated by the subset recurrence, using \(L(G)=L(\overline{G})\) to halve
the work. In the VT class the unique minimum is \(K_{8,8}\) at n=16 and
\(K_{9,9}\) at n=18. This settles the transitive case unconditionally.

We did not rerun geng enumeration (infeasible at these orders). We also did
not complete the paper’s §5 maximal-subgroup descent in \(S_{16}\) (116597
classes above T=1556): this environment has no GAP, which that walk needs.
Twin-free classes above the threshold were already checked in the paper
(128 iso types, none below \(K_{8,8}\)). So we do **not** prove Conjecture 1
at unrestricted n=16 or n=18. Has the §5 enumeration been extended to 16/18,
in particular to the 116553 twin-forcing classes (or equivalently blow-ups
of twin-free bases of order 9–15 with \(|\mathrm{Aut}|\ge 1556\))?

Code, exact JSON, and a longer write-up: `RESULTS.md` and `data/` in this
repository (C11/GMP binary `likelihood`, Python `code/independent.py`).
