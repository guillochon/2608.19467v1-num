# Exact graph-likelihood witnesses

This is a source-only C11/GMP evaluator for selected graph families in
Severini--Weisstein, *The minimum of the graph likelihood*
(arXiv:2608.19467v1).  Do **not** build it on the 1 GiB review host.

For a graph `G`, it evaluates the deletion recurrence

```text
F(empty) = 1
F(S) = sum(v in S) F(S-v) / binom(|S|-1, deg_{G[S]}(v))
L(G) = F(V) / (n! |Aut(G)|).
```

## Supported families

```sh
./likelihood bipartite 7 8
./likelihood c5blowup 3
./likelihood c5parts 2 3 4 3 2
```

`bipartite a b` uses the ordinary integer-scaled subset recurrence and is
intended for modest orders.  `c5parts a b c d e` is the blow-up of the five
cycle by five positive independent modules in cyclic order.

Unlike a subset implementation, the C5 mode uses the fact that the recurrence
is constant on induced subgraphs having the same five remaining module sizes.
It evaluates the exact recurrence on `(a+1)(b+1)(c+1)(d+1)(e+1)` quotient
states.  Thus balanced C5 blow-ups can be explored well beyond the 2^n range,
subject to GMP time and memory.  This is a structured-family evaluator, not an
exhaustive search over all graphs and not a proof of a global minimum.

The automorphism order used for C5 parts is

```text
(product_i m_i!) * |{dihedral symmetries of C5 preserving (m_0,...,m_4)}|.
```

This formula is valid because all parts are positive independent modules and
C5 has no twins.

## Build and regression test

A capable machine needs a C11 compiler and GMP development headers/library.

```sh
make
make check
```

`make check` verifies the two exact n=15 values stated in Theorem 5.  The
program emits progress information to stderr and the final rational result to
stdout.
