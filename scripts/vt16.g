# Generate all vertex-transitive graphs of degree n with |Aut| >= T
# as unions of pair-orbitals of TransitiveGroup(n,i) with Size >= T.
#
# Coverage: every vertex-transitive G has Aut(G) conjugate in S_n to some
# TransitiveGroup(n,i). If |Aut(G)| >= T then that representative has
# Size >= T, and G (after relabelling) is a union of its pair-orbitals.
# Over-generation is expected (a graph appears for every transitive
# subgroup of its Aut); we canonicalise later.
#
# Output: one graph6 string per line to stdout; log to stderr.

LoadPackage("transgrp");;
n := 16;
T := 1556;
logfile := "/home/james/src/graph-likelihood/data/task4_vt16_groups.log";
outpath := "/home/james/src/graph-likelihood/data/task4_vt16_raw.g6";

Graph6OfAdj := function(n, adj)
    local s, bits, b, i, j, k, c;
    s := [CHAR_INT(63 + n)];
    bits := [];
    for j in [1..n-1] do
        for i in [0..j-1] do
            Add(bits, adj[i+1][j+1]);
        od;
    od;
    k := 0;
    while k < Length(bits) do
        c := 0;
        for b in [1..6] do
            c := 2 * c;
            if k + b <= Length(bits) and bits[k+b] then
                c := c + 1;
            fi;
        od;
        Add(s, CHAR_INT(63 + c));
        k := k + 6;
    od;
    return s;
end;

PairOrbitals := function(G, n)
    local pairs, orbs, used, i, j, a, b, o, p, img, rep, adj0;
    pairs := [];
    used := List([1..n], x -> List([1..n], y -> false));
    orbs := [];
    for i in [1..n] do
        for j in [i+1..n] do
            if used[i][j] then
                continue;
            fi;
            o := [];
            for p in G do
                a := i^p; b := j^p;
                if a > b then
                    img := [b, a];
                else
                    img := [a, b];
                fi;
                if not used[img[1]][img[2]] then
                    used[img[1]][img[2]] := true;
                    Add(o, img);
                fi;
            od;
            Add(orbs, o);
        od;
    od;
    return orbs;
end;

# Faster pair-orbitals via orbit algorithm on a representative pair
PairOrbitalsFast := function(G, n)
    local used, orbs, i, j, stack, a, b, na, nb, p, gens, o, u, v, tmp;
    gens := GeneratorsOfGroup(G);
    used := List([1..n], x -> BlistList([1..n], []));
    orbs := [];
    for i in [1..n] do
        for j in [i+1..n] do
            if used[i][j] then
                continue;
            fi;
            o := [];
            stack := [[i, j]];
            used[i][j] := true;
            Add(o, [i, j]);
            while Length(stack) > 0 do
                u := stack[Length(stack)]; Remove(stack);
                for p in gens do
                    a := u[1]^p; b := u[2]^p;
                    if a > b then
                        tmp := a; a := b; b := tmp;
                    fi;
                    if not used[a][b] then
                        used[a][b] := true;
                        Add(o, [a, b]);
                        Add(stack, [a, b]);
                    fi;
                od;
            od;
            Add(orbs, o);
        od;
    od;
    return orbs;
end;

out := OutputTextFile(outpath, false);
log := OutputTextFile(logfile, false);
ng := NrTransitiveGroups(n);
PrintTo(log, "# n=", n, " T=", T, " NrTransitiveGroups=", ng, "\n");
kept := 0;
skipped_small := 0;
skipped_k := 0;
nraw := 0;
MAXK := 22;  # 2^22 = 4M graphs; larger k logged and skipped (none expected for |G|>=T)

for idx in [1..ng] do
    G := TransitiveGroup(n, idx);
    sz := Size(G);
    if sz < T then
        skipped_small := skipped_small + 1;
        continue;
    fi;
    orbs := PairOrbitalsFast(G, n);
    k := Length(orbs);
    kept := kept + 1;
    PrintTo(log, idx, " |G|=", sz, " k=", k, "\n");
    if k > MAXK then
        skipped_k := skipped_k + 1;
        PrintTo(log, "# SKIP k=", k, " group ", idx, "\n");
        continue;
    fi;
    m := 2^k;
    for mask in [0..m-1] do
        adj := List([1..n], x -> List([1..n], y -> false));
        for t in [0..k-1] do
            if RemInt(QuoInt(mask, 2^t), 2) = 1 then
                for e in orbs[t+1] do
                    adj[e[1]][e[2]] := true;
                    adj[e[2]][e[1]] := true;
                od;
            fi;
        od;
        g6 := Graph6OfAdj(n, adj);
        WriteLine(out, g6);
        nraw := nraw + 1;
    od;
    if RemInt(kept, 10) = 0 then
        Print("kept ", kept, " groups, raw graphs ", nraw, " last k=", k, " idx=", idx, "\n");
    fi;
od;

PrintTo(log, "# kept_groups=", kept, " skipped_small=", skipped_small,
        " skipped_k=", skipped_k, " raw_graphs=", nraw, "\n");
CloseStream(out);
CloseStream(log);
Print("DONE kept=", kept, " raw=", nraw, " skipped_small=", skipped_small, " skipped_k=", skipped_k, "\n");
QUIT;
