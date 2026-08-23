# Maximal-subgroup descent in S_n of conjugacy-class representatives
# with |H| >= T. Coverage (paper §5): every H <= S_n with |H| >= T lies
# in a chain S_n = H_0 > H_1 > ... > H_k = H of maximal subgroups, so
# appears (up to conjugacy in S_n) in this walk. Over-generation is
# harmless; conjugacy is taken in S_n, not in the parent.
#
# For each kept representative P we record:
#   index, Size, k_pair_orbitals, forced_twin (bool), n_graphs = 2^k
# Graphs with k <= MAXK are written as graph6; larger k are logged only.

n := 16;
T := 1556;
MAXK := 20;
outdir := "/home/james/src/graph-likelihood/data/";
logpath := Concatenation(outdir, "task4_s16_groups.log");
g6path := Concatenation(outdir, "task4_s16_raw.g6");
sumpath := Concatenation(outdir, "task4_s16_summary.txt");

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
            if k + b <= Length(bits) and bits[k + b] then
                c := c + 1;
            fi;
        od;
        Add(s, CHAR_INT(63 + c));
        k := k + 6;
    od;
    return s;
end;

PairOrbitalsFast := function(G, n)
    local used, orbs, i, j, stack, a, b, p, gens, o, u, tmp;
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

ForcedTwin := function(orbs, n)
    # true if some u < v have {u,w} and {v,w} in the same pair-orbit
    # for every w not in {u,v} (paper §5 forced-twin test).
    local label, i, j, t, e, u, v, w, ou, ov, a, b;
    label := List([1..n], x -> List([1..n], y -> 0));
    t := 0;
    for i in [1..Length(orbs)] do
        for e in orbs[i] do
            label[e[1]][e[2]] := i;
            label[e[2]][e[1]] := i;
        od;
    od;
    for u in [1..n] do
        for v in [u+1..n] do
            ou := true;
            for w in [1..n] do
                if w = u or w = v then
                    continue;
                fi;
                if label[u][w] <> label[v][w] then
                    ou := false;
                    break;
                fi;
            od;
            if ou then
                return true;
            fi;
        od;
    od;
    return false;
end;

IsConjugateSn := function(S, H, seen)
    local K;
    for K in seen do
        if Size(K) = Size(H) and RepresentativeAction(S, H, K) <> fail then
            return true;
        fi;
    od;
    return false;
end;

S := SymmetricGroup(n);
seen := [];
todo := [S];
log := OutputTextFile(logpath, false);
g6out := OutputTextFile(g6path, false);
PrintTo(log, "# S_", n, " maximal-subgroup descent, T=", T, "\n");

nraw := 0;
nskipk := 0;
nforced := 0;
ntwinfree := 0;
processed := 0;

while Length(todo) > 0 do
    G := todo[Length(todo)]; Remove(todo);
    if Size(G) < T then
        continue;
    fi;
    if IsConjugateSn(S, G, seen) then
        continue;
    fi;
    Add(seen, G);
    processed := processed + 1;
    orbs := PairOrbitalsFast(G, n);
    k := Length(orbs);
    ft := ForcedTwin(orbs, n);
    if ft then
        nforced := nforced + 1;
    else
        ntwinfree := ntwinfree + 1;
    fi;
    PrintTo(log, processed, " |G|=", Size(G), " k=", k, " forced_twin=", ft, "\n");
    if k <= MAXK and not ft then
        # Twin-free classes are few (paper: 44). Expand their orbital graphs.
        # Twin-forcing classes are logged; expanding all 116k of them is the
        # remaining n=16 search (paper left it open). Over-generation of
        # twin-free graphs is harmless.
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
            WriteLine(g6out, Graph6OfAdj(n, adj));
            nraw := nraw + 1;
        od;
    else
        nskipk := nskipk + 1;
        PrintTo(log, "# SKIP k=", k, " |G|=", Size(G), "\n");
    fi;
    if RemInt(processed, 5) = 0 then
        Print("classes=", processed, " todo=", Length(todo), " raw=", nraw,
              " last|G|=", Size(G), " k=", k, "\n");
    fi;
    # children: maximal subgroups of G of order >= T
    maxs := MaximalSubgroupClassReps(G);
    for M in maxs do
        if Size(M) >= T then
            Add(todo, M);
        fi;
    od;
od;

PrintTo(log, "# classes=", Length(seen), " twin_force=", nforced,
        " twin_free=", ntwinfree, " skip_k=", nskipk, " raw=", nraw, "\n");
PrintTo(sumpath, "classes=", Length(seen), "\n",
        "twin_force=", nforced, "\n",
        "twin_free=", ntwinfree, "\n",
        "skip_k=", nskipk, "\n",
        "raw_graphs=", nraw, "\n");
CloseStream(log);
CloseStream(g6out);
Print("DONE classes=", Length(seen), " raw=", nraw, " forced=", nforced,
      " twinfree=", ntwinfree, " skipk=", nskipk, "\n");
QUIT;
