# Twin-free orbital graphs of primitive permutation groups of degree 16
# with |G| >= T. This is a proper subset of the S16 maximal-subgroup walk
# (every primitive G appears in that walk); it is only an early harvest.

n := 16;
T := 1556;
outdir := "/home/james/src/graph-likelihood/data/";
g6path := Concatenation(outdir, "task4_s16_primitive_twinfree.g6");
logpath := Concatenation(outdir, "task4_s16_primitive.log");

Graph6OfAdj := function(nn, adj)
    local s, bits, b, i, j, k, c;
    s := [CHAR_INT(63 + nn)];
    bits := [];
    for j in [1..nn-1] do
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

PairOrbitalsFast := function(G, nn)
    local used, orbs, i, j, stack, a, b, p, gens, o, u, tmp;
    gens := GeneratorsOfGroup(G);
    used := List([1..nn], x -> BlistList([1..nn], []));
    orbs := [];
    for i in [1..nn] do
        for j in [i+1..nn] do
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

ForcedTwin := function(orbs, nn)
    local label, i, e, u, v, w, ok;
    label := List([1..nn], x -> List([1..nn], y -> 0));
    for i in [1..Length(orbs)] do
        for e in orbs[i] do
            label[e[1]][e[2]] := i;
            label[e[2]][e[1]] := i;
        od;
    od;
    for u in [1..nn] do
        for v in [u+1..nn] do
            ok := true;
            for w in [1..nn] do
                if w <> u and w <> v and label[u][w] <> label[v][w] then
                    ok := false;
                    break;
                fi;
            od;
            if ok then
                return true;
            fi;
        od;
    od;
    return false;
end;

ExpandOrbitalGraphs := function(orbs, nn, path)
    local k, m, mask, adj, t, e, nout;
    k := Length(orbs);
    m := 2^k;
    nout := 0;
    for mask in [0..m-1] do
        adj := List([1..nn], x -> List([1..nn], y -> false));
        for t in [0..k-1] do
            if RemInt(QuoInt(mask, 2^t), 2) = 1 then
                for e in orbs[t+1] do
                    adj[e[1]][e[2]] := true;
                    adj[e[2]][e[1]] := true;
                od;
            fi;
        od;
        AppendTo(path, Graph6OfAdj(nn, adj), "\n");
        nout := nout + 1;
    od;
    return nout;
end;

PrintTo(g6path, "");
PrintTo(logpath, "# primitive groups degree ", n, " |G|>=", T, "\n");
nprim := NrPrimitiveGroups(n);
nabove := 0;
ntf := 0;
nraw := 0;
t0 := Runtime();
for i in [1..nprim] do
    G := PrimitiveGroup(n, i);
    if Size(G) < T then
        continue;
    fi;
    nabove := nabove + 1;
    orbs := PairOrbitalsFast(G, n);
    k := Length(orbs);
    ft := ForcedTwin(orbs, n);
    AppendTo(logpath, i, " |G|=", Size(G), " k=", k, " forced_twin=", ft, "\n");
    Print("prim ", i, "/", nprim, " |G|=", Size(G), " k=", k, " ft=", ft, "\n");
    if (not ft) and k <= 16 then
        ntf := ntf + 1;
        nraw := nraw + ExpandOrbitalGraphs(orbs, n, g6path);
    fi;
od;
AppendTo(logpath, "# nprim=", nprim, " above_T=", nabove, " twinfree=", ntf,
         " raw=", nraw, " t_ms=", Runtime()-t0, "\n");
Print("DONE nprim=", nprim, " above_T=", nabove, " twinfree=", ntf,
      " raw=", nraw, " t_s=", QuoInt(Runtime()-t0, 1000), "\n");
QUIT;
