# Maximal-subgroup descent in S_n for conjugacy-class representatives
# with |H| >= T. Coverage (paper §5): every H <= S_n with |H| >= T lies
# in a chain S_n = H_0 > H_1 > ... > H_k = H of maximal subgroups, so
# appears (up to conjugacy in S_n) in this walk. Over-generation is
# harmless; conjugacy is taken in S_n, not in the parent.
#
# Pass 1: log every class (size, k, forced_twin, generators). Expand
# twin-free orbital graphs with k <= MAXK (paper: 44 classes). Twin-forcing
# classes are logged so they can be expanded later; skipping them is how
# the paper left unrestricted n=16 open.

n := 16;
T := 1556;
MAXK := 16;
outdir := "/home/james/src/graph-likelihood/data/";
logpath := Concatenation(outdir, "task4_s16_groups.log");
g6path := Concatenation(outdir, "task4_s16_twinfree.g6");
genspath := Concatenation(outdir, "task4_s16_gens.g");
sumpath := Concatenation(outdir, "task4_s16_summary.txt");
ckptpath := Concatenation(outdir, "task4_s16_checkpoint.txt");

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

# Cheap conjugacy invariant (no pair-orbitals). Element-order histogram
# for small groups almost always separates classes, so IsConjugate is rare.
Fingerprint := function(G, nn)
    local ol, ordh;
    ol := Collected(List(Orbits(G, [1..nn]), Length));
    if Size(G) <= 20000 then
        ordh := Collected(List(G, Order));
    else
        ordh := Collected(List(GeneratorsOfGroup(G), Order));
    fi;
    return Concatenation(String(Size(G)), ":", String(ol), ":", String(ordh));
end;

# Store generating sets, not live groups: 50k live subgroups of S16
# otherwise grow to many GB (observed 6 GB at 13k classes).
PackG := function(G)
    return rec(gens := List(GeneratorsOfGroup(G)), sz := Size(G));
end;

UnpackG := function(r)
    local H;
    H := Group(r.gens);
    SetSize(H, r.sz);
    return H;
end;

HeapPush := function(heap, r)
    local i, p, t;
    Add(heap, r);
    i := Length(heap);
    while i > 1 do
        p := QuoInt(i, 2);
        if heap[p].sz >= heap[i].sz then
            break;
        fi;
        t := heap[p]; heap[p] := heap[i]; heap[i] := t;
        i := p;
    od;
end;

HeapPop := function(heap)
    local r, i, n, l, b, t;
    r := heap[1];
    heap[1] := heap[Length(heap)];
    Remove(heap);
    n := Length(heap);
    i := 1;
    while true do
        l := 2 * i;
        b := i;
        if l <= n and heap[l].sz > heap[b].sz then
            b := l;
        fi;
        if l + 1 <= n and heap[l + 1].sz > heap[b].sz then
            b := l + 1;
        fi;
        if b = i then
            break;
        fi;
        t := heap[i]; heap[i] := heap[b]; heap[b] := t;
        i := b;
    od;
    return r;
end;

HashStr := function(s)
    local h, c;
    h := 0;
    for c in s do
        h := (h * 131 + INT_CHAR(c)) mod 20011;
    od;
    return h + 1;
end;

NHASH := 20011;
bins := List([1..NHASH], x -> []);  # each: list of [key, [groups]]

AlreadySeen := function(S, G, key)
    local h, bucket, entry, K;
    h := HashStr(key);
    bucket := bins[h];
    for entry in bucket do
        if entry[1] = key then
            for K in entry[2] do
                if IsConjugate(S, G, UnpackG(K)) then
                    return true;
                fi;
            od;
            Add(entry[2], PackG(G));
            return false;
        fi;
    od;
    Add(bins[h], [key, [PackG(G)]]);
    return false;
end;

ExpandOrbitalGraphs := function(orbs, nn, g6path)
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
        AppendTo(g6path, Graph6OfAdj(nn, adj), "\n");
        nout := nout + 1;
    od;
    return nout;
end;

PrintTo(logpath, "# S_", n, " maximal-subgroup descent, T=", T, "\n");
PrintTo(g6path, "");
PrintTo(genspath, "# generators of each class representative; read with Read()\n");
PrintTo(ckptpath, "");

S := SymmetricGroup(n);
todo := [];
HeapPush(todo, PackG(S));
nraw := 0;
nskipk := 0;
nforced := 0;
ntwinfree := 0;
nfailedmax := 0;
ndup := 0;
processed := 0;
t0 := Runtime();

# Mark S_n itself as seen so we do not re-enqueue it.
AlreadySeen(S, S, Fingerprint(S, n));

while Length(todo) > 0 do
    G := UnpackG(HeapPop(todo));
    if Size(G) < T then
        continue;
    fi;
    orbs := PairOrbitalsFast(G, n);
    k := Length(orbs);
    processed := processed + 1;
    ft := ForcedTwin(orbs, n);
    if ft then
        nforced := nforced + 1;
    else
        ntwinfree := ntwinfree + 1;
    fi;
    gens := GeneratorsOfGroup(G);
    AppendTo(logpath, processed, " |G|=", Size(G), " k=", k,
             " forced_twin=", ft, " norb=", Length(Orbits(G, [1..n])),
             " ngens=", Length(gens), "\n");
    if not ft then
        AppendTo(genspath, "CLASS", processed, " := rec(idx:=", processed,
                 ", size:=", Size(G), ", k:=", k, ", forced_twin:=", ft,
                 ", gens:=", gens, ");\n");
    fi;
    if (not ft) and k <= MAXK then
        nraw := nraw + ExpandOrbitalGraphs(orbs, n, g6path);
    else
        nskipk := nskipk + 1;
    fi;
    if RemInt(processed, 1) = 0 and processed <= 20 then
        Print("classes=", processed, " todo=", Length(todo), " raw=", nraw,
              " twinfree=", ntwinfree, " forced=", nforced, " dup=", ndup,
              " last|G|=", Size(G), " k=", k, " ft=", ft,
              " t_s=", QuoInt(Runtime()-t0, 1000), "\n");
    elif RemInt(processed, 25) = 0 then
        Print("classes=", processed, " todo=", Length(todo), " raw=", nraw,
              " twinfree=", ntwinfree, " forced=", nforced, " dup=", ndup,
              " last|G|=", Size(G), " k=", k, " ft=", ft,
              " t_s=", QuoInt(Runtime()-t0, 1000), "\n");
        PrintTo(ckptpath, "classes=", processed, " todo=", Length(todo),
                " twinfree=", ntwinfree, " forced=", nforced,
                " failedmax=", nfailedmax, " dup=", ndup, " raw=", nraw,
                " t_ms=", Runtime()-t0, "\n");
    fi;
    caught := CALL_WITH_CATCH(MaximalSubgroupClassReps, [G]);
    if caught[1] = false then
        nfailedmax := nfailedmax + 1;
        AppendTo(logpath, "# FAILED MaximalSubgroupClassReps |G|=", Size(G), "\n");
        Print("# FAILED maxsub |G|=", Size(G), " ", caught[2], "\n");
    else
        for M in caught[2] do
            if Size(M) < T then
                continue;
            fi;
            if AlreadySeen(S, M, Fingerprint(M, n)) then
                ndup := ndup + 1;
            else
                HeapPush(todo, PackG(M));
            fi;
        od;
    fi;
od;

PrintTo(sumpath,
        "classes=", processed, "\n",
        "twin_force=", nforced, "\n",
        "twin_free=", ntwinfree, "\n",
        "skip_or_forced=", nskipk, "\n",
        "raw_twinfree_graphs=", nraw, "\n",
        "failed_maxsub=", nfailedmax, "\n",
        "dup_enqueue=", ndup, "\n",
        "time_ms=", Runtime()-t0, "\n");
AppendTo(logpath, "# classes=", processed, " twin_force=", nforced,
         " twin_free=", ntwinfree, " skip=", nskipk, " raw=", nraw,
         " failed_maxsub=", nfailedmax, " dup=", ndup,
         " time_ms=", Runtime()-t0, "\n");
Print("DONE classes=", processed, " raw=", nraw, " forced=", nforced,
      " twinfree=", ntwinfree, " failedmax=", nfailedmax, " dup=", ndup,
      " t_s=", QuoInt(Runtime()-t0, 1000), "\n");
QUIT;
