# Maximal-subgroup descent in S_n for conjugacy-class representatives
# with |H| >= T. Coverage (paper §5): every H <= S_n with |H| >= T lies
# in a chain S_n = H_0 > H_1 > ... > H_k = H of maximal subgroups, so
# appears (up to conjugacy in S_n) in this walk. Over-generation
# is harmless; conjugacy is taken in S_n, not in the parent.
#
# Pass 1: log every class (size, k, forced_twin, generators). Expand
# twin-free orbital graphs with k <= MAXK (paper: 44 classes).
#
# Seen-set key (no IsConjugate): transgrp id of each transitive
# constituent, |G|, pair-orbital size signature, and 1-dim row
# signatures of the pair-orbit colouring. For transitive groups of
# degree 16 this is a complete conjugacy invariant. For intransitive
# groups it is a complete invariant of the 2-closure up to rare
# 1-dim collisions; orbital graphs of a 2-closure are exactly what
# we expand. The old generator-order fingerprint was not conjugacy
# invariant and both over-generated and forced millions of
# IsConjugate(S16, ·, ·) calls (memory crash at -o 2g).
#
# Load as a library without running:
#   S16WalkAsLibrary := true; Read("scripts/s16_walk.g");
#
# Resume: if data/s16_resume/meta.g exists, continue from it and append
# to logs. Delete that directory (or set S16_FRESH := true) for a new
# descent from S_n. Optional: S16_OUTDIR, S16_MAXCLASSES, S16_SAVE_EVERY,
# S16_SAVE_MS.

LoadPackage("transgrp");;

n := 16;
T := 1556;
MAXK := 16;
if not IsBound(S16_OUTDIR) then
    S16_OUTDIR := "/home/james/src/graph-likelihood/data/";
fi;
if not IsBound(S16_SAVE_EVERY) then
    S16_SAVE_EVERY := 200;          # dump every N processed classes
fi;
if not IsBound(S16_SAVE_MS) then
    S16_SAVE_MS := 120000;          # and at least every 2 minutes
fi;
if not IsBound(S16_FRESH) then
    S16_FRESH := false;
fi;
if not IsBound(S16_MAXCLASSES) then
    S16_MAXCLASSES := 0;
fi;
S16_META := rec();
S16_TODO := [];

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

# Conjugacy-invariant (transitive: complete via transgrp). Also used
# as the 2-closure key together with OrbitalKey.
ConstituentParts := function(G, nn)
    local orbs, o, A, parts, tid, caught;
    orbs := Orbits(G, [1..nn]);
    parts := [];
    for o in orbs do
        if Length(o) = 1 then
            Add(parts, [1, 1, 0]);
        else
            A := Action(G, o);
            caught := CALL_WITH_CATCH(TransitiveIdentification, [A]);
            if caught[1] then
                tid := caught[2];
            else
                tid := -1;
            fi;
            Add(parts, [Length(o), Size(A), tid]);
        fi;
    od;
    Sort(parts);
    return parts;
end;

LabelMatrix := function(orbs, nn)
    local label, i, e;
    label := List([1..nn], x -> List([1..nn], y -> 0));
    for i in [1..Length(orbs)] do
        for e in orbs[i] do
            label[e[1]][e[2]] := i;
            label[e[2]][e[1]] := i;
        od;
    od;
    return label;
end;

# Pair-orbit *sizes* (not discovery-order labels) so the key is
# conjugacy-invariant. Transitive degree-16 classes are unique via
# transgrp; intransitive 2-closures are separated by constituent ids
# plus the 1-dim orbital-size signature.
OrbitalKey := function(G, nn)
    local orbs, label, parts, osz, sizes, rowsigs, v, w, row;
    orbs := PairOrbitalsFast(G, nn);
    osz := List(orbs, Length);
    label := LabelMatrix(orbs, nn);
    parts := ConstituentParts(G, nn);
    sizes := Collected(osz);
    rowsigs := [];
    for v in [1..nn] do
        row := [];
        for w in [1..nn] do
            if w <> v and label[v][w] > 0 then
                Add(row, osz[label[v][w]]);
            fi;
        od;
        Add(rowsigs, Collected(row));
    od;
    Sort(rowsigs);
    return Concatenation(String(Size(G)), ":c=", String(parts),
                         ":k=", String(Length(orbs)),
                         ":osz=", String(sizes),
                         ":row=", String(rowsigs));
end;

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
    local r, i, nloc, l, b, t;
    r := heap[1];
    heap[1] := heap[Length(heap)];
    Remove(heap);
    nloc := Length(heap);
    i := 1;
    while true do
        l := 2 * i;
        b := i;
        if l <= nloc and heap[l].sz > heap[b].sz then
            b := l;
        fi;
        if l + 1 <= nloc and heap[l + 1].sz > heap[b].sz then
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
        h := (h * 131 + INT_CHAR(c)) mod 1000003;
    od;
    return h + 1;
end;

NHASH := 1000003;
bins := List([1..NHASH], x -> []);

AlreadySeen := function(G)
    local key, h, bucket;
    key := OrbitalKey(G, n);
    h := HashStr(key);
    bucket := bins[h];
    if key in bucket then
        return true;
    fi;
    Add(bins[h], key);
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

StripNL := function(s)
    local nloc;
    nloc := Length(s);
    while nloc > 0 and (s[nloc] = '\n' or s[nloc] = '\r') do
        nloc := nloc - 1;
    od;
    return s{[1..nloc]};
end;

LoadSeenKeys := function(path)
    local f, line, key, h;
    bins := List([1..NHASH], x -> []);
    if not IsExistingFile(path) then
        return 0;
    fi;
    f := InputTextFile(path);
    while true do
        line := ReadLine(f);
        if line = fail then
            break;
        fi;
        key := StripNL(line);
        if Length(key) = 0 then
            continue;
        fi;
        h := HashStr(key);
        Add(bins[h], key);
    od;
    CloseStream(f);
    return Sum(bins, Length);
end;

SaveResume := function(todo, processed, nraw, nskipk, nforced, ntwinfree,
                       nfailedmax, ndup, elapsed_ms)
    local tmp, final, f, r, bucket, key, nseen;
    tmp := Concatenation(S16_OUTDIR, "s16_resume.tmp/");
    final := Concatenation(S16_OUTDIR, "s16_resume/");
    Exec(Concatenation("rm -rf ", tmp));
    Exec(Concatenation("mkdir -p ", tmp));
    PrintTo(Concatenation(tmp, "meta.g"),
            "S16_META := rec(processed:=", processed,
            ", nraw:=", nraw,
            ", nskipk:=", nskipk,
            ", nforced:=", nforced,
            ", ntwinfree:=", ntwinfree,
            ", nfailedmax:=", nfailedmax,
            ", ndup:=", ndup,
            ", elapsed_ms:=", elapsed_ms,
            ", ntodo:=", Length(todo), ");\n");
    f := OutputTextFile(Concatenation(tmp, "todo.g"), false);
    SetPrintFormattingStatus(f, false);
    PrintTo(f, "S16_TODO := [\n");
    for r in todo do
        AppendTo(f, r, ",\n");
    od;
    AppendTo(f, "];\n");
    CloseStream(f);
    nseen := 0;
    f := OutputTextFile(Concatenation(tmp, "seen.txt"), false);
    SetPrintFormattingStatus(f, false);
    for bucket in bins do
        for key in bucket do
            AppendTo(f, key, "\n");
            nseen := nseen + 1;
        od;
    od;
    CloseStream(f);
    Exec(Concatenation("rm -rf ", final, " && mv ", tmp, " ", final));
    Print("saved resume classes=", processed, " todo=", Length(todo),
          " seen=", nseen, "\n");
end;

TryLoadResume := function()
    local final, metafile, todofile, seenfile, nseen;
    final := Concatenation(S16_OUTDIR, "s16_resume/");
    metafile := Concatenation(final, "meta.g");
    todofile := Concatenation(final, "todo.g");
    seenfile := Concatenation(final, "seen.txt");
    if S16_FRESH then
        return fail;
    fi;
    if not IsExistingFile(metafile) or not IsExistingFile(todofile) then
        return fail;
    fi;
    Read(metafile);
    Read(todofile);
    nseen := LoadSeenKeys(seenfile);
    Print("RESUME classes=", S16_META.processed, " todo=", Length(S16_TODO),
          " seen=", nseen, " elapsed_s=", QuoInt(S16_META.elapsed_ms, 1000), "\n");
    return true;
end;

S16WalkMain := function()
    local S, todo, nraw, nskipk, nforced, ntwinfree, nfailedmax, ndup,
          processed, t0, G, orbs, k, ft, gens, caught, M, tmark, tseen, tmax,
          logpath, g6path, genspath, sumpath, ckptpath, saved_ms, lastsave,
          elapsed, do_save, resumed;

    logpath := Concatenation(S16_OUTDIR, "task4_s16_groups.log");
    g6path := Concatenation(S16_OUTDIR, "task4_s16_twinfree.g6");
    genspath := Concatenation(S16_OUTDIR, "task4_s16_gens.g");
    sumpath := Concatenation(S16_OUTDIR, "task4_s16_summary.txt");
    ckptpath := Concatenation(S16_OUTDIR, "task4_s16_checkpoint.txt");

    tseen := 0;
    tmax := 0;
    t0 := Runtime();
    lastsave := t0;
    resumed := TryLoadResume();

    if resumed = true then
        todo := S16_TODO;
        processed := S16_META.processed;
        nraw := S16_META.nraw;
        nskipk := S16_META.nskipk;
        nforced := S16_META.nforced;
        ntwinfree := S16_META.ntwinfree;
        nfailedmax := S16_META.nfailedmax;
        ndup := S16_META.ndup;
        saved_ms := S16_META.elapsed_ms;
        AppendTo(logpath, "# RESUME processed=", processed, " todo=",
                 Length(todo), "\n");
    else
        bins := List([1..NHASH], x -> []);
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
        saved_ms := 0;
        tmark := Runtime();
        AlreadySeen(S);
        tseen := tseen + Runtime() - tmark;
    fi;

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
        elapsed := saved_ms + Runtime() - t0;
        if RemInt(processed, 1) = 0 and processed <= 20 then
            Print("classes=", processed, " todo=", Length(todo), " raw=", nraw,
                  " twinfree=", ntwinfree, " forced=", nforced, " dup=", ndup,
                  " last|G|=", Size(G), " k=", k, " ft=", ft,
                  " t_s=", QuoInt(elapsed, 1000),
                  " seen_ms=", tseen, " max_ms=", tmax, "\n");
        elif RemInt(processed, 25) = 0 then
            Print("classes=", processed, " todo=", Length(todo), " raw=", nraw,
                  " twinfree=", ntwinfree, " forced=", nforced, " dup=", ndup,
                  " last|G|=", Size(G), " k=", k, " ft=", ft,
                  " t_s=", QuoInt(elapsed, 1000),
                  " seen_ms=", tseen, " max_ms=", tmax, "\n");
            PrintTo(ckptpath, "classes=", processed, " todo=", Length(todo),
                    " twinfree=", ntwinfree, " forced=", nforced,
                    " failedmax=", nfailedmax, " dup=", ndup, " raw=", nraw,
                    " t_ms=", elapsed, "\n");
        fi;
        tmark := Runtime();
        caught := CALL_WITH_CATCH(MaximalSubgroupClassReps, [G]);
        tmax := tmax + Runtime() - tmark;
        if caught[1] = false then
            nfailedmax := nfailedmax + 1;
            AppendTo(logpath, "# FAILED MaximalSubgroupClassReps |G|=", Size(G), "\n");
            Print("# FAILED maxsub |G|=", Size(G), " ", caught[2], "\n");
        else
            for M in caught[2] do
                if Size(M) < T then
                    continue;
                fi;
                tmark := Runtime();
                if AlreadySeen(M) then
                    tseen := tseen + Runtime() - tmark;
                    ndup := ndup + 1;
                else
                    tseen := tseen + Runtime() - tmark;
                    HeapPush(todo, PackG(M));
                fi;
            od;
        fi;
        do_save := (processed <= 20) or
                   (Runtime() - lastsave >= S16_SAVE_MS) or
                   (RemInt(processed, S16_SAVE_EVERY) = 0 and
                    Runtime() - lastsave >= 60000);
        if S16_MAXCLASSES > 0 and processed >= S16_MAXCLASSES then
            do_save := true;
        fi;
        if do_save then
            SaveResume(todo, processed, nraw, nskipk, nforced, ntwinfree,
                       nfailedmax, ndup, saved_ms + Runtime() - t0);
            lastsave := Runtime();
        fi;
        if S16_MAXCLASSES > 0 and processed >= S16_MAXCLASSES then
            Print("STOP S16_MAXCLASSES=", S16_MAXCLASSES, "\n");
            return;
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
            "time_ms=", saved_ms + Runtime() - t0, "\n");
    AppendTo(logpath, "# classes=", processed, " twin_force=", nforced,
             " twin_free=", ntwinfree, " skip=", nskipk, " raw=", nraw,
             " failed_maxsub=", nfailedmax, " dup=", ndup,
             " time_ms=", saved_ms + Runtime() - t0, "\n");
    Print("DONE classes=", processed, " raw=", nraw, " forced=", nforced,
          " twinfree=", ntwinfree, " failedmax=", nfailedmax, " dup=", ndup,
          " t_s=", QuoInt(saved_ms + Runtime() - t0, 1000), "\n");
    Exec(Concatenation("rm -rf ", S16_OUTDIR, "s16_resume ",
                       S16_OUTDIR, "s16_resume.tmp"));
end;

if not IsBound(S16WalkAsLibrary) then
    S16WalkMain();
fi;
