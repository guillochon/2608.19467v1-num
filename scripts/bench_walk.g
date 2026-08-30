# Microbench: costs inside the S16 walk. Does not start the walk.
LoadPackage("transgrp");;

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

OldFingerprint := function(G, nn)
    local ol, ordh;
    ol := Collected(List(Orbits(G, [1..nn]), Length));
    if Size(G) <= 20000 then
        ordh := Collected(List(G, Order));
    else
        ordh := Collected(List(GeneratorsOfGroup(G), Order));
    fi;
    return Concatenation(String(Size(G)), ":", String(ol), ":", String(ordh));
end;

ConstituentKey := function(G, nn)
    local orbs, o, A, parts;
    orbs := Orbits(G, [1..nn]);
    parts := [];
    for o in orbs do
        if Length(o) = 1 then
            Add(parts, [1, 1, 0]);
        else
            A := Action(G, o);
            Add(parts, [Length(o), Size(A), TransitiveIdentification(A)]);
        fi;
    od;
    Sort(parts);
    return Concatenation(String(Size(G)), ":", String(parts));
end;

OrbitalKey := function(G, nn)
    local orbs, sizes;
    orbs := PairOrbitalsFast(G, nn);
    sizes := Collected(List(orbs, Length));
    return Concatenation(ConstituentKey(G, nn), ":k=", String(Length(orbs)),
                         ":osz=", String(sizes));
end;

TimeIt := function(label, nrep, fn)
    local t, i, dummy;
    t := Runtime();
    for i in [1..nrep] do
        dummy := fn();
    od;
    t := Runtime() - t;
    Print(label, " n=", nrep, " total_ms=", t, " per_ms=",
          Float(t) / nrep, "\n");
    return dummy;
end;

nn := 16;
S := SymmetricGroup(nn);

# Typical intransitive layers from the walk.
samples := [
    rec(name := "S8xS8", G := Group(Concatenation(
        GeneratorsOfGroup(SymmetricGroup([1..8])),
        GeneratorsOfGroup(SymmetricGroup([9..16]))))),
    rec(name := "S6xS5fix", G := Group(Concatenation(
        GeneratorsOfGroup(SymmetricGroup([1..6])),
        GeneratorsOfGroup(SymmetricGroup([7..11]))))),
    rec(name := "S6wrS2", G := WreathProduct(SymmetricGroup(6), SymmetricGroup(2))),
    rec(name := "AGL42", G := PrimitiveGroup(16, 11)),
    rec(name := "C2wrS8", G := WreathProduct(SymmetricGroup(2), SymmetricGroup(8))),
];

for s in samples do
    G := s.G;
    # embed into S16 if needed
    if NrMovedPoints(G) < 16 then
        # already on 1..n
    fi;
    Print("==== ", s.name, " |G|=", Size(G), " norb=",
          Length(Orbits(G, [1..nn])), " ====\n");
    Print("  oldfp=", OldFingerprint(G, nn), "\n");
    Print("  const=", ConstituentKey(G, nn), "\n");
    TimeIt(Concatenation("  oldfp ", s.name), 50, function()
        return OldFingerprint(G, nn);
    end);
    TimeIt(Concatenation("  const ", s.name), 50, function()
        return ConstituentKey(G, nn);
    end);
    TimeIt(Concatenation("  pairorb ", s.name), 20, function()
        return PairOrbitalsFast(G, nn);
    end);
    TimeIt(Concatenation("  orbkey ", s.name), 20, function()
        return OrbitalKey(G, nn);
    end);
    t := Runtime();
    maxs := MaximalSubgroupClassReps(G);
    Print("  maxsub nreps=", Length(maxs), " ms=", Runtime()-t, "\n");
    if Length(maxs) >= 2 then
        H := maxs[1];
        K := maxs[Minimum(2, Length(maxs))];
        t := Runtime();
        c := IsConjugate(S, H, K);
        Print("  IsConjugate two maxsubs ", c, " ms=", Runtime()-t, "\n");
        t := Runtime();
        c := IsConjugate(S, H, H^Random(S));
        Print("  IsConjugate conjugate copy ", c, " ms=", Runtime()-t, "\n");
    fi;
od;

# Collision demo: old fingerprint on conjugate copies
G := samples[2].G;
keys_old := [];
keys_new := [];
for i in [1..20] do
    H := G^Random(S);
    Add(keys_old, OldFingerprint(H, nn));
    Add(keys_new, ConstituentKey(H, nn));
od;
Print("oldfp unique among 20 conjugates: ", Length(Set(keys_old)), "\n");
Print("const unique among 20 conjugates: ", Length(Set(keys_new)), "\n");
QUIT;
