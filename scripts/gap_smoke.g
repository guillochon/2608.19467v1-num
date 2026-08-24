LoadPackage("transgrp");;
Print("GAP ", GAPInfo.Version, "\n");
Print("NrTransitiveGroups(16)=", NrTransitiveGroups(16), "\n");
Print("NrTransitiveGroups(18)=", NrTransitiveGroups(18), "\n");
S := SymmetricGroup(16);;
Print("computing MaximalSubgroupClassReps(S16)...\n");
t := Runtime();
maxs := MaximalSubgroupClassReps(S);
Print("nreps=", Length(maxs), " time_ms=", Runtime()-t, "\n");
for M in maxs do
  Print("  |M|=", Size(M), " trans=", IsTransitive(M, [1..16]),
        " norb=", Length(Orbits(M, [1..16])), "\n");
od;
Print("SMOKE_OK\n");
QUIT;
