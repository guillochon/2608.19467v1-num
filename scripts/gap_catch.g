Print("CALL_WITH_CATCH ", IsBoundGlobal("CALL_WITH_CATCH"), "\n");
S := SymmetricGroup(4);;
A := AlternatingGroup(4);;
B := A^((1,2,3,4));
Print("IsConjugate subgroups ", IsConjugate(S, A, B), "\n");
r := CALL_WITH_CATCH(MaximalSubgroupClassReps, [S]);
Print("catch ok=", r[1], " nreps=", Length(r[2]), "\n");
Print("OK\n");
QUIT;
