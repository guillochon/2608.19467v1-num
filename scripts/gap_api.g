# Probe GAP APIs used by the S16 walk.
Print("CHAR_INT ", IsBoundGlobal("CHAR_INT"), " CHAR_INT ", IsBoundGlobal("CHAR_INT"), "\n");
Print("WriteLine ", IsBoundGlobal("WriteLine"), "\n");
Print("BlistList ", IsBoundGlobal("BlistList"), " BlistList ", IsBoundGlobal("BlistList"), "\n");
Print("MaximalSubgroupClassReps ", IsBoundGlobal("MaximalSubgroupClassReps"), "\n");
Print("MaximalSubgroupClassReps ", IsBoundGlobal("MaximalSubgroupClassReps"), "\n");
Print("NewDictionary ", IsBoundGlobal("NewDictionary"), "\n");
S := SymmetricGroup(8);;
H := AlternatingGroup(8);;
Print("IsConjugate ", IsConjugate(S, H, H^Random(S)), "\n");
d := NewDictionary("", false);;
AddDictionary(d, "abc", [1,2]);
Print("lookup ", LookupDictionary(d, "abc"), " missing ", LookupDictionary(d, "zzz"), "\n");
Print("API_OK\n");
QUIT;
