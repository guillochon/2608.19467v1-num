#!/bin/bash
set -euo pipefail
cd /home/james/src/graph-likelihood
mkdir -p data logs
echo "=== subset vs rec6 at n=16,18 ==="
./likelihood bipartite 8 8
./likelihood kbal 16
./likelihood bipartite 9 9
./likelihood kbal 18

echo "=== C tower m=1..12 (mpq grid) ==="
/usr/bin/time -f 'elapsed %e sec maxrss %M KB' bash -c '
for m in $(seq 1 12); do
  echo -n "m=$m "
  ./likelihood c5blowup $m
done
' | tee data/task2_C_tower_m1_12.txt

echo "=== independent tower m=1..12 ==="
python3 code/independent.py tower --mmin 1 --mmax 12 --out data/task2_py_tower_m1_12.json

echo "=== C tower m=13..16 ==="
/usr/bin/time -f 'elapsed %e sec maxrss %M KB' bash -c '
for m in $(seq 13 16); do
  echo -n "m=$m "
  ./likelihood c5blowup $m
done
' | tee data/task2_C_tower_m13_16.txt
