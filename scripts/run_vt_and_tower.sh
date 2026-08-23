#!/bin/bash
set -euo pipefail
cd /home/james/src/graph-likelihood
mkdir -p data/census logs
NAUTY=/home/james/src/nauty2_8_9
export PATH="$NAUTY:$PATH"

echo "==== download Holt-Royle VT census ===="
cd data/census
if [[ ! -f alltrans16.tar ]]; then
  curl -L --fail --retry 5 -o alltrans16.tar \
    "https://zenodo.org/records/4010122/files/alltrans16.tar?download=1"
fi
if [[ ! -f alltrans18.tar ]]; then
  curl -L --fail --retry 5 -o alltrans18.tar \
    "https://zenodo.org/records/4010122/files/alltrans18.tar?download=1"
fi
ls -l alltrans16.tar alltrans18.tar
mkdir -p n16 n18
if [[ ! -f n16/.unpacked ]]; then
  tar xf alltrans16.tar -C n16
  gunzip -f n16/*.gz || true
  touch n16/.unpacked
fi
if [[ ! -f n18/.unpacked ]]; then
  tar xf alltrans18.tar -C n18
  gunzip -f n18/*.gz || true
  touch n18/.unpacked
fi
echo "n16 files $(ls n16 | wc -l) lines $(cat n16/* 2>/dev/null | grep -cv -e '^>' -e '^$' || true)"
echo "n18 files $(ls n18 | wc -l) lines $(cat n18/* 2>/dev/null | grep -cv -e '^>' -e '^$' || true)"
cd /home/james/src/graph-likelihood

echo "==== C tower timings ===="
if [[ ! -f logs/task2_C_tower_times.txt ]]; then
  /usr/bin/time -f 'C_m1_12 elapsed %e sec maxrss %M KB' bash -c '
    for m in $(seq 1 12); do echo -n "m=$m "; ./likelihood c5blowup $m; done
  ' > logs/task2_C_m1_12_rerun.txt 2> logs/task2_C_m1_12_time.txt
  /usr/bin/time -f 'C_m13_16 elapsed %e sec maxrss %M KB' bash -c '
    for m in $(seq 13 16); do echo -n "m=$m "; ./likelihood c5blowup $m; done
  ' > logs/task2_C_m13_16_rerun.txt 2> logs/task2_C_m13_16_time.txt
  cat logs/task2_C_m1_12_time.txt logs/task2_C_m13_16_time.txt > logs/task2_C_tower_times.txt
fi
cat logs/task2_C_tower_times.txt

echo "==== VT n=16 ===="
python3 scripts/vt_census.py data/census/n16 --n 16 --T 1556 \
  --labelg "$NAUTY/labelg" --g6aut code/g6aut \
  --provenance "Holt-Royle vertex-transitive census, Zenodo 10.5281/zenodo.4010122, arXiv:1811.09015; disconnected included (all valencies 0..15)" \
  --out data/task4_vt16.json

echo "==== VT n=18 ===="
python3 scripts/vt_census.py data/census/n18 --n 18 --T 462 \
  --labelg "$NAUTY/labelg" --g6aut code/g6aut \
  --provenance "Holt-Royle vertex-transitive census, Zenodo 10.5281/zenodo.4010122, arXiv:1811.09015; disconnected included (all valencies 0..17)" \
  --out data/task4_vt18.json

echo DONE
