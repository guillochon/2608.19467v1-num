#!/bin/bash
set -euo pipefail
cd /home/james/src/graph-likelihood
mkdir -p data/census logs
export NAUTY=/home/james/src/nauty2_8_9
export PATH="$NAUTY:$PATH"

echo "==== $(date -u) host ====" | tee logs/run_followup.txt
lscpu | sed -n '1,20p' | tee -a logs/run_followup.txt
grep -E 'MemTotal|MemAvailable' /proc/meminfo | tee -a logs/run_followup.txt

echo "==== nauty gtools ====" | tee -a logs/run_followup.txt
if [[ ! -x $NAUTY/labelg || ! -x $NAUTY/countg ]]; then
  make -C "$NAUTY" -j"$(nproc)" labelg countg listg shortg complg 2>&1 | tee -a logs/nauty_make.txt | tail -40
fi
ls -l "$NAUTY/labelg" "$NAUTY/countg" "$NAUTY/geng" "$NAUTY/listg" | tee -a logs/run_followup.txt

echo "==== exact T ====" | tee -a logs/run_followup.txt
python3 scripts/exact_T.py 2>&1 | tee -a logs/exact_T.txt

echo "==== task3 fit residuals ====" | tee -a logs/run_followup.txt
python3 scripts/task3_fit.py 2>&1 | tee logs/task3_fit.txt

echo "==== C tower wall times (n=80 is m=16) ====" | tee -a logs/run_followup.txt
if [[ ! -f logs/task2_C_tower_times.txt ]]; then
  /usr/bin/time -f 'C_m1_12 elapsed %e sec maxrss %M KB' bash -c '
    for m in $(seq 1 12); do echo -n "m=$m "; ./likelihood c5blowup $m; done
  ' > logs/task2_C_m1_12_rerun.txt 2> logs/task2_C_m1_12_time.txt
  /usr/bin/time -f 'C_m13_16 elapsed %e sec maxrss %M KB' bash -c '
    for m in $(seq 13 16); do echo -n "m=$m "; ./likelihood c5blowup $m; done
  ' > logs/task2_C_m13_16_rerun.txt 2> logs/task2_C_m13_16_time.txt
  cat logs/task2_C_m1_12_time.txt logs/task2_C_m13_16_time.txt | tee logs/task2_C_tower_times.txt
fi
cat logs/task2_C_tower_times.txt | tee -a logs/run_followup.txt

echo "==== Holt-Royle census download ====" | tee -a logs/run_followup.txt
cd data/census
if [[ ! -f alltrans16.tar ]]; then
  curl -L --fail --retry 5 -o alltrans16.tar \
    "https://zenodo.org/records/4010122/files/alltrans16.tar?download=1"
fi
if [[ ! -f alltrans18.tar ]]; then
  curl -L --fail --retry 5 -o alltrans18.tar \
    "https://zenodo.org/records/4010122/files/alltrans18.tar?download=1"
fi
ls -l alltrans16.tar alltrans18.tar | tee -a ../../logs/run_followup.txt
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
echo "n16 files $(ls n16 | wc -l) lines $(cat n16/* 2>/dev/null | grep -v '^>' | grep -v '^$' | wc -l)" | tee -a ../../logs/run_followup.txt
echo "n18 files $(ls n18 | wc -l) lines $(cat n18/* 2>/dev/null | grep -v '^>' | grep -v '^$' | wc -l)" | tee -a ../../logs/run_followup.txt
cd /home/james/src/graph-likelihood

echo "==== evaluate VT n=16 ====" | tee -a logs/run_followup.txt
python3 scripts/vt_census.py data/census/n16 --n 16 --T 1556 \
  --labelg "$NAUTY/labelg" --countg "$NAUTY/countg" \
  --provenance "Holt-Royle vertex-transitive census, Zenodo 10.5281/zenodo.4010122, arXiv:1811.09015; includes disconnected (all valencies 0..15)" \
  --out data/task4_vt16.json 2>&1 | tee logs/vt16_eval.txt

echo "==== GAP install attempt ====" | tee -a logs/run_followup.txt
if command -v gap >/dev/null 2>&1; then
  echo "gap already present $(command -v gap)" | tee -a logs/run_followup.txt
elif sudo -n true 2>/dev/null; then
  sudo apt-get update -qq
  sudo DEBIAN_FRONTEND=noninteractive apt-get install -y gap gap-transgrp 2>&1 | tail -30 | tee -a logs/run_followup.txt
else
  echo "no passwordless sudo; GAP not installed" | tee -a logs/run_followup.txt
fi

if command -v gap >/dev/null 2>&1; then
  echo "==== GAP S16 walk (background) + VT orbital cross-check ====" | tee -a logs/run_followup.txt
  gap -q scripts/s16_walk.g > logs/s16_walk.stdout 2> logs/s16_walk.stderr &
  echo "s16_walk pid $!" | tee -a logs/run_followup.txt
  gap -q scripts/vt16.g > logs/vt16_gap.stdout 2> logs/vt16_gap.stderr &
  echo "vt16.g pid $!" | tee -a logs/run_followup.txt
fi

echo "==== n=16 VT eval done; n=18 starts after ====" | tee -a logs/run_followup.txt
python3 scripts/vt_census.py data/census/n18 --n 18 --T 462 \
  --labelg "$NAUTY/labelg" --countg "$NAUTY/countg" \
  --provenance "Holt-Royle vertex-transitive census, Zenodo 10.5281/zenodo.4010122, arXiv:1811.09015; includes disconnected (all valencies 0..17)" \
  --out data/task4_vt18.json 2>&1 | tee logs/vt18_eval.txt

echo "ALL FOREGROUND STEPS DONE $(date -u)" | tee -a logs/run_followup.txt
