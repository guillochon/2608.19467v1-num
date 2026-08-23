#!/bin/bash
set -euo pipefail
ROOT=/home/james/src/graph-likelihood
cd "$ROOT"
mkdir -p data logs code
cp -f likelihood.c Makefile code/ 2>/dev/null || true
{
  echo "=== host ==="
  date -u +%Y-%m-%dT%H:%M:%SZ
  uname -a
  lscpu | sed -n '1,20p'
  grep -E 'MemTotal|MemAvailable' /proc/meminfo
  echo "=== git ==="
  git rev-parse HEAD 2>/dev/null || echo 'no git yet'
  echo "=== sha256 sources ==="
  sha256sum likelihood.c Makefile code/independent.py
} | tee logs/hardware.txt

echo "=== make ==="
make clean
make
make check

echo "=== C extra checkpoints ==="
./likelihood c5blowup 1
./likelihood c5blowup 2
./likelihood c5parts 3 3 3 3 4
./likelihood kbal 15
./likelihood bipartite 7 8

echo "=== independent python checkpoints ==="
python3 -m py_compile code/independent.py
python3 code/independent.py check
