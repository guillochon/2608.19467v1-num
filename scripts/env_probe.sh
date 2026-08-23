#!/bin/bash
set -euo pipefail
cd /home/james/src/graph-likelihood
date -u +%Y-%m-%dT%H:%M:%SZ
uname -a
echo '--- cpu ---'
lscpu | sed -n '1,25p'
echo '--- mem ---'
grep -E 'MemTotal|MemAvailable' /proc/meminfo
echo '--- tools ---'
python3 --version
python3 - <<'PY'
try:
    import gmpy2
    print('gmpy2', gmpy2.version())
except Exception as e:
    print('gmpy2 unavailable:', e)
PY
which gcc make git python3
gcc --version | head -1
command -v geng || true
command -v dreadnaut || true
command -v nauty-geng || true
ls -la
