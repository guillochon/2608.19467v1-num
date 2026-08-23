#!/bin/bash
# Environment setup + downloads. Run from repo root.
set -euo pipefail
cd /home/james/src/graph-likelihood
mkdir -p data/census logs /tmp/nauty-build

echo "=== python ==="
python3 --version

echo "=== nauty tools ==="
NAUTY=/home/james/src/nauty2_8_9
if [[ ! -x "$NAUTY/labelg" ]] || [[ ! -x "$NAUTY/countg" ]]; then
  echo "building nauty gtools (labelg/countg)..."
  (cd "$NAUTY" && make -j"$(nproc)" labelg countg complg shortg 2>&1 | tail -30)
fi
ls -l "$NAUTY/labelg" "$NAUTY/countg" "$NAUTY/geng" "$NAUTY/listg" 2>/dev/null || true
export PATH="$NAUTY:$PATH"

echo "=== gap ==="
if ! command -v gap >/dev/null 2>&1; then
  echo "GAP not in PATH; trying apt..."
  if sudo -n true 2>/dev/null; then
    sudo apt-get update -qq
    sudo apt-get install -y gap gap-transgrp gap-libs 2>&1 | tail -20
  else
    echo "NO_PASSWORDLESS_SUDO"
    apt-cache search gap-transgrp || true
    apt-cache search transgrp || true
  fi
fi
command -v gap || true
gap -q -c 'Print(GAPInfo.Version,"\n"); Print(NrTransitiveGroups(16),"\n"); QUIT;' 2>/dev/null || true

echo "=== download Holt-Royle VT census n=16,18 ==="
mkdir -p data/census
cd data/census
if [[ ! -f alltrans16.tar ]]; then
  curl -L --fail -o alltrans16.tar \
    "https://zenodo.org/records/4010122/files/alltrans16.tar?download=1"
fi
if [[ ! -f alltrans18.tar ]]; then
  curl -L --fail -o alltrans18.tar \
    "https://zenodo.org/records/4010122/files/alltrans18.tar?download=1"
fi
ls -l alltrans16.tar alltrans18.tar
if [[ ! -d n16 ]]; then
  mkdir -p n16 n18
  tar xf alltrans16.tar -C n16
  tar xf alltrans18.tar -C n18
  gunzip -f n16/*.gz n18/*.gz || true
fi
echo "n16 files:"; ls n16 | head
echo "n16 line counts:"; wc -l n16/* 2>/dev/null | tail -5
echo "n18 files:"; ls n18 | head
echo "DONE setup"
