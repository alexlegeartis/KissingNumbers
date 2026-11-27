#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python tail26.py
python real26.py
python arc26.py
python pin_angles.py
python relax26.py
