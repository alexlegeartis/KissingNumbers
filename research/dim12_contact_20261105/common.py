"""Shared paths and single-threaded numerical setup for the contact campaign."""
import os
for key in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    os.environ[key] = '1'
from pathlib import Path
import sys
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parent
OLD = ROOT.parent / 'dim12_quaternionic_20261104'
AUDIT = ROOT.parent / 'dim12_paper_audit_20261104'
sys.path.insert(0, str(OLD))
from search import normalize, metrics, relax
from continuation import optimize
from verify import verify

def save_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2) + '\n')
    tmp.replace(path)

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def reserve():
    if 127 in os.sched_getaffinity(0):
        raise RuntimeError('Run with taskset -c 0-126; reserve CPU 127.')
