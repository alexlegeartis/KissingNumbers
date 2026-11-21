#!/usr/bin/env python3
"""Reproduce the dimension-88 re-marked chord and sparse construction."""
from fractions import Fraction as Q
from itertools import combinations, product
from pathlib import Path
import importlib.util
import json
import sys

HERE = Path(__file__).resolve().parent
V = HERE.parent
sys.path.insert(0, str(V / "dim86-87-axial-pless"))
import axial_pless as axial

def local_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

chord = local_module("chord", V / "dim86-87-ternary-chord-sections/verify.py")
sqs = local_module("sqs", V / "dim90-91-a2-rotated-square/support.py")
FIXED = [0, 3]
AXES = [1, 2, 24, 28]
FULL = 42
DIM = 88
PREVIOUS = 13029906432

def no_occupation(W, i):
    return (((W[:, 2*(i//24)] | W[:, 2*(i//24)+1]) >> axial.coordinate_bit(i)) & 1) == 0

def results():
    assert __debug__, "Run without -O"
    assert chord.hashlib.sha256(chord.QR.read_bytes()).hexdigest() == chord.SOURCES[chord.QR]
    G = chord.np.asarray(json.loads(chord.QR.read_text())["generator_Q48"], dtype="i8")
    assert G.shape == (24, 48) and set(G.flat) == {0, 1, 2}
    assert chord.rank3(G) == 24 and (G @ G.T % 3 == 0).all()
    assert len(set(FIXED + AXES)) == len(FIXED) + len(AXES) == 6
    assert chord.rank3(G[:, FIXED + AXES]) == 6
    assert 48 - len(FIXED + AXES) == FULL
    norm = FULL + len(AXES)*Q(3,4)
    assert norm == 45 and 2*FULL + len(AXES) == DIM
    assert all(-norm/2 <= norm - Q(3,2)*d <= norm/2 for d in range(15, 47, 3))
    dense = 2 * 2**len(AXES) * 3**(24 - len(FIXED + AXES))

    counts = axial.phase_counts()
    axial.verify_geometry()
    _, W = axial.build_outer()
    for i in FIXED:
        W = W[no_occupation(W, i)]
    H = axial.census(W, AXES)
    assert len(W) == 194304 and H.tolist() == [37468, 80272, 58152, 16784, 1628]
    sparse = int(H @ counts[:len(H)])

    # Induced Steiner quadruple supports on the full triangle planes.
    supports = [row for row in sqs.supports() if max(row) < FULL]
    triples = [triple for block in supports for triple in combinations(block, 3)]
    assert len(triples) == len(set(triples))
    phases = [s for s in range(256) if s.bit_count() % 2 == 0]
    assert len(phases) == 128 and all((x ^ y).bit_count() >= 2 for x, y in combinations(phases, 2))
    squares = len(supports) * len(phases)
    deltas = [min((a-b)%360, (b-a)%360)
              for a,b in product(range(15,360,90), range(0,360,60))]
    assert min(deltas) == 15 and 12**2*3 < (45-24)**2  # square/Pless
    assert 4**2 < 45  # square/QR
    assert 2*3 < norm/4  # D88 root/QR or Pless by Cauchy--Schwarz
    assert 15 - len(AXES) > 4 and FULL > 4  # support sizes separate the layers
    # D88 roots have exactly two signed unit entries; distinct roots have dot<=1.
    # Square coordinates have absolute value at most 1/sqrt(2) in the same frame.
    roots = 2 * DIM * (DIM - 1)
    total = dense + sparse + squares + roots
    return dict(dimension=DIM, lower_bound=total, previous=PREVIOUS, gain=total-PREVIOUS,
                dense=dense, sparse=sparse, square_supports=len(supports), squares=squares,
                roots=roots, histogram=H.tolist())

if __name__ == "__main__":
    print(json.dumps(results(), indent=2))
