#!/usr/bin/env python3
"""Check the marked dense section and replay the axial Pless construction."""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
import json

import axial_pless as axial

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "chord_proof", HERE.parent / "dim86-87-ternary-chord-sections" / "verify.py"
)
chord = importlib.util.module_from_spec(spec)
spec.loader.exec_module(chord)

if __name__ == "__main__":
    assert __debug__, "Run without -O"
    assert chord.hashlib.sha256(chord.QR.read_bytes()).hexdigest() == chord.SOURCES[chord.QR]
    G = chord.np.asarray(json.loads(chord.QR.read_text())["generator_Q48"], dtype="i8")
    assert G.shape == (24, 48) and set(G.flat) == {0, 1, 2}
    assert chord.rank3(G) == 24 and (G @ G.T % 3 == 0).all()
    for dimension, (axes, _) in axial.AXES.items():
        assert (0 in axes) == (dimension == 87)
        assert len(set(axes) | {0}) == 9
        assert chord.rank3(G[:, sorted(set(axes) | {0})]) == 9
        assert 48 - len(set(axes) | {0}) == 39
        assert 2 * 39 + len(axes) == dimension
    norm = 39 + 8 * Q(3, 4)
    assert norm == 45
    assert axial.DENSE == 2 * 2**8 * 3**(24-9)
    assert all(-norm/2 <= norm - Q(3, 2)*d <= norm/2 for d in range(15, 48, 3))
    assert norm * (8 * Q(3, 4)) < (norm/2)**2  # dense / pure tail
    assert 15 < 39  # the sparse layer vanishes on at least one full plane
    axial.main()
