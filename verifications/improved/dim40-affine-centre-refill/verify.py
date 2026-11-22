"""Verify the rank-40 affine centre and the individual shell additions."""

from collections import Counter
from importlib.util import module_from_spec, spec_from_file_location
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
AFFINE = HERE.parent / "dim40-affine-shells" / "source"
GLUING = HERE.parent / "dim41-46-p48-heads" / "coordinated45"
CORE_COUNT = 1084590


def load_module(path):
    spec = spec_from_file_location(path.stem, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def witness():
    data = json.loads((HERE / "certificate.json").read_text())
    assert data["dimension"] == 40
    assert type(data["centre"]) is int and 0 < data["centre"] < 6561
    for key in ("refill", "norm36_refill"):
        assert all(len(row) == 40 and all(type(x) is int and abs(x) < 128 for x in row)
                   for row in data[key])
    with np.load(HERE / "auxiliary.npz", allow_pickle=False) as source:
        assert set(source.files) == {"norm36", "norm30"}
        aux = [source["norm36"], source["norm30"]]
    assert all(rows.dtype == np.int8 and rows.ndim == 2 and rows.shape[1] == 40
               for rows in aux)
    return (data, aux, np.array(data["refill"], dtype=np.int64),
            np.array(data["norm36_refill"], dtype=np.int64))


def verify():
    assert __debug__, "Run without Python -O"
    data, aux, refill, extra = witness()
    centre = data["centre"]
    gram = np.array(json.loads((AFFINE / "gram.json").read_text()), dtype=np.int64)
    pinned = np.load(GLUING / "body_complement.npz", allow_pickle=False)
    assert np.array_equal(gram, pinned["q"])
    gluing = load_module(GLUING / "verify_inherited.py")
    geometry, residues, mu, tiny = gluing.geometry()
    digits, powers = geometry["digits"], geometry["powers"]
    assert mu[centre] == 6
    expected = np.unique(np.r_[centre, ((digits[centre] + digits[mu == 2]) % 3) @ powers])
    medium = expected[mu[expected] == 6]
    shortest = expected[tiny[expected] > 0]
    assert len(medium) == 82 and len(shortest) == 72
    assert Counter(mu[shortest]) == {8: 72}

    lookup = {tuple(row): i for i, row in enumerate(residues)}
    for v, norm, wanted in zip(aux, (36, 30), (
        {int(c): 495 for c in medium},
        {int(c): int(tiny[c]) for c in shortest},
    )):
        assert len(np.unique(v, axis=0)) == len(v)
        assert np.all((v @ gram * v).sum(axis=1) == norm)
        assert Counter(lookup[tuple(row)] for row in v % 3) == wanted

    # In units where body squared norms are r, the gluing bound on the
    # doubled product is r+s-18+mu(c-d). Distinct bodies in the same class
    # differ by a nonzero rank-40 lattice vector.
    classes = np.r_[centre, medium, shortest]
    radius = np.repeat([18, 12, 10], [1, len(medium), len(shortest)])
    delta = ((digits[classes, None] - digits[classes]) % 3) @ powers
    doubled = radius[:, None] + radius - 18 + mu[delta]
    assert np.all((doubled <= 0) | (doubled**2 <= radius[:, None] * radius))

    assert refill.shape == (10, 40)
    assert len(np.unique(refill, axis=0)) == len(refill)
    assert np.all((refill @ gram * refill).sum(axis=1) == 6)
    dots = refill @ gram @ refill.T
    assert int(dots[np.triu_indices(len(refill), 1)].max()) <= 3
    for v, norm in zip(aux, (36, 30)):
        products = v @ gram @ refill.T
        assert np.all((products <= 0) | (4 * products**2 <= 6 * norm))

    # For x in the core and y in this refill, xQy is a multiple of three.
    # A forbidden product is therefore at least 12. Then x-3y lies in
    # the centre class and has norm at most 36. The class minimum bound
    # is 36, so all possible conflicts come from the 495 centre bodies.
    assert not np.any((residues[centre] @ gram) % 3)
    first = aux[0][(aux[0] % 3 == residues[centre]).all(axis=1)]
    assert first.shape == (495, 40)
    assert not np.any(first @ gram @ refill.T == -6)
    assert extra.shape == (3, 40)
    assert len(np.unique(extra, axis=0)) == len(extra)
    assert np.all((extra @ gram * extra).sum(axis=1) == 36)
    owner = lookup[tuple(extra[0] % 3)]
    assert owner == 3006 and owner not in classes
    assert np.all(extra % 3 == residues[owner])
    assert int((extra @ gram @ extra.T)[np.triu_indices(len(extra), 1)].max()) <= 18
    for v, norm in zip([*aux, refill], (36, 30, 6)):
        products = extra @ gram @ v.T
        assert np.all((products <= 0) | (4 * products**2 <= 36 * norm))

    # If a norm-36 addition y conflicts with a core row x, the product
    # is at least 23, so x-y has squared norm at most 44. Its class
    # minimum is 30 and all norms in the class are congruent modulo 18.
    # This leaves only the norm-30 difference shell.
    difference = int(((digits[centre] - digits[owner]) % 3) @ powers)
    assert mu[difference] == 8
    assert not np.any((residues[difference] @ gram) % 3)
    raw = np.load(GLUING / "tiny_bodies.npz", allow_pickle=False)["H"].astype(np.int64)
    lifted = raw @ pinned["b"].T @ gluing.inverse(gram, 3)
    assert not np.any(lifted % 288)
    lifted //= 288
    assert np.array_equal(lifted @ pinned["b"], raw)
    differences = lifted[(lifted % 3 == residues[difference]).all(axis=1)]
    assert len(differences) == tiny[difference] > 0
    assert not np.any(extra @ gram @ differences.T == -6)

    total = CORE_COUNT + sum(map(len, aux)) + len(refill) + len(extra)
    assert total == data["total"] == 1126327
    return {
        "status": "PASS",
        "dimension": 40,
        "points": total,
        "centre": centre,
        "core": CORE_COUNT,
        "medium": len(aux[0]),
        "tiny": len(aux[1]),
        "minimum_shell_refill": len(refill),
        "norm36_refill": len(extra),
        "family_angle_checks": len(classes) * (len(classes) + 1) // 2,
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
