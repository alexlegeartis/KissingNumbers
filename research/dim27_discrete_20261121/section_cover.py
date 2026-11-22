"""Verify the fractional clique cover for low-cost imports from 28 dimensions."""
from argparse import ArgumentParser
import json

import numpy as np

from source import HERE, ROOT, keys, leech_shell, load_points


def catalogue(source):
    package = ROOT / "verifications/improved/dim28-takhanov-root-exchange/data"
    with np.load(package / "coordinates28_223528.npz") as archive:
        points = archive["points"].astype(np.int16)
    native = np.all(points[:, 24:] == 0, axis=1)
    bodies = np.unique(points[~native, :24], axis=0)
    shell = 2 * leech_shell().astype(np.int16)
    old_equator = source[np.all(source[:, 24:] == 0, axis=1), :24]
    missing = shell[~np.isin(keys(shell), keys(old_equator))]
    imports = [np.pad(missing, ((0, 0), (0, 3)))]
    # Each 28D body already coexists with every original 28D equatorial root.
    # The other equatorial roots, plus the nonequatorial source, suffice to price it.
    relevant = np.r_[source[np.any(source[:, 24:] != 0, axis=1)],
                     np.pad(old_equator[~np.isin(keys(old_equator), keys(points[native, :24]))],
                            ((0, 0), (0, 3)))]
    matrix = relevant[:, :24].astype(np.float32).T
    tail_types = np.unique(source[:, 24:], axis=0)
    required = 128 - np.sum(bodies * bodies, axis=1)
    for norm in np.unique(required):
        group = bodies[required == norm]
        tails = tail_types[np.sum(tail_types * tail_types, axis=1) == norm]
        limits = 64 - tails @ relevant[:, 24:].T
        for start in range(0, len(group), 64):
            block = group[start:start + 64]
            products = block.astype(np.float32) @ matrix
            for tail, cap in zip(tails, limits):
                keep = block[np.sum(products > cap, axis=1) <= 12]
                imports.append(np.c_[keep, np.tile(tail, (len(keep), 1))])
    rows = np.concatenate(imports)
    return rows[~np.isin(keys(rows), keys(source))]


def verify(check_catalogue):
    source = load_points().astype(np.int16)
    with np.load(HERE / "data/d27_section_clique_cover.npz") as data:
        new = data["points"]
        anchors = data["anchors"]
        weights = data["numerators"].astype(np.int64)
        denominator = int(data["denominator"])
        indices, offsets = data["indices"], data["offsets"]

    assert len(np.unique(keys(new))) == len(new)
    assert np.all(np.sum(new.astype(np.int64) ** 2, axis=1) == 128)
    assert not np.any(np.isin(keys(new), keys(source)))
    assert denominator > 0 and np.all(weights > 0)
    assert len(offsets) == len(anchors) + 1 == len(weights) + 1
    assert offsets[0] == 0 and offsets[-1] == len(indices) and np.all(np.diff(offsets) > 0)
    assert np.all((indices >= 0) & (indices < len(new)))
    by_anchor, inverse = np.unique(keys(anchors), return_inverse=True)
    assert np.all(np.isin(by_anchor, keys(source)))
    spent = np.zeros(len(by_anchor), dtype=np.int64)
    np.add.at(spent, inverse, weights)
    assert np.all(spent <= denominator)
    covered = np.zeros(len(new), dtype=np.int64)
    for i, anchor in enumerate(anchors):
        members = indices[offsets[i]:offsets[i + 1]]
        assert len(np.unique(members)) == len(members)
        points = new[members].astype(np.int64)
        assert np.all(points @ anchor > 64) and np.all(points @ points.T > 64)
        covered[members] += weights[i]
    assert np.all(covered >= denominator)

    if check_catalogue:
        expected = catalogue(source)
        assert np.array_equal(np.sort(keys(new)), np.sort(keys(expected)))
    return dict(candidates=len(new), deleted_roots=int(np.all(new[:, 24:] == 0, axis=1).sum()),
                cliques=len(anchors), denominator=denominator, maximum_union_size=len(source),
                catalogue_rebuilt=check_catalogue)


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--rebuild-catalogue", action="store_true")
    print(json.dumps(verify(parser.parse_args().rebuild_catalogue), indent=2))
