"""Self-contained exact verification of the signed-octad kissing constructions.

Python 3.9+, standard library only. Run `python verify.py` beside the certificates.
Base coordinates have squared norm 32; additional heavy-pair heads have
squared norm n+16. Normalize each vector by its own length.
No search result, floating-point tolerance or external code is trusted.
"""
import hashlib
import itertools
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def pc(x):
    return bin(x).count('1')


def verify(path):
    started = time.monotonic()
    raw = path.read_bytes()
    c = json.loads(raw)
    n = c['dimension']
    require(type(n) is int and 32 <= n <= 35, 'Dimension')
    require(c['format'] == 'signed_octad_kissing_v1', 'Format')
    generators = c['dense_generator']
    require(len(generators) == 17 and all(type(g) is int and 0 < g < 1 << 32 for g in generators), 'Generator format')
    words = {0}
    for g in generators:
        words |= {w ^ g for w in words}
    require(len(words) == 131072, 'Dense code dimension')
    distance = min(pc(w) for w in words if w)
    require(distance >= 8, 'Dense code minimum distance')
    shift = c.get('dense_coset_shift', 0)
    require(type(shift) is int and 0 <= shift < 1 << 32, 'Coset shift')
    full = c['full_even_supports']
    lines = [tuple(row) for row in c['exceptional_antipodal_lines']]
    require(len(set(full)) == len(full), 'Repeated full support')
    require(all(type(m) is int and 0 <= m < 1 << n and pc(m) == 8 for m in full), 'Full support format')
    require(len(set(lines)) == len(lines), 'Repeated exceptional line')
    require(all(len(row) == 2 for row in lines), 'Line format')
    for m, s in lines:
        require(type(m) is int and 0 <= m < 1 << n and pc(m) == 8, 'Line support')
        require(type(s) is int and s >= 0 and s & m == s, 'Line signs')
        require(s & (m & -m) == 0, 'Antipodal representative is not canonical')
    exceptional_supports = {m for m, s in lines}
    require(not exceptional_supports.intersection(full), 'Middle-layer duplicates')
    full_max = max((pc(m & k) for m, k in itertools.combinations(full, 2)), default=0)
    require(full_max <= 4, 'Full-support conflict')
    cross_max = max((pc(m & k) for m in exceptional_supports for k in full), default=0)
    require(cross_max <= 4, 'Exceptional/full conflict')
    signed_max = 0
    for (m, s), (k, t) in itertools.combinations(lines, 2):
        overlap = m & k
        dot = pc(overlap) - 2*pc(overlap & (s ^ t))
        signed_max = max(signed_max, abs(dot))
        require(abs(dot) <= 4, 'Exceptional signed conflict')
    # All points have norm squared 32.  Distinct vectors within a full octad
    # differ on >=2 signs, hence their inner product is <=32-8*2=16.
    # Full signed octads and exceptional lines have entries +/-2, roots +/-4.
    bounds = dict(dense_dense=32-2*distance, full_same_support=16,
                  full_different_supports=4*full_max, exceptional_full=4*cross_max,
                  exceptional_exceptional=4*signed_max, antipodal_pair=-32,
                  roots_roots=16, dense_middle=16, dense_roots=8, middle_roots=16)
    require(max(bounds.values()) <= 16, 'Pair bound')
    extension = c.get('heavy_pair_extension')
    head_checks = None
    if extension is not None:
        a, b = extension['coordinates']
        require(32 <= a < b < n, 'Heavy coordinates must lie outside the dense support')
        require(extension['removed_root_antipodal_lines'] == 1 and
                extension['added_head_antipodal_lines'] == 2 and
                extension['head_squared_norm'] == n+16, 'Heavy-pair metadata')
        heads = []
        for body_sign in (-1, 1):
            for tail_sign in (-1, 1):
                q = [body_sign] * n
                q[a], q[b] = 3*tail_sign, -3*tail_sign
                heads.append(q)
        require(all(sum(v*v for v in q) == n+16 for q in heads), 'Head norms')
        cube_weights = [pc(w ^ shift) for w in words]
        cube_max = max(abs(32-2*w) for w in cube_weights)
        middle_max = 0
        for q in heads:
            for m in full:
                entries = [q[j] for j in range(n) if m >> j & 1]
                value = sum(abs(v) for v in entries)
                if sum(v < 0 for v in entries) % 2:
                    value -= 2*min(abs(v) for v in entries)
                middle_max = max(middle_max, 2*value)
            for m, s in lines:
                value = 2*sum(q[j]*(-1 if s >> j & 1 else 1) for j in range(n) if m >> j & 1)
                middle_max = max(middle_max, abs(value))
        root_max = 0
        removed_count = 0
        for i, j in itertools.combinations(range(n), 2):
            for si, sj in itertools.product((-1, 1), repeat=2):
                if (i, j) == (a, b) and si == -sj:
                    removed_count += 1
                    continue
                root_max = max(root_max, *(4*(si*q[i]+sj*q[j]) for q in heads))
        require(removed_count == 2, 'Root deletion count')
        for value in (cube_max, middle_max, root_max):
            require(4*value*value <= 32*(n+16), 'Head/base normalized inner product')
        head_max = max(sum(x*y for x, y in zip(q, r)) for q, r in itertools.combinations(heads, 2))
        require(2*head_max <= n+16, 'Head/head normalized inner product')
        head_checks = dict(head_count=4, removed_roots=2, head_squared_norm=n+16,
                           dense_coset_minimum_weight=min(cube_weights), dense_coset_maximum_weight=max(cube_weights),
                           head_dense_dot_max=cube_max, head_middle_dot_max=middle_max,
                           head_root_dot_max=root_max, head_head_dot_max=head_max)
    count = len(words) + 128*len(full) + 2*len(lines) + 2*n*(n-1) + (2 if extension else 0)
    require(type(c['points']) is int and count == c['points'], 'Point count')
    baseline = c['comparison']['kissing_points']
    require(count > baseline, 'No improvement over supplied benchmark')
    return dict(dimension=n, verified_points=count, benchmark=baseline, gain=count-baseline,
                base_squared_norm=32, max_normalized_inner_product='1/2',
                full_supports=len(full), exceptional_lines=len(lines),
                dense_code_size=len(words), dense_code_minimum_distance=distance,
                base_pair_inner_product_upper_bounds=bounds, heavy_pair_checks=head_checks,
                certificate_sha256=hashlib.sha256(raw).hexdigest(),
                seconds=time.monotonic()-started, all_checks_pass=True)


def main():
    paths = [Path(x) for x in sys.argv[1:]] or sorted(ROOT.glob('certificate*.json'))
    require(bool(paths), 'No certificates')
    results = [verify(path) for path in paths]
    for result in results:
        print(json.dumps(result, sort_keys=True))
    recorded = [{k: v for k, v in r.items() if k != 'seconds'} for r in results]
    (ROOT/'verification.json').write_text(json.dumps(recorded, indent=2)+'\n')


if __name__ == '__main__':
    main()
