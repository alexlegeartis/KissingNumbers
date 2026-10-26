"""Construct several compatible four-head families using free octad parities."""
import argparse
import itertools
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OLD = ROOT.parent / 'dim32_35_20261026' / 'result'
pc = lambda x: bin(x).count('1')


def kernel(rows, width):
    pivots = {}
    for row in rows:
        while row:
            bit = row.bit_length()-1
            if bit in pivots:
                row ^= pivots[bit]
            else:
                pivots[bit] = row
                break
    basis = []
    for bit in range(width):
        if bit in pivots:
            continue
        v = 1 << bit
        for p, row in sorted(pivots.items()):
            if pc(v & row) % 2:
                v ^= 1 << p
        basis.append(v)
    return basis


def span(basis):
    result = [0]
    for b in basis:
        result += [x ^ b for x in result]
    return result


def construct(baseline=None, output=None):
    old = json.loads((baseline or OLD/'certificate35.json').read_text())
    generators = old['dense_generator']
    supports = old['full_even_supports']
    extras = (0, 1, 2)
    dense_coordinates = list(range(3, 35))
    pairs = list(itertools.combinations(extras, 2))
    triple_rows = [s >> 3 for s in supports if s & 7 == 7]
    equations = [sum((pc(row & g) % 2) << j for j, g in enumerate(generators))
                 for row in triple_rows]
    coefficient_basis = kernel(equations, len(generators))
    code_basis = []
    for coefficient in coefficient_basis:
        word = 0
        for j, g in enumerate(generators):
            if coefficient >> j & 1:
                word ^= g
        code_basis.append(word)
    candidates = [c for c in span(code_basis) if 12 <= pc(c) <= 20]
    print('triple constraints', len(equations), 'kernel dimension', len(code_basis))
    for c1, c2 in itertools.combinations(candidates, 2):
        if not 12 <= pc(c1 ^ c2) <= 20:
            continue
        bodies = [0, c1, c2]
        heads = []
        for (a, b), body in zip(pairs, bodies):
            for delta in (-1, 1):
                q = [1] * 35
                for j in range(32):
                    q[j+3] = -1 if body >> j & 1 else 1
                q[a], q[b] = 3*delta, -3*delta
                heads.append(q)
        if any(2*abs(sum(x*y for x, y in zip(a, b))) > 51
               for a, b in itertools.combinations(heads, 2)):
            continue
        signed_supports = []
        for support in supports:
            parities = [pc((support >> 3) & body) % 2
                        for pair, body in zip(pairs, bodies)
                        if all(support >> j & 1 for j in pair)]
            assert len(set(parities)) <= 1
            signed_supports.append([support, parities[0] if parities else 0])
        certificate = dict(
            format='layered_kissing_v2', dimension=35,
            points=131072+128*len(supports)+2*35*34-6+12,
            comparison=old['comparison'],
            dense_generator=generators,
            dense_coset_shift=old['dense_coset_shift'],
            dense_coordinates=dense_coordinates,
            full_parity_supports=signed_supports,
            removed_root_lines=[[a, b, -1] for a, b in pairs],
            head_antipodal_lines=heads,
            construction=dict(extra_coordinates=list(extras), body_codewords=bodies,
                              coefficient_kernel_dimension=len(code_basis)))
        output = output or ROOT/'result'/'certificate35.json'
        output.parent.mkdir(exist_ok=True, parents=True)
        output.write_text(json.dumps(certificate, separators=(',', ':'))+'\n')
        print('candidate points', certificate['points'], 'body weights', [pc(c) for c in bodies])
        return certificate
    raise RuntimeError('No compatible body choices')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    construct(args.baseline, args.output)
