"""Construct SQS(68) explicitly and verify the enlarged 66--68D layers."""
from itertools import combinations
import json
from math import comb
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent


def blocks68():
    q, infinity = 67, 67
    inverse = [0]+[pow(x, q-2, q) for x in range(1, q)]
    root = next(x for x in range(q) if (x*x-x+1) % q == 0)
    seed = [infinity, 0, 1, root]
    blocks = set()
    matrices = 0
    for a in range(q):
        for c in range(q):
            lead = a if a else c
            if not lead or lead > q//2:
                continue
            for k in range(q):
                if a:
                    b, d = k, (1+k*c)*inverse[a] % q
                else:
                    b, d = -inverse[c] % q, k
                assert (a*d-b*c) % q == 1
                matrices += 1
                image = []
                for x in seed:
                    num, den = (a, c) if x == infinity else ((a*x+b) % q, (c*x+d) % q)
                    image.append(num*inverse[den] % q if den else infinity)
                assert len(set(image)) == 4
                blocks.add(tuple(sorted(image)))
    assert matrices == q*(q*q-1)//2
    return sorted(blocks)


def verify_blocks(blocks):
    assert len(blocks) == len(set(blocks)) == 12529
    triples = [t for block in blocks for t in combinations(block, 3)]
    assert len(triples) == len(set(triples)) == comb(68, 3)
    incidence = np.zeros((len(blocks), 68), dtype=np.float32)
    for i, block in enumerate(blocks):
        assert len(block) == len(set(block)) == 4 and 0 <= min(block) <= max(block) < 68
        incidence[i, block] = 1
    # Independent coordinate replay; every partial sum is an integer <=68.
    with threadpool_limits(limits=3):
        for start in range(0, len(blocks), 256):
            end = min(start+256, len(blocks))
            dots = incidence[start:end] @ incidence[:end].T
            dots[np.arange(end-start), start+np.arange(end-start)] = 0
            assert np.all(dots <= 2)
    return len(triples)


def main():
    if not __debug__:
        raise SystemExit('Run without -O: assertions are required.')
    blocks = blocks68()
    triples = verify_blocks(blocks)
    parent = json.loads((HERE.parent / 'dim61_96/expanded_support_layers.json').read_text())
    dimensions = {}
    for n in (66, 67, 68):
        count = sum(max(b) < n for b in blocks)
        assert count == {66:11088, 67:11792, 68:12529}[n]
        johnson = n*((n-1)*((n-2)//2)//3)//4
        if n in (67, 68):
            assert count == johnson
        old = parent['dimensions'][str(n)]
        terms = dict(old['terms'], weight4=16*count)
        assert terms['dense64'] == 1 << 28 and terms['poles'] == 2*n
        assert terms['weight16'] == old['support16_count']*2048
        dimensions[str(n)] = dict(terms=terms, total=sum(terms.values()),
                                 previous_bound=old['total'], gain=sum(terms.values())-old['total'],
                                 support4_count=count, fixed_support_size_packing_bound=johnson)
    # Norm64 amplitudes: dense1, support16=2, support4=4, axial8.
    assert 4*4*4 == 64 and 2*4*4 == 32 and 4*4*2 == 32
    assert 1*4*4 <= 32 and 8*4 == 32
    report = dict(status='Verified: complete PSL(2,67) orbit, every triple and independent coordinate intersections',
                  design=dict(points=68, blocks=len(blocks), triples=triples, seed=[67,0,1,30]),
                  dimensions=dimensions,
                  input_scope='Dense and weight-16 layers and their code-distance proof are unchanged; all new weight-four supports are constructed and checked here.')
    (HERE / 'netto.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
