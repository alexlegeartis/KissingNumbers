"""Exhaust all 651 rank-four periods inside one recovered 64-word kernel.

For each 64-word quotient, enumerate independent dominating sets of at most
four vertices that contain zero. Translation makes the latter condition safe.
This is an exact finite enumeration, not a claim about asymmetric hole sets.
"""
import itertools, json, time
import numpy as np
from model import C, N, P, DATA, labels_from_holes


def span_basis(words):
    span, basis = [0], []
    for x in words:
        x = int(x)
        if x in span: continue
        basis.append(x); span += [y ^ x for y in span]
    return np.array(span), basis


if __name__ == '__main__':
    started = time.monotonic()
    H0 = np.load(DATA / 'recovery.npz')['holes'][0]
    K, basis = span_basis(C[H0])
    assert len(K) == 64 and len(basis) == 6
    periods = sorted({tuple(sorted((a, b, a ^ b))) for a in range(1,64)
                      for b in range(a+1,64)})
    assert len(periods) == 651
    lookup = {int(w): i for i,w in enumerate(C)}
    nonlinear, records, total, minimum = [], [], 0, 4
    for idx, (a,b,_) in enumerate(periods):
        T = K[[not ((a & i).bit_count()%2 or (b & i).bit_count()%2) for i in range(64)]]
        assert len(T) == 16
        rep = np.min(C[:,None] ^ T, axis=1)
        reps, orbit = np.unique(rep, return_inverse=True)
        QN = orbit[N[[lookup[int(w)] for w in reps]]]
        assert len(reps) == 64 and reps[0] == 0
        closed = [sum(1 << int(v) for v in set([u, *QN[u]])) for u in range(64)]
        allbits = (1 << 64) - 1
        solutions, nonlin = 0, 0
        # 0 is selected. Pick the other holes in increasing order, independent.
        available = allbits ^ closed[0]
        def enumerate_sets(chosen, candidates, covered):
            if covered == allbits:
                yield chosen
                return
            if len(chosen) == 4: return
            while candidates:
                bit = candidates & -candidates; candidates ^= bit
                v = bit.bit_length()-1
                yield from enumerate_sets(chosen+[v], candidates & ~closed[v], covered | closed[v])
        for chosen in enumerate_sets([0], available, closed[0]):
            solutions += 1; total += 1; minimum = min(minimum, len(chosen))
            xor = 0
            for j in chosen: xor ^= int(reps[j])
            # Four cosets form a linear subspace precisely when their sum is in T.
            nonlinear_flag = len(chosen) == 4 and xor not in T
            if nonlinear_flag or len(chosen) < 4:
                H = np.isin(orbit, chosen)
                labels_from_holes(H)
                nonlin += 1
                if len(nonlinear) < 128: nonlinear.append(H)
        records.append(dict(period=idx, dual=[a,b], solutions_containing_zero=solutions,
                            nonlinear_solutions=nonlin))
        if idx % 100 == 0:
            print(json.dumps(dict(event='progress', periods=idx+1, solutions=total,
                                  nonlinear_saved=len(nonlinear), seconds=time.monotonic()-started)),flush=True)
    if nonlinear: np.save(DATA / 'nonlinear_quotient_holes.npy', nonlinear)
    report = dict(periods=651, quotient_words=64, maximum_quotient_holes=4,
                  minimum_found_quotient_holes=minimum, lifted_minimum=16*minimum,
                  solutions_containing_zero=total,
                  nonlinear_solutions=sum(r['nonlinear_solutions'] for r in records),
                  seconds=time.monotonic()-started, records=records)
    (DATA/'quotients_result.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='records'}),flush=True)
