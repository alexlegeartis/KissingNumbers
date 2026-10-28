"""Malformed generators and supports must fail in both independent proofs."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def module(name):
    spec = importlib.util.spec_from_file_location('z8_'+name, ROOT/'source'/(name+'.py'))
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def main():
    if not __debug__:
        raise SystemExit('Run without -O.')
    first, second = module('replay'), module('verify_independent')
    with np.load(ROOT/'data/supports.npz',allow_pickle=False) as z:
        original = dict(z)
    index = int(np.flatnonzero(original['dual_witnesses'][:,1])[0])
    sample = {k:v.copy() for k,v in original.items()}
    for key in ('supports_low','supports_high','dual_witnesses'):
        sample[key] = sample[key][[index]]
    checkers = (first.verify_supports, second.support_proof)
    for checker in checkers:
        checker(sample)
    cases = {}
    missing = {k:v.copy() for k,v in sample.items()}
    missing['dual_witnesses'][0,1] = 0
    cases['missing_trace_check'] = missing
    overlap = {k:v.copy() for k,v in sample.items()}
    overlap['dual_witnesses'][0,1] = overlap['dual_witnesses'][0,0]
    cases['overlapping_trace_checks'] = overlap
    wrong = {k:v.copy() for k,v in sample.items()}
    mask = int(wrong['supports_low'][0]) | (int(wrong['supports_high'][0]) << 64)
    d = int(wrong['dual_witnesses'][0,0])
    row = next(r for r in range(24) if (d >> (2*r)) & 3 and bin((mask >> (4*r)) & 15).count('1') == 2)
    mask ^= 15 << (4*row)
    wrong['supports_low'][0] = mask & (2**64-1)
    wrong['supports_high'][0] = mask >> 64
    cases['wrong_trace_parity'] = wrong
    duplicate = {k:v.copy() for k,v in sample.items()}
    for key in ('supports_low','supports_high','dual_witnesses'):
        duplicate[key] = np.repeat(duplicate[key],2,axis=0)
    cases['duplicate_support'] = duplicate
    rejected = {}
    for name, candidate in cases.items():
        rejected[name] = []
        for checker in checkers:
            try:
                checker(candidate)
            except AssertionError:
                rejected[name].append(checker.__module__)
            else:
                raise AssertionError(f'{checker.__module__} accepted {name}')
    raw = json.loads((ROOT/'data/z8_generator.json').read_text())
    bad = deepcopy(raw)
    bad['generator'][0][1] = (bad['generator'][0][1]+4) % 8
    rejected['incorrect_ring_lift'] = []
    for checker in (first.verify_z8, second.parent_proof):
        try:
            checker(bad)
        except AssertionError:
            rejected['incorrect_ring_lift'].append(checker.__module__)
        else:
            raise AssertionError('An incorrect stored ring lift was accepted')
    print(json.dumps(dict(status='PASS', positive_support_controls=2,
                         rejected_cases=rejected),indent=2))


if __name__ == '__main__':
    main()
