"""Positive controls and deliberately invalid supports for both exact replays."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main():
    scripts = ['mixed_pair_certificate.py', 'verify_support_bits.py']
    cases = {}
    with np.load(ROOT/'source/dim61_96/mixed_pair_certificate.npz') as data:
        original = dict(data)
    duplicate = {k:v.copy() for k,v in original.items()}
    for key in ('supports_low', 'supports_high'):
        duplicate[key][1] = duplicate[key][0]
    cases['duplicate_support'] = duplicate
    overlap = {k:v.copy() for k,v in original.items()}
    mask = int(original['supports_low'][0]) | (int(original['supports_high'][0]) << 64)
    for r in range(24):
        active = [4*r+c for c in range(4) if mask & (1 << (4*r+c))]
        inactive = [4*r+c for c in range(4) if not mask & (1 << (4*r+c))]
        if len(active) == 2:
            changed = mask ^ (1 << active[0]) ^ (1 << inactive[0])
            break
    else:
        raise RuntimeError('No two-cell row for the intersection negative control')
    assert bin(changed & mask).count('1') == 23
    overlap['supports_low'][1] = changed & ((1 << 64)-1)
    overlap['supports_high'][1] = changed >> 64
    cases['intersection_23'] = overlap
    weight = {k:v.copy() for k,v in original.items()}
    removed = mask & (mask-1)
    weight['supports_low'][0] = removed & ((1 << 64)-1)
    weight['supports_high'][0] = removed >> 64
    cases['weight_23'] = weight
    report = {'status':'PASS', 'positive_controls':[], 'rejected_cases':{}}
    with tempfile.TemporaryDirectory(prefix='kissing_negative_controls_') as directory:
        source = Path(directory)/'source'
        shutil.copytree(ROOT/'source', source, ignore=shutil.ignore_patterns('__pycache__'))
        folder = source/'dim61_96'
        for script in scripts:
            subprocess.run([sys.executable, str(folder/script)], cwd=folder,
                           check=True, capture_output=True, text=True)
            report['positive_controls'].append(script)
        for name, candidate in cases.items():
            np.savez_compressed(folder/'mixed_pair_certificate.npz', **candidate)
            rejected = []
            for script in scripts:
                result = subprocess.run([sys.executable, str(folder/script)], cwd=folder,
                                        capture_output=True, text=True)
                if result.returncode == 0 or 'AssertionError' not in result.stderr:
                    raise RuntimeError(f'{script} did not reject {name} as a bad certificate')
                rejected.append(script)
            report['rejected_cases'][name] = rejected
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
