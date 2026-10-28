"""Isolated replay of the larger 61--62D improvements; optional C++ cross-check."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[2] / 'common'))
from checkpoint import check_manifest, replay, temporary_source
from final import results


def validate(source):
    path = Path('significant49_72/compressed_layers.json')
    fresh = json.loads((source/path).read_text())
    assert fresh == json.loads((ROOT/'source'/path).read_text())
    totals = {int(d): row['total'] for d, row in fresh['dimensions'].items()}
    current = {row['dimension']: row['total'] for row in results()}
    assert all(totals[d] < current[d] for d in totals)
    return totals


def cpp_crosscheck(source):
    compiler = shutil.which('g++')
    if not compiler:
        raise SystemExit('--cpp requires g++.')
    folder = source/'dim61_96'
    executable = folder/('enumerate64.exe' if sys.platform == 'win32' else 'enumerate64')
    subprocess.run([compiler, '-O3', '-mpopcnt', '-mbmi',
                    str(folder/'enumerate_binary.cpp'), '-o', str(executable)], check=True)
    sys.path.insert(0, str(source/'dim49plus'))
    from binary_62_code import build_code64_rows
    rows = build_code64_rows()
    generator = folder/'independent64.generator.txt'
    generator.write_text('64 28\n'+'\n'.join(f'{row:x} 0' for row in rows)+'\n')
    selected = folder/'independent64.weight16.bin'
    output = subprocess.check_output([str(executable), str(generator), str(selected)], text=True)
    actual = json.loads(output)
    report = json.loads((source/'significant49_72/compressed_layers.json').read_text())
    assert (actual['n'], actual['k']) == (64, 28)
    assert actual['histogram'] == report['parent_weight_histogram']
    raw = np.fromfile(selected, dtype=np.uint64).reshape(-1, 2)
    assert raw.shape == (30828, 2) and not np.any(raw[:, 1])
    mask = sum(1 << c for c in report['removed_coordinates'])
    supports = np.sort(raw[(raw[:, 0] & np.uint64(mask)) == 0, 0]).astype('<u8')
    assert len(supports) == report['supports']
    assert hashlib.sha256(supports.tobytes()).hexdigest() == report['supports_sha256']
    print('PASS: independent C++ Gray-code enumeration agrees on all 268,435,456 words and 9,673 supports.', flush=True)


def main():
    check_manifest(ROOT)
    with temporary_source(ROOT) as source:
        replay(source, 'significant49_72/compressed_layers.py')
        totals = validate(source)
        if '--cpp' in sys.argv:
            cpp_crosscheck(source)
    print(json.dumps(dict(status='PASS', dimensions=totals, cpp='--cpp' in sys.argv)))


if __name__ == '__main__':
    main()
