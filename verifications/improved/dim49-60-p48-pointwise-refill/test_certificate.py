"""Semantic mutation controls for the final head certificates and file manifest."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import tempfile
import numpy as np
import verify
import verify_heads

D = Path(__file__).resolve().parent


def main():
    if not __debug__:
        raise SystemExit('Run without -O')
    base, gram, _ = verify_heads.source(D / 'source')
    original = D / 'heads49.npz'
    valid = verify_heads.certificate(original, base, gram, expected_classes=1)
    with np.load(original, allow_pickle=False) as z:
        pristine = {name: z[name].copy() for name in z.files}
    passed = []
    with tempfile.TemporaryDirectory(prefix='p48-negative-controls-') as directory:
        directory = Path(directory)
        def reject(name, change, expected_classes=1):
            data = {k: v.copy() for k, v in pristine.items()}
            change(data)
            path = directory / (name + '.npz')
            np.savez_compressed(path, **data)
            try:
                verify_heads.certificate(path, base, gram, expected_classes)
            except (ValueError, AssertionError, OverflowError):
                passed.append(name)
            else:
                raise AssertionError('Malformed certificate accepted: ' + name)
        reject('wrong_norm', lambda z: z['extra_points'].__setitem__((0, 0), int(z['extra_points'][0, 0]) + 100))
        reject('unsafe_integer', lambda z: z['extra_points'].__setitem__((0, 0), 2**40))
        reject('nonzero_padding', lambda z: z['retained_bits'].__setitem__((0, -1), int(z['retained_bits'][0, -1]) | 1))
        reject('nonisometric_transform', lambda z: z['transforms'].__setitem__((0, 0, 0), int(z['transforms'][0, 0, 0]) + 1))
        reject('wrong_partition_end', lambda z: z['extra_offsets'].__setitem__(-1, len(z['extra_points']) + 1))
        reject('wrong_tail_class_count', lambda z: None, expected_classes=2)
        reject('changed_gram', lambda z: z['gram'].__setitem__((0, 0), int(z['gram'][0, 0]) + 1))
        reject('changed_base', lambda z: z['base'].__setitem__((0, 0), int(z['base'][0, 0]) + 1))
        def duplicate(z):
            z['extra_points'] = np.concatenate((z['extra_points'], z['extra_points'][:1]))
            z['extra_offsets'][-1] += 1
        reject('duplicate_added_line', duplicate)
        def duplicate_retained(z):
            mask = np.unpackbits(z['retained_bits'], axis=1)[0, :len(base)]
            row = base[np.flatnonzero(mask)[0]]
            z['extra_points'] = np.concatenate((z['extra_points'], row[None, :]))
            z['extra_offsets'][-1] += 1
        reject('retained_extra_collision', duplicate_retained)
        raw = b'fixed proof input\n'
        (directory / 'input.bin').write_bytes(raw)
        (directory / 'MANIFEST.json').write_text(json.dumps({'files': {'input.bin': {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}}}))
        assert verify.verify_manifest(directory) == 1
        (directory / 'input.bin').write_bytes(b'changed proof input\n')
        try:
            verify.verify_manifest(directory)
        except ValueError:
            passed.append('file_hash_corruption')
        else:
            raise AssertionError('Corrupted input accepted')
    run = subprocess.run([sys.executable, '-O', str(D / 'verify.py')], capture_output=True, text=True)
    assert run.returncode != 0 and 'Run without -O' in run.stdout + run.stderr
    passed.append('disabled_assertions_rejected')
    report = dict(status='PASS', valid_control_sha256=valid['certificate_sha256'],
                  rejection_controls=passed, count=len(passed),
                  scope='Controls demonstrate rejection paths; the full positive replay remains required.')
    print(json.dumps(report, indent=2))
    return report


if __name__ == '__main__':
    main()
