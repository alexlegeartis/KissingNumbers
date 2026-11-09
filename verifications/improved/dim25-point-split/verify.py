"""Check every package hash and replay the exact proof in an isolated directory."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

from final import results

ROOT = Path(__file__).resolve().parent


def check_manifest():
    manifest = ROOT / 'data/SHA256SUMS.txt'
    listed = set()
    for line in manifest.read_text().splitlines():
        digest, name = line.split('  ', 1)
        path = (ROOT / name).resolve()
        if (ROOT not in path.parents or path == manifest or name in listed
                or hashlib.sha256(path.read_bytes()).hexdigest() != digest):
            raise ValueError('Invalid package manifest entry: ' + name)
        listed.add(name)
    actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*')
              if p.is_file() and '__pycache__' not in p.parts and p != manifest}
    if actual != listed:
        raise ValueError('Package manifest does not cover exactly the package files')


def main(output=None):
    check_manifest()
    expected = results()[0]
    with tempfile.TemporaryDirectory(prefix='kissing25_exact_') as directory:
        source = Path(directory)
        for name in ('certificate25.json', 'verify_general.py'):
            shutil.copyfile(ROOT / 'source' / name, source / name)
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
                   MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
        subprocess.run([sys.executable, '-B', 'verify_general.py', 'certificate25.json',
                        '--output', 'fresh.json'], cwd=source, env=env, check=True,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        receipt = json.loads((source / 'fresh.json').read_text())
    if not receipt['verified'] or receipt['packing_size'] != expected['total']:
        raise ValueError('Exact proof did not establish the registered count')
    receipt['certificate'] = 'source/certificate25.json'
    if output is not None:
        output.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'status': 'PASS', 'results': [expected], 'verification': receipt}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='Optional fresh receipt outside the package')
    main(parser.parse_args().output)
