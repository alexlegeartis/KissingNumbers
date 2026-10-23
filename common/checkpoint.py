"""Integrity checks and isolated replays for the October certificate packages."""
from contextlib import contextmanager
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def check_manifest(root):
    if not __debug__:
        raise SystemExit('Run without -O: certificate assertions must be enabled.')
    root = Path(root).resolve()
    manifest = root / 'data/SHA256SUMS.txt'
    listed = set()
    for line in manifest.read_text(encoding='utf-8').splitlines():
        digest, name = line.split('  ', 1)
        path = (root / name).resolve()
        assert root in path.parents and path != manifest and name not in listed
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, name
        listed.add(name)
    actual = {p.relative_to(root).as_posix() for p in root.rglob('*')
              if p.is_file() and '__pycache__' not in p.parts and p != manifest}
    assert listed == actual, 'manifest coverage'


@contextmanager
def temporary_source(root):
    with tempfile.TemporaryDirectory(prefix='kissing_checkpoint_') as directory:
        target = Path(directory) / 'source'
        shutil.copytree(Path(root) / 'source', target,
                        ignore=shutil.ignore_patterns('__pycache__'))
        yield target


def replay(source, script, *args):
    path = Path(source) / script
    env = dict(os.environ, OPENBLAS_NUM_THREADS='4', OMP_NUM_THREADS='4', MKL_NUM_THREADS='4')
    subprocess.run([sys.executable, '-B', str(path), *args], cwd=path.parent,
                   env=env, check=True)
