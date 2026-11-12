"""Resolve and hash-check the frozen source files used in this campaign."""
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PACKINGSTAR_DEFAULT = REPO.parent / 'PackingStar'
SOURCE_AUDIT = json.loads((HERE/'reports/source_audit.json').read_text())
SOURCE_ROWS = {row['dimension']: row for row in SOURCE_AUDIT['dimensions']}
REPOSITORY_INPUTS = json.loads((HERE/'reports/repository_inputs.json').read_text())

def checked(path, expected):
    path = Path(path)
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != expected:
        raise ValueError(f'Input hash mismatch: {path}')
    return path

def packingstar_file(directory, dimension):
    row = SOURCE_ROWS[dimension]
    return checked(Path(directory)/row['path'], row['sha256'])

def repository_file(dimension):
    row = REPOSITORY_INPUTS['files'][str(dimension)]
    return checked(REPO/row['path'], row['sha256'])

def fresh_leech_shell():
    path = REPO/'verifications/improved/dim26-27-midpoint-replacements/source/verify.py'
    spec = importlib.util.spec_from_file_location('_reclaim_leech_generator', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.leech_shell().astype('int16')
