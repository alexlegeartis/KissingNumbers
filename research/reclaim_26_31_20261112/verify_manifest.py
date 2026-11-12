"""Check complete SHA-256 coverage of this frozen research campaign."""
from pathlib import Path
import hashlib
import json

root = Path(__file__).resolve().parent
expected = json.loads((root/'manifest.json').read_text())
actual = {
    p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
    for p in root.rglob('*')
    if p.is_file() and '__pycache__' not in p.parts
    and p != root/'manifest.json'
}
if actual != expected:
    missing = sorted(set(expected)-set(actual))
    extra = sorted(set(actual)-set(expected))
    changed = sorted(k for k in expected.keys() & actual.keys()
                     if expected[k] != actual[k])
    raise ValueError({'missing': missing, 'extra': extra, 'changed': changed})
print(f'PASS: all {len(actual)} research files match the complete manifest')
