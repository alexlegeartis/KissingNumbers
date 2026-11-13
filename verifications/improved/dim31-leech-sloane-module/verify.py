"""Replay the complete geometry, exact module and pinned-source exchange."""
from pathlib import Path
import hashlib
import json
from verify_geometry import verify
from verify_module import replay
from negative_controls import main as negative_controls

ROOT = Path(__file__).resolve().parent


def check_manifest():
    expected = json.loads((ROOT / 'data/manifest.json').read_text())
    actual = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in ROOT.rglob('*') if p.is_file()
              and '__pycache__' not in p.parts and p.name != 'manifest.json'}
    if actual != expected:
        raise ValueError('Complete package manifest')


if __name__ == '__main__':
    check_manifest()
    geometry = verify(ROOT / 'data/coordinates31_251395.npz')
    expected = json.loads((ROOT / 'data/full_geometry.json').read_text())
    if {k: v for k, v in geometry.items() if k != 'seconds'} != {
            k: v for k, v in expected.items() if k != 'seconds'}:
        raise ValueError('Complete geometry receipt')
    module = replay()
    expected_module = json.loads((ROOT / 'data/module_verification.json').read_text())
    if json.loads(json.dumps({k: v for k, v in module.items() if k != 'seconds'})) != {k: v for k, v in expected_module.items() if k != 'seconds'}:
        raise ValueError('Independent integer module receipt')
    controls = negative_controls()
    expected_controls = json.loads((ROOT / 'data/negative_controls.json').read_text())
    if controls != expected_controls:
        raise ValueError('Negative-control receipt')
    print(json.dumps(dict(status='PASS', geometry=geometry, module=module,
                          negative_controls=controls), indent=2))
