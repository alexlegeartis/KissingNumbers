"""Exact translation-kernel and antipodality classification of this witness."""
from pathlib import Path
import collections
import hashlib
import json

ROOT = Path(__file__).resolve().parent
if not __debug__:
    raise SystemExit('Run without -O: assertions are required.')
certificate = ROOT / 'certificate32_352992_fixed_body_explicit.json'
raw = certificate.read_bytes()
assert hashlib.sha256(raw).hexdigest() == 'ab06963895d074d8dd6162bea2298796391dc876330ae997a5324adfb6638c5e'
document = json.loads(raw)
fibers = collections.defaultdict(set)
for support, negative in document['exceptional_signed_points']:
    fibers[support].add(negative)
rows = []
for support, signs in sorted(fibers.items()):
    bits = [1 << j for j in range(32) if support & (1 << j)]
    translations = [sum(bits[j] for j in range(8) if label & (1 << j)) for label in range(256)]
    kernel = [t for t in translations if {n ^ t for n in signs} == signs]
    rows.append(dict(support=support, size=len(signs), missing_antipodes=sum((n ^ support) not in signs for n in signs), translation_kernel_size=len(kernel), affine=len(kernel) == len(signs)))
observed = dict(certificate_file=certificate.name, certificate_sha256=hashlib.sha256(raw).hexdigest(), points=sum(r['size'] for r in rows), supports=len(rows), missing_antipodes=sum(r['missing_antipodes'] for r in rows), affine_fibers=sum(r['affine'] for r in rows), nonlinear_fibers=sum(not r['affine'] for r in rows), partial_antipodality_fibers=sum(0 < r['missing_antipodes'] < r['size'] for r in rows), fibers=rows)
assert observed == json.loads((ROOT / 'fiber_classification.json').read_text())
assert (observed['affine_fibers'], observed['nonlinear_fibers'], observed['partial_antipodality_fibers']) == (14, 15, 15)
assert all(0 < r['missing_antipodes'] < r['size'] for r in rows if not r['affine'])
print(json.dumps(dict(status='PASS', affine_fibers=14, nonlinear_fibers=15, partial_antipodality_fibers=15)))
