"""Replay both exact delta verifiers; optionally use independent brute-force audits."""
from pathlib import Path
import subprocess,sys,json,hashlib
HERE=Path(__file__).resolve().parent
if not __debug__:raise RuntimeError('Run without -O.')
for name,digest in json.loads((HERE/'source_sha256.json').read_text()).items():
 assert hashlib.sha256((HERE.parents[2]/name).read_bytes()).hexdigest()==digest, name
script='verify_independent.py' if '--independent' in sys.argv else 'verify.py'
for dim in ('dim61-62','dim66'):
 subprocess.run([sys.executable,str(HERE/dim/script)],check=True)
print('PASS: 61D 88979562; 62D 122543436; 66D 342291652')
