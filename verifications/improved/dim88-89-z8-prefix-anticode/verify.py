from pathlib import Path
import json,subprocess,sys,hashlib
import parent_independent
P=Path(__file__).resolve().parent
raw=json.loads((P/'z8_generator.json').read_text())
_,report=parent_independent.parent_proof(raw)
assert report==json.loads((P/'parent_receipt.json').read_text())['z8']
subprocess.run([sys.executable,str(P/'check.py')],check=True)
subprocess.run([sys.executable,str(P/'mixed89/check.py')],check=True)
assert json.loads((P/'claims.json').read_text())=={'88':12884901888,'89':13958643714}
print('PASS complete 88D parent/fiber/geometry proof; 89D mixed-latitude proof')
