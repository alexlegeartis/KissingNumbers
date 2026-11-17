from pathlib import Path
import subprocess,sys,json
P=Path(__file__).resolve().parent
subprocess.run([sys.executable,str(P.parent/'dim93-96-three-crosspolytope-alphabets/verify.py')],check=True)
for f in ['z8_fourrow92_partitions_20261117/verify.py','z8_fourrow92_independent_20261117/audit.py']:
 subprocess.run([sys.executable,str(P/'research'/f)],check=True,cwd=P)
assert json.loads((P/'claims.json').read_text())=={'92':29797911144}
print('PASS: independently audited four-row 92D partition-head construction')
