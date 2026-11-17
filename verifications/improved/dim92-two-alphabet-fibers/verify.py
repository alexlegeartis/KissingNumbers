from pathlib import Path
import subprocess,sys,json
P=Path(__file__).resolve().parent
subprocess.run([sys.executable,str(P.parent/'dim93-96-three-crosspolytope-alphabets/verify.py')],check=True)
for f in ['ring2_final92_95_20261117/verify_dense.py','ring2_final92_95_20261117/verify_sparse.py','ring2_final92_95_20261117/verify_rs17.py','ring2_92_independent_20261117/audit.py']:
 subprocess.run([sys.executable,str(P/'research'/f)],check=True,cwd=P)
assert json.loads((P/'claims.json').read_text())=={'92':26038868280}
print('PASS: approved 92D two-alphabet fiber construction')
