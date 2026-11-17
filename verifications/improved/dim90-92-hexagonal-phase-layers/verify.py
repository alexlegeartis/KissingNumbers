"""Invoke the original independent checks and verify the published90–92 totals."""
from pathlib import Path
import subprocess,sys,json,numpy as np
P=Path(__file__).resolve().parent
for path in ['hexagon_sparse15_independent_20261117/verify.py','hexagon_sparse15_phase7812_20261117/verify.py','hexagon_sparse15_independent_20261117/verify_phase7812.py']:
 subprocess.run([sys.executable,str(P/path)],cwd=P,check=True)
G=np.array(json.loads((P/'hexagon_sparse15_phase7812_20261117/outer_certificate.json').read_text())['outer_generator'])
S=np.zeros((24,24),int);S[0,1:]=1;S[1:,0]=2
for i in range(1,24):
 for j in range(1,24):
  if i!=j:S[i,j]=1 if pow((i-j)%23,11,23)==1 else 2
assert np.array_equal(G,np.concatenate([np.eye(24,dtype=int),S],axis=1))
Q=P.parent/'dim68-91-lattice-code-continuation'
subprocess.run([sys.executable,str(Q/'ternary_qr48_independent_20261117/replay_generator.py')],check=True)
c=json.loads((Q/'ternary_qr48_20261117/certificate.json').read_text());assert c['counts']['90']==2*3**21
H=np.array([[2,0],[1,1],[-1,1],[-2,0],[-1,-1],[1,-1]]);Z=H@np.diag([1,3])@H.T
assert np.all(np.diag(Z)==4)and np.all(Z[np.triu_indices(6,1)]<=2)
base=2*3**21+130944*7812
assert json.loads((P/'claims.json').read_text())=={'90':base,'91':base+2,'92':base+6}
print('PASS: approved90–92 hexagonal phase-layer bounds',base,base+2,base+6)
