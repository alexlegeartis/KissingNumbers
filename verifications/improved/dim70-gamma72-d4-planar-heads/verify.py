from pathlib import Path
import json,hashlib,subprocess,sys,tempfile,runpy,numpy as np
if not __debug__:raise RuntimeError('Do not use Python -O')
p=Path(__file__).resolve().parent
for path,h in json.load(open(p/'MANIFEST.json')).items():assert hashlib.sha256((p/path).read_bytes()).hexdigest()==h,path
subprocess.run([sys.executable,str(p/'check_geometry.py')],check=True)
q=p.parent/'dim68-91-lattice-code-continuation/public_gamma72_provenance_20261117'
# Rebuild all completion weights from explicit integer witnesses and the public theorem.
saved=sys.argv;sys.argv=[str(q/'verify.py')]
try:g=runpy.run_path(str(q/'verify.py'))
finally:sys.argv=saved
src=g['s'];D=np.load(p.parent/'dim68-91-lattice-code-continuation/mixed_shell71_20261117/leech_structure.npz')
assert all(np.array_equal(D[k],src[v])for k,v in [('gram','F'),('a','A'),('b','B'),('t','T')])
with tempfile.TemporaryDirectory() as tmp:
 a=Path(tmp);g['N'].astype('<u2').tofile(a/'fiber8_table.bin')
 subprocess.run(['g++','-O3','-std=c++17',str(p/'count128.cpp'),'-o',str(a/'count128')],check=True)
 subprocess.run([sys.executable,str(p/'replay.py'),str(a)],check=True)
 r=json.load(open(a/'RESULT.json'));assert r['bound']==json.load(open(p/'claims.json'))['70']
print('PASS: source provenance, all 121 populations and all same/cross-fiber geometry')
