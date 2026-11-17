"""Invoke the unchanged frozen independent68/69D checker."""
from pathlib import Path
import subprocess,sys,json,numpy as np,sympy as sp,hashlib
P=Path(__file__).resolve().parent;A=P/'d4_recovery_20261117/independent_audit';S=P/'mixed_shell71_20261117'
subprocess.run(['g++','-O3','-std=c++17',str(A/'count128.cpp'),'-o',str(A/'count128')],check=True)
Q=P/'public_gamma72_provenance_20261117'
subprocess.run([sys.executable,str(Q/'verify.py')],check=True)
source=json.loads((Q/'source.json').read_text())
D=np.load(S/'leech_structure.npz');F=D['gram'];T=D['t'];G=np.array(json.loads((S/'catalogue_gram72.json').read_text()))
assert all(np.array_equal(D[k],np.array(source[v])) for k,v in [('gram','F'),('a','A'),('b','B'),('t','T')]) and np.array_equal(G,np.array(source['G']))
assert np.array_equal(2*G,T@np.kron(np.eye(3,dtype=int),F)@T.T)
L=np.load(S/'leech4.npz');K=L['K'];H=K@F@K.T
assert abs(sp.Matrix(K).det())==1 and np.array_equal(L['H'],H)
tokens=(S/'leech4.input').read_text().split();assert list(map(int,tokens[:3]))==[24,4,0]
assert float(tokens[3])==0 and int(tokens[4])==0
assert np.array_equal(np.array(tokens[5:581],dtype=np.int64).reshape(24,24),H)
assert not any(map(int,tokens[581:605]))
assert list(map(int,tokens[605:629]))==[sum((int(v)&1)<<j for j,v in enumerate(row)) for row in K]
subprocess.run([str(A/'count128'),str(S/'leech4.input'),'raw',str(S/'replayed_keys.bin'),'8'],check=True)
x=np.fromfile(S/'replayed_keys.bin',dtype=np.uint32);y=np.fromfile(S/'min4keys.bin',dtype=np.uint32)
assert len(x)==196560 and np.array_equal(np.sort(x),np.sort(y))
subprocess.run([sys.executable,str(S/'build_table.py')],check=True)
assert hashlib.sha256((S/'fiber8_table.bin').read_bytes()).hexdigest()=='356cce48ef46e88bf114740aeb1a3afeddc092d1afa326f5225044341bee2037'
for d in (68,69):subprocess.run([sys.executable,str(A/'replay_scan.py'),str(d)],check=True)
