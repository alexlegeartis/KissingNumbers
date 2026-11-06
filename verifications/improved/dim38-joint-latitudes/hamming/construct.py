from pathlib import Path
import json,numpy as np,importlib.util,hashlib
P=Path(__file__).parent;S=P.parent;old=json.loads((S/'refill/certificate.json').read_text());spec=importlib.util.spec_from_file_location('baseline_coordinates',S/'refill/verify_coordinates.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);A,_=m.reconstruct();roots=[]
for i,j in [(14,15),(16,17),(18,19),(20,21)]:
 for sign in [-1,1]:
  u=np.zeros(24,np.int64);u[i]=4;u[j]=4*sign;roots.append(u)
roots=np.array(roots);ids=[int(np.flatnonzero(np.all(A==u,axis=1))[0]) for u in roots];B=[];Y=[];metadata=[]
for half in range(2):
 for word in range(128):
  color=0
  for k in range(7):
   if word>>k&1:color^=k+1
  y=np.zeros(14,np.int64);y[7*half:7*half+7]=[1-2*((word>>k)&1) for k in range(7)]
  for sign in [-1,1]:B.append(sign*roots[color]);Y.append(y.copy());metadata.append([half,word,color,sign])
import argparse
parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');args=parser.parse_args()
if args.check:
 z=np.load(P/'new_points.npz');assert np.array_equal(z['body_integer'],B) and np.array_equal(z['head_integer'],Y) and np.array_equal(z['metadata'],metadata);c=json.loads((P/'certificate.json').read_text());assert c['body_roots']==roots.tolist() and c['body_line_indices']==ids and c['removed_cap_body_lines']==sorted(set(old['body_line_indices'])|set(ids));print('PASS: exact reconstruction of 512 Hamming points and 22 cap owners.')
else:
 np.savez_compressed(P/'new_points.npz',body_integer=np.array(B),head_integer=np.array(Y),metadata=np.array(metadata));print('Wrote new_points.npz; certificate hash must be regenerated before verification.')
