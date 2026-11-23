from pathlib import Path
from hashlib import sha1,sha256
from fractions import Fraction as Q
import base64,zlib,struct,json,re
import numpy as np
assert __debug__
p=Path(__file__).resolve().parent
def pin(rel,h):
 b=(p.parent/rel).read_bytes()
 assert sha1(f'blob {len(b)}\0'.encode()+b).hexdigest()==h,rel
 return b
r=json.loads(pin('dim86-87-90-rotated-plane-trades/claims.json','f1e8bcb12be6b4dbcceea4a663e752a878cb0e91'))
a=json.loads(pin('dim86-87-remarked-square/claims.json','aaa0c8ee25ddcb542644f70dd3f9b0e024981c57'))
b=json.loads(pin('dim88-full-pless-chord-marks/claims.json','8dd77d3f09a22d2edc205653358d299ef2cbf015'))
s=pin('dim90-92-hexagonal-phase-layers/hexagon_sparse15_phase7812_20261117/PROOF.md','c65378853f0ec3ebe420ebabd2c7894fb6551662')
q=int(re.search(rb'\|45\|0,1,2\|(\d+)\|',s)[1])
base={d['dimension']:d['lower_bound']for d in r['dimensions']};base[88]=b['lower_bound']
occ={d['dimension']:d['histogram'][0]for d in a};occ.update({88:b['histogram'][0],90:q})
raw=zlib.decompress(base64.b64decode((p/'phase.b64').read_bytes()))
v=np.array([n for n, in struct.iter_unpack('<I',raw)],dtype='i8')
assert len(v)==len(set(v))==7906 and max(v)<3**15
digits=v[:,None]//(3**np.arange(14,-1,-1,dtype='i8'))%3
words=(b'\n'.join(bytes((r+48).astype('u1'))for r in digits)+b'\n')
assert sha256(words).hexdigest()=='220a8d800eaaf391898b424e5e27f54fe957c4f0548d6ce980e9b488637ab369'
A=(digits==1)@(1<<np.arange(15,dtype='i4'))
B=(digits==2)@(1<<np.arange(15,dtype='i4'))
P=np.array([i.bit_count()for i in range(1<<15)],dtype='u1')
mini=min(int(P[(A[i]^A[i+1:])|(B[i]^B[i+1:])].min())for i in range(7905))
assert mini==5
assert 3*15==9*mini==45
assert Q(45,2)**2==Q(45*45,4)  # mixed against dense, both squared norm 45
assert Q(7,4)**2>3 and 24+12*Q(7,4)<=45
assert (Q(15,2))**2==Q(45*5,4) and 27<Q(45*3,4) and 12<Q(45*2,4)
out={'status':'PASS','lower_bounds':{str(n):base[n]+occ[n]*(len(v)-7812)for n in sorted(base)}}
assert out==json.loads((p/'claims.json').read_text())
print(json.dumps(out))
