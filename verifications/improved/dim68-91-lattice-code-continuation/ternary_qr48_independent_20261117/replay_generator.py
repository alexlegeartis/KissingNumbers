"""Exact finite checks, not an independent minimum-distance computation.
The [48,24,15] minimum distance is the established QR48 theorem; see proof.md.
"""
import json, hashlib
from pathlib import Path
import sympy as sp
P=Path(__file__).resolve().parent
x=sp.symbols('x')
g=x**23+x**20-x**18+x**17-x**15+x**12+x**10-x**9-x**7-x**6-x**5-x**2-x-1
h=x**23+x**22+x**21+x**18+x**17+x**16+x**14-x**13-x**11+x**8-x**6+x**5-x**3-1
assert sp.Poly((x-1)*g*h-(x**47-1),x, modulus=3).is_zero
f=sp.Poly((x-1)*g,x)
G=[[int(f.nth(j-i))%3 if j>=i else 0 for j in range(47)]+[0] for i in range(23)]+[[1]*48]
def rref(a):
 a=[r[:] for r in a]; piv=[]; k=0
 for j in range(len(a[0])):
  v=next((i for i in range(k,len(a)) if a[i][j]%3),None)
  if v is None: continue
  a[k],a[v]=a[v],a[k]; inv=pow(a[k][j],-1,3);a[k]=[z*inv%3 for z in a[k]]
  for i in range(len(a)):
   if i!=k:
    b=a[i][j];a[i]=[(v-b*w)%3 for v,w in zip(a[i],a[k])]
  piv.append(j);k+=1
  if k==len(a):break
 return a[:k],piv
def shorten(G,s):
 # Gaussian elimination puts pivots in first s columns; discard those rows and columns.
 R,piv=rref(G);assert piv[:s]==list(range(s))
 return [r[s:] for r in R[s:]]
assert len(rref(G)[0])==24
assert all(sum(a*b for a,b in zip(u,v))%3==0 for u in G for v in G)
Gs={str(s):shorten(G,s) for s in (3,4)}
for s in (3,4): assert len(rref(Gs[str(s)])[0])==24-s
assert len(rref([r[:2] for r in Gs['3']])[0])==2
# Gleason enumerator, assuming established minimum >=15.
t=sp.symbols('t'); F=1+8*t; H=t*(1-t)**3
W=sp.Poly(F**12-96*F**9*H+2400*F**6*H**2-12704*F**3*H**3+96*H**4,t)
weights={3*i:int(W.nth(i)) for i in range(17)}
assert sum(weights.values())==3**24
assert all(weights[i]==0 for i in (3,6,9,12))
data={'generator_Q48':G,'shortened_generators':Gs,'weight_enumerator_assuming_Q48_theorem':weights,'counts':{'88':2*3**20,'89':10*3**19,'90':2*3**21,'91':2*3**21}}
out=P/'certificate.json';out.write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps({'rank':24,'self_dual':True,'shortened_ranks':{'45':21,'44':20},'counts':data['counts'],'A45':weights[45],'sha256_certificate':hashlib.sha256(out.read_bytes()).hexdigest()},indent=2))
