"""Exact 66D pentacode-Netto signed support and partial-phase replay."""
from itertools import combinations
from fractions import Fraction as Q
import json,gzip
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from prep import input66

def phase():
 seeds=['01112','03110','21310','21132'];gray=[0,1,3,2]
 P={sum(gray[(a*int(s[(j+k)%5]))%4]<<(2*j)for j in range(5))>>1 for s in seeds for k in range(5)for a in (-1,1)}
 assert len(P)==40 and min((a^b).bit_count()for a,b in combinations(P,2))==3
 return tuple(sorted(P))
def positions(b):
 return [j for j in range(66)if b>>j&1]
def pair128(A):
 return np.array([(int(x)&((1<<64)-1),int(x)>>64) for x in A],dtype=np.uint64)
def pops(x):
 fn=getattr(np,'bitwise_count',None)
 if fn is not None:return fn(x)
 x=x-((x>>np.uint64(1))&np.uint64(0x5555555555555555))
 x=(x&np.uint64(0x3333333333333333))+((x>>np.uint64(2))&np.uint64(0x3333333333333333))
 x=(x+(x>>np.uint64(4)))&np.uint64(0x0f0f0f0f0f0f0f0f)
 return (x*np.uint64(0x0101010101010101))>>np.uint64(56)
def phase_assignments(S,T):
 P=phase()
 marked=(1,2,5,6); used=set(marked);m=sum(1<<j for j in marked)
 assert len({p&m for p in P})==12
 U=pair128(S);V=pair128(T);R=[]
 for i,q in enumerate(S):
  lo,hi=U[i]
  jj=np.flatnonzero(((V[:,0]&lo)==V[:,0])&((V[:,1]&hi)==V[:,1]))
  assert len(jj)<=1
  t=T[int(jj[0])] if len(jj) else 0
  j=positions(q)
  if t:
   pr=positions(t);qt=positions(q^t)
   j=[None]*9
   for h,p in zip(marked,pr):j[h]=p
   for h,p in zip((v for v in range(9)if v not in used),qt):j[h]=p
   assert sorted(j)==positions(q)
  R.append((t,tuple(sum(1<<j[h]for h in range(9)if p>>h&1)for p in P)))
 return R
def main():
 assert __debug__
 W,T,prior=input66()
 new=list(map(int,json.loads((HERE/'supports.json').read_text())))
 S=prior+new
 phase_path=HERE/'partial.json.gz'
 extra=json.loads(gzip.decompress(phase_path.read_bytes())) if phase_path.is_file() else []
 allsup=S+[int(x[0]) for x in extra]
 assert all(s>>64 for s in new+[int(x[0]) for x in extra])
 assert len(allsup)==len(set(allsup)) and all(0<s<(1<<66)and s.bit_count()==9 for s in allsup)
 B=pair128(allsup);W=pair128(W);U0=pair128(T)
 # Base and partial supports are both freshly checked against all old weight-16 supports.
 for i in range(0,len(B),64):
  Y=B[i:i+64]
  Q0=pops(Y[:,None,0]&W[:,0])+pops(Y[:,None,1]&W[:,1])
  assert int(Q0.max())<=6
 assert all((x&y).bit_count()<=4 for x,y in combinations(S,2))
 R=phase_assignments(S,T)
 for i,(s,sg) in enumerate(extra,start=len(S)):
  s=int(s);lo,hi=B[i]
  jj=np.flatnonzero(((U0[:,0]&lo)==U0[:,0])&((U0[:,1]&hi)==U0[:,1]))
  assert len(jj)<=1
  t=T[int(jj[0])] if len(jj) else 0
  assert isinstance(sg,list) and sg and all(isinstance(p,int) and 0<=p<512 for p in sg)
  assert len(sg)==len(set(sg)) and all((p^q).bit_count()>=3 for p,q in combinations(sg,2))
  jj=positions(s)
  assert len(jj)==9
  R.append((t,tuple(sum(1<<j for h,j in enumerate(jj)if p>>h&1)for p in sg)))
 retired={};hits=0;total=0
 for s,(t,signs) in zip(allsup,R):
  assert all(p&s==p for p in signs)
  total+=len(signs)
  if not t:continue
  hits+=1
  lost={p&t for p in signs}
  retired.setdefault(t,set()).update(lost)
  assert len(retired[t])<=16
 relevant=0;comparisons=0
 for start in range(0,len(B),128):
  stop=min(len(B),start+128)
  Y=B[start:stop]
  I=pops(Y[:,None,0]&B[:stop,0])+pops(Y[:,None,1]&B[:stop,1])
  a,b=np.where((I>=5)&(np.arange(start,stop)[:,None]>np.arange(stop)[None,:]))
  for aa,bb in zip(a,b):
   i=int(aa)+start;j=int(bb)
   v=allsup[i]&allsup[j];n=v.bit_count()
   p={s&v for s in R[i][1]};q={s&v for s in R[j][1]}
   assert not (p&q)
   if n>=7:
    for j0 in positions(v):
     assert not any(s^(1<<j0) in p for s in q)
   relevant+=1;comparisons+=len(p)*len(q)
 # All raw dots are compared at squared norm 64.
 a=Q(8,3)
 assert 9*a*a==64 and (9-2*3)*a*a<=32
 assert 4*a*a<=32 and 6*2*a==32
 assert 9*a<=32 and 8*a<=32
 assert 3*4*a==32 and 4*4*a>32 and (4-2)*4*a<=32
 # Distinct new points have nine amplitudes 8/3; they are not old-layer points.
 assert all(v*a*a<=32 for v in range(-9,5))and 5*a*a>32
 old=342297172;gone=sum(len(c)for c in retired.values())
 assert len(prior)*40==5520
 gain=total-5520-gone
 assert gain>0
 return dict(dimension=66,old=old,total=old+gain,gain=gain,selected_supports=len(new),
             partial_supports=len(extra),one_netto=hits,removed_netto_supports=len(retired),
             removed_sign_points=gone,added_sign_points=total-5520,
             stock_weight16=len(W),netto_four_supports=len(T),overlap_pair_count=relevant,
             projected_sign_pair_count=comparisons)
if __name__=='__main__':print(json.dumps(main(),indent=2))

