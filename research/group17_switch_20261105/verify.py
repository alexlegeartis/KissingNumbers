"""Exact verification and inequivalence of two 5730-point 17D codes."""
from pathlib import Path
from collections import Counter
import numpy as np
z=np.load(Path(__file__).with_name('coordinates.npz'));X=z['source'].astype(np.int64);Y=z['switched'].astype(np.int64);Z=z['switched_twice'].astype(np.int64);g=z['weights'].astype(np.int64);norm=int(z['norm']);assert X.shape==Y.shape==Z.shape==(5730,17) and tuple(g)==(1,)*16+(2,) and norm==93312
profiles=[]
for name,P in [('source',X),('switched',Y),('switched_twice',Z)]:
 assert len(np.unique(P,axis=0))==len(P) and np.all(((P*P)*g).sum(1)==norm)
 h=Counter()
 for st in range(0,len(P),96):
  G=(P[st:st+96]*g)@P.T;ix=np.arange(st,min(st+96,len(P)));off=G.copy();off[np.arange(len(ix)),ix]=0;assert np.max(off)<=norm//2
  for row in G:
   value,count=np.unique(row,return_counts=True);h[tuple(zip(map(int,value),map(int,count)))]+=1
 profiles.append(h);print(name,'count',len(P),'per-point profile classes',len(h),flush=True)
sx=set(map(tuple,X));sy=set(map(tuple,Y));assert len(sx&sy)==5602 and len(sx-sy)==128 and len(sy-sx)==128
def overall(h):
 result=Counter()
 for profile,multiplicity in h.items():
  for dot,count in profile: result[dot]+=multiplicity*count
 return result
assert overall(profiles[0])==overall(profiles[1])==overall(profiles[2])
for i in range(3):
 for j in range(i):
  assert profiles[i]!=profiles[j];w=min(k for k in profiles[i].keys()|profiles[j].keys() if profiles[i][k]!=profiles[j][k]);print('inequivalence',i,j,'witness multiplicities',profiles[i][w],profiles[j][w],flush=True)
sz=set(map(tuple,Z));assert len(sx&sz)==5474 and len(sy&sz)==5602
# Check the second switch applied without the first, so independence is explicit.
D=np.array(sorted((sx-(sy-sz))|(sz-sy)),dtype=np.int64);assert D.shape==(5730,17)
for st in range(0,len(D),96):
 G=(D[st:st+96]*g)@D.T;ix=np.arange(st,min(st+96,len(D)));G[np.arange(len(ix)),ix]=0;assert np.max(G)<=norm//2
print('both switches are independently applicable',flush=True)
# Each chosen switched point is among explicitly stated rational group images.
u=np.zeros(17,np.int64);v=u.copy();w=u.copy();u[[0,16]]=1;v[[1,16]]=1;w[[2,3,4,5]]=1;gen=[u,v,w];image=set()
for word in [[],[0],[1],[2],[0,1],[0,2],[1,2],[0,1,2]]:
 P=X.copy()
 for j in word:
  t=gen[j];den=int((t*g)@t);numer=2*np.outer(P@(t*g),t);assert not np.any(numer%den);P=P-numer//den
 image.update(map(tuple,P))
assert sy<=image and sz<=image
print('two independent 128-for-128 switches; all chosen points belong to the eight stated images; no new kissing bound',flush=True)
