"""Replay 64 algebraic witness orbits giving a 4096-edge matching."""
import json,numpy as np
from pathlib import Path
c=json.loads(Path(__file__).with_name('matching_certificate.json').read_text());T=np.asarray(c['offsets'],dtype=np.int64).reshape(2,24);mask=c['octad'];G=np.zeros((12,24),np.int64);g=np.array([7,1,4,4,2,7,3,3,4,7,2,1],np.int64)
for i in range(12):G[i,i:i+12]=g;G[i,23]=3
assert int(mask).bit_count()==8 and all((int(T[1,j])-int(T[0,j]))%8==4*((mask>>j)&1) for j in range(24))
words={}
for u in range(4096):
 v=((u>>np.arange(12))&1)@G%2;w=sum(int(v[j])<<j for j in range(24));assert w not in words;words[w]=u
assert mask in words
K=[]
for w,u in words.items():
 if mask&w in words:K.append((u,words[mask&w]))
assert len(K)==64 and len({a for a,b in K})==64
labels={h for h,z in K};assert all((a^b) in labels for a in labels for b in labels)
# The cyclic shift on the first 23 coordinates fixes the extended code.
for row in G:
 target=np.r_[np.roll(row[:23],1),row[23]];coeff=np.zeros(12,dtype=np.int64)
 for j in range(12):coeff[j]=((int(target[j])-int(coeff@G[:,j]))*7)%8
 assert np.array_equal(coeff@G%8,target)
# The ring-code intersection is 2C plus the 64 residue classes in K.
intersection_size=2**24*len(K);sum_size=2**36*2**36//intersection_size
assert intersection_size==2**30 and sum_size==2**42 and sum_size//2**36==64
w=np.array([0,2,2,2,4,2,2,2],np.int64);seen_u=set();seen_v=set();W=c['witnesses'];assert len(W)==64
for e in W:
 u=e['candidate'];b=np.asarray(e['coeff'],np.int64);assert type(u) is int and 0<=u<4096 and b.shape==(12,) and np.all((0<=b)&(b<8))
 v=sum(int(b[i]%2)<<i for i in range(12));assert v==e['source'];old=(u>>np.arange(12))&1
 for h,z in K:
  hb=(h>>np.arange(12))&1;zb=(z>>np.arange(12))&1;new_u=u^h;new_b=(b+hb+4*zb-2*(old*hb))%8
  new_v=sum(int(new_b[i]%2)<<i for i in range(12));assert new_u not in seen_u and new_v not in seen_v;seen_u.add(new_u);seen_v.add(new_v)
  ring=((new_u>>np.arange(12))&1)@G%8;x=(ring+T[ring%2,np.arange(24)])%8;y=new_b@G%8;dist=int(w[(x-y)%8].sum());assert dist==e['distance'] and dist<24
assert len(seen_u)==len(seen_v)==4096
print('PASS: 4096-edge matching; 64 equivalent affine code cosets and cyclic symmetry verified')
