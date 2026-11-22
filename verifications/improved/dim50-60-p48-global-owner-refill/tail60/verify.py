#!/usr/bin/env python3
"""Exact P48 pole splitting, no optimization."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from math import lcm
import argparse,json,gzip,numpy as np

HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--siblings',type=Path,default=HERE.parent)
args=parser.parse_args()
P=args.siblings/'dim49-60-p48-pointwise-refill'
D=args.siblings/'dim49-60-p48-signed-refill'

def pin(path, digest):assert sha256(path.read_bytes()).hexdigest()==digest,path;return path
head_file=pin(P/'heads60.npz','cc3a2e07fbe3f5f91c33d61ba69102057c5fd8f90dc509d3c5d83d29adac9b25')
gram_file=pin(P/'source/p48p_gram.npy','01a5d107e4424092e5300e201592aaff62804ee27eeb1d23b6e516249897d3e1')
delta_file=pin(D/'p48_signed_delta60.json.gz','6e516dedf722e13c6c660a873285332d73ee0f813ad8cfd39e2d97c5715612ac')
G=np.load(gram_file,allow_pickle=False).astype('i8')
U=np.load(head_file,allow_pickle=False)
config=json.loads(gzip.decompress(delta_file.read_bytes()))
data=json.loads((HERE/'witness.json').read_text())
W=data['frame']
sparse=np.asarray(data['sparse'],dtype='i8')
up={r['color']:[None if q is None else np.asarray(q,dtype='i8')for q in r['tetrahedra']]for r in data['upgrade']}
upall=[q for r in up.values()for q in r if q is not None]
assert len(up)==len(data['upgrade'])
assert config['dimension']==60 and config['bound']==58631571
assert G.shape==(48,48) and np.array_equal(U['gram'],G)
Hs=np.array(W,dtype='i8')
V=U['base'].astype('i8');Q=U['extra_points'].astype('i8');T=U['transforms'].astype('i8')
O=U['extra_offsets'];bits=np.unpackbits(U['retained_bits'],axis=1)[:,:len(V)].astype(bool)
C={x['class_index']:x['add']for x in config['changes']}
omega=np.random.default_rng(41).integers(-10**7,10**7,48,dtype='i8')
def m(a):return max(map(abs,map(int,a.flat)),default=0)
assert all(v.dtype.kind in 'iu' for v in (G,V,Q,T))
def ip(a,b):
 assert a.shape[-1]==b.shape[0] and a.shape[-1]*m(a)*m(b)<2**62
 return a@b
assert np.array_equal(ip(ip(Hs,G),Hs.T),6*np.eye(len(W),dtype='i8'))
assert sparse.shape==(4,4,48) and all(len(v)==8 for v in up.values()) and all(v.shape==(4,48)for v in upall)
allheads=np.concatenate([Hs,sparse.reshape(-1,48),*upall],axis=0)
def projective(a):
    a=tuple(map(int,a));return min(a,tuple(-x for x in a))
assert len(set(map(projective,allheads)))==len(allheads)
assert np.all(np.sum(ip(allheads,G)*allheads,axis=1)==6)
target=set(np.abs(ip(allheads,omega)).tolist())
for j,Uj in enumerate(T):
 Cj=ip(Uj,omega)
 for X in (V[bits[j]],Q[O[j]:O[j+1]],np.array(C.get(j,[]),dtype='i8').reshape(-1,48)):
  assert target.isdisjoint(np.abs(ip(X,Cj)).tolist()), ('ownership',j)

from itertools import product
from collections import Counter
from functools import reduce
from operator import xor
assert Hs.shape==(len(W),48) and 1<=len(W)<=16
supports=[]
for j in (0,6):
 for k in (0,2,4):
  supports.append(tuple(j+i for i in range(6) if i not in (k,k+1)))
for c in range(5):
 M=[(c,5),((c+1)%5,(c-1)%5),((c+2)%5,(c-2)%5)]
 for e,f in product(M,repeat=2): supports.append(tuple(e)+tuple(j+6 for j in f))
D=[]
for i in range(12):
 for a in (-2,2):
  v=[0]*12;v[i]=a;D.append(v)
assert len(supports)==51 and len(set(tuple(sorted(e)) for e in supports))==51
for e in supports:
 for signs in product((-1,1),repeat=4):
  v=[0]*12
  for i,a in zip(e,signs,strict=True):v[i]=a
  D.append(v)
D=np.array(D,dtype='i8'); M=D@D.T
assert len(D)==840 and np.all(np.diag(M)==4)
assert np.all(M[np.triu_indices(840,1)]<=2)
assert all(tuple(-v) in set(map(tuple,D)) for v in D)
assert np.max(np.abs(D).sum(1))<=4
balls={}
for w in range(1<<12):
 if w.bit_count()%2==0:
  u=reduce(xor,(j for j in range(12) if (w>>j)&1),0)
  balls.setdefault(u,[]).append(w)
assert set(balls)==set(range(16))
for words in balls.values():
 assert len(words)==128 and all((a^b).bit_count()>=4 for a,b in combinations(words,2))
# All new points are (±u/2, sqrt(3) z/2), with u a unit unowned head
# and z a unit sign tail. Old poles have |z·pole|<=1/sqrt(3) by L1.
# Old lifts score <=sqrt(3)/8+1/4 < 1/2 by ownership.
# Same new color scores <=1/4+(3/4)(1/3); different colors:
# orthogonal heads and distance >=2 imply score <=(3/4)(2/3).
assert F(1,4)+F(3,4)*F(1,3)==F(3,4)*F(2,3)==F(1,2)
# Four palettes indexed by the complement-pair class of a triple (Q3 modulo
# the all-one word). A triple chooses one vertex of each marked matching pair.
def headgram(A,B):return ip(ip(A,G),B.T)
def tetra(X):assert np.array_equal(headgram(X,X),8*np.eye(4,dtype='i8')-2)
for j,X in enumerate(sparse):
    tetra(X)
    for Y in sparse[:j]:
        assert np.max(abs(headgram(X,Y)))<=1
triples={}
for b in (0,6):
    for w in range(8):
        T=frozenset(b+2*j+((w>>j)&1) for j in range(3))
        c=min(w,w^7)
        triples[T]=c
assert len(triples)==16 and set(triples.values())==set(range(4))
for t,c in triples.items():
    for v,e in triples.items():
        if t==v:continue
        assert len(t&v)<=(2 if c!=e else 0)
    assert all(len(t&set(e))<=2 for e in supports)
all_triples=set(frozenset(t) for t in combinations(range(12),3))
assert set(triples)==all_triples-{frozenset(t)for e in supports for t in combinations(e,3)}
# Old-pole cap |z·p| <= 1/sqrt(3). For sparse-lift cross pairs use
# root means sqrt(2)<99/70, sqrt(30)<137/25.
assert F(5,24)<F(1,4) and F(99,70)**2>2 and F(137,25)**2>30
assert 3*F(99,70)/16+F(137,25)/24<F(1,2)
# Sparse-sparse, same tail: tetrahedron. Same support and equal tail-sign
# parity: Hamming at least 2. Different parity: the two signed head banks
# are opposite tetrahedra, whose max cross normalized product is 1/3.
assert -F(3,8)/3+F(5,8)==F(1,2)
assert F(3,8)-F(5,8)/3<F(1,2)
assert F(3,8)/3+F(5,8)/3<F(1,2)
assert F(3,8)/6+F(5,8)*F(2,3)<F(1,2) and F(3,8)<F(1,2)
# Sparse-dense: any 3-supported tail has absolute product <=1/2 with
# a unit uniform cube word. Any distinct signed P48 heads have product <=1/2.
# For arbitrary two altitudes Cauchy gives a combined cap 1/2.
# Independent holdout exact all tail words: subcosets used for the heavy
# cube layer and the modified 16-class construction.
L=set(data['subcode'])
assert len(L)==16 and 0 in L and L<=set(balls[0])
assert all((a^b) in L for a in L for b in L)
assert all(x==0 or x.bit_count()>=6 for x in L)
assert set(up)<=set(balls) and len(W)==16
assert all(set(balls[i])==set(w^balls[i][0] for w in balls[0]) for i in balls)
part={}
for i in up:
    remaining=set(balls[i]); palette=[]
    while remaining:
        z=min(remaining)
        coset={z^w for w in L}
        assert coset<=remaining
        palette.append(coset)
        remaining-=coset
    part[i]=palette
for i,A in up.items():
    for X in A:
        if X is None:continue
        tetra(X)
        low=[j for j in balls if j!=i and (j not in up or any(u is None for u in up[j]))]
        assert all(np.max(abs(headgram(X,Hs[j:j+1])))<=1 for j in low)
    palette=part[i]
    assert len(palette)==8
    for C in palette:
        assert all((x^y).bit_count()>=6 for x,y in combinations(C,2))
    for j,B in up.items():
        if i<=j:continue
        for a,C in enumerate(part[i]):
            for b,D in enumerate(part[j]):
                if A[a] is not None and B[b] is not None and any((u^v).bit_count()==2 for u in C for v in D):
                    assert np.max(abs(headgram(A[a],B[b])))<=1
# Larger same-parity distances d>=4 across two heavy colors are safe for
# distinct heads by the existing head-shell bound 1/2.
# Dense-heavy pure pole, old lifted stock (rational squared upper comparisons).
assert F(3,5)/3<F(1,4)
assert F(2191,2000)**2>F(6,5) and F(223607,100000)**2>5
assert F(2191,2000)/4+F(223607,100000)/10<F(1,2)
# Mixed dense hemispheres from different parity-equal syndromes use
# max tail product 2/3, max head product 1/6. sqrt(10)/60 + sqrt(5)/5 < 1/2.
assert F(316228,100000)**2>10
assert F(316228,100000)/60+F(223607,100000)/5<F(1,2)
# Dense-heavy pairs in same word / same subcoset / original syndrome / different syndromes.
assert -F(2,5)/3+F(3,5)<F(1,2) # same word, different tetra directions
assert F(2,5)<F(1,2) # same subcoset, different tails at Hamming >=6
assert F(2,5)/2+F(3,5)/3<F(1,2) # different subcosets, original syndrome equal
assert F(2,5)/6+F(3,5)*F(2,3)<F(1,2) # different original syndromes
assert len(set(map(projective,sparse.reshape(-1,48))))==16
removed_eq=32+4*len(upall)
add=4096+32*len(upall)+16*8*4
new=config['bound']-1-removed_eq+add
assert new==58636146+28*len(upall)
print(json.dumps(dict(dimension=60,previous=58635154,frame_colors=16,heavy_syndromes=len(up),heavy_subcosets=len(upall),
    retired_poles_net=1,equator_retirements=removed_eq,new_rays=add,lower_bound=new)))
