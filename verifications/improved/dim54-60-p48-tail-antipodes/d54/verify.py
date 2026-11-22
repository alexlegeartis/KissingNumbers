#!/usr/bin/env python3
"""Exact 54D minuscule-tail increment over the signed P48 witness."""
from pathlib import Path
from itertools import combinations, product
from hashlib import sha256
import argparse, gzip, json, numpy as np
HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--siblings',type=Path,default=HERE.parent)
s=parser.parse_args().siblings
P=s/'dim49-60-p48-pointwise-refill';C=s/'dim49-60-p48-signed-refill'
def checked(p,h):
    assert sha256(p.read_bytes()).hexdigest()==h,p
    return p
X=np.load(checked(P/'heads54.npz','4fe02e7dd69ce14f81833ac2d15f2ec8665125ad1c06bd4c15a994d579f8e112'),allow_pickle=False)
G=np.load(checked(P/'source/p48p_gram.npy','01a5d107e4424092e5300e201592aaff62804ee27eeb1d23b6e516249897d3e1'),allow_pickle=False).astype('i8')
assert np.array_equal(G,X['gram'])
cfg=json.loads(gzip.decompress(checked(C/'p48_signed_delta54.json.gz','9622fce8977fe03d3f4fdd565f7173ea7696b3e38d02c905a5e25ae6a7080ce6').read_bytes()))
assert cfg['dimension']==54 and cfg['bound']==52976623
H=np.array(json.loads((HERE/'witness.json').read_text())['head_triple'],dtype='i8')
assert H.shape==(3,48) and not H.sum(axis=0).any()
def mm(a,b):
    ma=max(map(abs,map(int,a.flat)),default=0)
    mb=max(map(abs,map(int,b.flat)),default=0)
    assert a.shape[-1]==b.shape[0] and 48*ma*mb<2**61
    return a@b
assert np.array_equal(mm(mm(H,G),H.T),9*np.eye(3,dtype='i8')-3*np.ones((3,3),dtype='i8'))
V=X['base'].astype('i8'); Q=X['extra_points'].astype('i8'); T=X['transforms'].astype('i8')
O=X['extra_offsets']; keep=np.unpackbits(X['retained_bits'],axis=1)[:,:len(V)].astype(bool)
adds={int(e['class_index']):e['add']for e in cfg['changes']}
omega=np.random.default_rng(109).integers(-10**7,10**7,48,dtype='i8')
targ=set(map(int,abs(mm(H,omega))))
for j,Tj in enumerate(T):
    p=mm(Tj,omega)
    for rows in (V[keep[j]],Q[O[j]:O[j+1]],np.array(adds.get(j,[]),dtype='i8').reshape(-1,48)):
        assert targ.isdisjoint(abs(mm(rows,p)).tolist()),('old ownership',j)
# Reconstruct the SAME E8, E7, E6 integer roots as the pinned predecessor.
R=[]
for i,j in combinations(range(8),2):
    for a,b in product((-2,2),repeat=2):
        v=[0]*8;v[i],v[j]=a,b;R.append(v)
R.extend(v for v in product((-1,1),repeat=8) if v.count(-1)%2==0)
R=np.array(R,dtype='i8');E7=R[R[:,6]+R[:,7]==0];E6=E7[E7[:,5]+E7[:,7]==0]
assert len(E7)==126 and len(E6)==72 and np.all((E6*E6).sum(1)==8)
w=np.array([0,0,0,0,0,2,-1,1],dtype='i8')
A=3*E7[E7@w==4]-2*w
assert A.shape==(27,8) and len(np.unique(A,axis=0))==27 and not any(A[:,6]+A[:,7]) and not any(A[:,5]+A[:,7])
assert np.max(abs(A@E6.T))==12 and np.all((A*A).sum(1)==48)
AA=A@A.T;i,j=np.triu_indices(27,1)
assert max(AA[i,j])==12 and min(AA[i,j])==-24
# Head triples normalized by sqrt6. Tails normalized by sqrt48.
# New ray: (u/sqrt3,sqrt(2/3)z); use -H with -A.
from fractions import Fraction as F
assert F(2,3)*F(12**2,48*8)==F(1,4)
assert F(1,3)+F(2,3)*F(12,48)==F(1,2)
assert -F(1,3)/2+F(2,3)==F(1,2)
assert F(1,3)/2+F(2,3)*F(24,48)==F(1,2)
# To all old equator heads except removed identical one: cap 1/(2sqrt3).
# Old lifts: the head term <=1/4, tail term <=1/4 by the root bound above.
# New and retained lifts have distinct head/tail norms; different new tuples distinct.
print(json.dumps({'dimension':54,'previous':cfg['bound'],'add':2*27*3,'remove':6,'lower_bound':cfg['bound']+2*27*3-6}))
