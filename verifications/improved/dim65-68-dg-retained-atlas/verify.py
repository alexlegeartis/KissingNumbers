#!/usr/bin/env python3
"""Exact 65–68D affine-atlas and split DG composition."""
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
from subprocess import run, check_output
from dg import code
import gzip
import json
import sys

assert __debug__
P=Path(__file__).resolve().parent
R=P.parent
signed=R/'dim63-67-signed-overlap'
original=R/'dim63-68-affine-hyperplane-heads'
atlas=R/'dim65-68-complete-hyperplane-atlas'
opt=R/'dim65-optimal-quadruple-layer'
grouped=R/'dim65-66-grouped-four-site-designs'
nested=R/'dim66-68-nested-affine-face-trades'

def pin(p,h,gz=False):
    b=p.read_bytes()
    if gz:b=gzip.decompress(b)
    assert sha256(b).hexdigest()==h,p
    return b
sources={
 'cubic80.txt': pin(signed/'cubic80.txt','5173adcf040a2c2978522f781aef88cab6cbfad7c7bcde6f1c99584554f32c8f'),
 'signed_witness.txt': pin(signed/'signed_witness.txt.gz','5173fabd124bd23e077d00c956534e7f56543101b28903cbac017ed252c34be7',True),
 'tagged_witness.txt': pin(original/'tagged_witness.txt.gz','cb51472a7bcbe09f85728b54f74ff1780afaf3a817f4bc442dd7d8b65ba26fb2',True),
 'hyperplane_heads.txt': pin(atlas/'hyperplane_heads.txt','e615223cf61ab75b2a923cc3ecc68b314f6a826c10428b67b5ff10b9da7e21fb'),
 'joint_witness.txt': pin(atlas/'joint_witness.txt','311ebf0e8de26703fc4f3b790a0f8d81e2f5c183eb154dd2131d530714b409ec'),
 'restored68.txt': pin(P/'restored68.txt','66e410c145caafbca1fa019dd4fc88ee60be320310f82c894f29c6636857ae6d'),
 'assignment.txt': pin(P/'assignment.txt','12b23c1814c6b4879195c12619b899a77d4a2f614f4d462c79ec65be06bbade7'),
}
previous={
 65:json.loads(pin(opt/'claims.json','a6ab64b9f793c54568be50a995ced72e58754698d41e2375728c7c8b9733c160'))['lower_bounds']['65'],
 66:json.loads(pin(grouped/'claims.json','2b7246370e971ac89dd0f4ba1303f39b308efa1aabe85aa61f7f6c5534365ee8'))['lower_bounds']['66'],
}
old=json.loads(pin(nested/'claims.json','ae0f187af9eb9c32d8841d5106c50f60a640a8ade1f0142993d3a96d71a8a30f'))['lower_bounds']
previous.update({n:old[str(n)] for n in (67,68)})
ph=json.loads(pin(P/'phase68.json','580fbcc1ebe442ea2283de43e84d67651694d7cd389b8bb21a10b0bb1d07eebb'))
assert (ph['N'],ph['E'],ph['rejected'])==(67923,31512,[])
X=int(ph['edgebits'],16);assert X>>ph['E']==0
sources['phase.txt']=(''.join(str(X>>i&1) for i in range(ph['E']))+'\n').encode()
with TemporaryDirectory() as t:
    T=Path(t)
    for name,b in sources.items():(T/name).write_bytes(b)
    z=T/'verify'
    run(['g++','-O3','-march=native','-mbmi2','-std=c++20',str(P/'verify_atlas.cpp'),'-o',str(z)],check=True)
    s=check_output([str(z)],cwd=T,text=True)
    print(s,end='')
rows=[tuple(map(int,l.split()[1:])) for l in s.splitlines() if l.startswith('ATLAS ')]
assert [n for n,_,_ in rows]==[65,66,67,68]
C,A=code()
assert len(C)==1024 and len(A)==128
def mon(u):return sum(1<<i for i in range(64)if i&u==u)
def span(iter):
    s={0}
    for x in iter:s|={y^x for y in s}
    return s
assert set(C)<=span(mon(u)for u in range(64)if u.bit_count()==2)
assert set(A)==span(mon(u)for u in range(64)if u.bit_count()<=1)
owner=32*len(C)*len(A);assert owner==4194304
a,b=Q(28,5),Q(29,20)
N=62+a*a+b*b
assert a>b>1 and N==Q(38185,400)
def gate(d,m=N,k=N):
    assert d<=0 or 4*d*d<=m*k,(d,m,k)
for dot in (62-2*a*b,16+a*a+b*b-2,16+2*(a*b+1),32+2*(a+b)):
    gate(dot)
gate(32+a+b,N,64)
for dot in (12+a+b,13+a,15+b,16):gate(dot,N,16)
gate(12+a+b,N,16)  # nonzero tail sparse: body has 14 positions.
gate(2+a+b,N,4)  # any 4-site sign support, irrespective of its tail.
gate(0,N,1)  # remaining poles are in the tail only.
# Atlas body amplitude sqrt(2), donor tail zero. At d=4, zero or
# negative change at the large site; at d>=5, total change <=a+b.
for dot in (24+b-1,22+a+b):
    assert 8*dot*dot<=N*72
sf={65:10741,66:11418,67:11792,68:12529}  # pinned predecessor four-site sets.
h={65:2,66:4,67:6,68:24}
values={}
for n,banks,heads in rows:
    assert heads==44*h[n]*2048
    neutral=80*2**22+banks*2048+sf[n]*16+2*n
    assert neutral+126*h[n]*2048-352*2048*(n==68)==previous[n]
    values[str(n)]=neutral+heads+owner-128
out={'status':'PASS','lower_bounds':values}
assert out==json.loads((P/'claims.json').read_text())
print(json.dumps(out,sort_keys=True))
