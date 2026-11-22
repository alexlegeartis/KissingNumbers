#!/usr/bin/env python3
"""Exact Pless weight-15 census and QR chord-rank replay."""
from hashlib import sha1
import json
import os
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

P=Path(__file__).resolve().parent
QR='dim68-91-lattice-code-continuation/ternary_qr48_20261117/certificate.json'
OUT='dim90-92-hexagonal-phase-layers/hexagon_sparse15_phase7812_20261117/outer_certificate.json'
INC='dim88-89-signed-row-designs/claims.json'
SOURCES={QR:'0bd3ae5f20417ca0785fca0994f23146605e1531',
         OUT:'d03993eaee704683e6189ff301df6ce13930dc7d',
         INC:'1d5278cb39b4ff76de27ccbd04e493f41bf45681',
         'dim88-89-signed-row-designs/verify.py':'e78cc0f43716a58ebd254455058f378329b2b912',
         'dim88-axial-pless/verify.py':'c0dbe398a28857ba38e1661e90143e26387afc88'}
def pin(s):
    raw=(P.parent/s).read_bytes()
    assert sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()==SOURCES[s],s
    return raw
def rank(A):
    A=[a[:] for a in A]
    r=0
    for j in range(len(A[0])):
        v=next((i for i in range(r,len(A)) if A[i][j]),None)
        if v is None:continue
        A[v],A[r]=A[r],[(A[v][j]*x)%3 for x in A[v]]
        for i in range(r+1,len(A)):
            m=A[i][j]
            A[i]=[(x-m*y)%3 for x,y in zip(A[i],A[r])]
        r+=1
    return r
def main():
    assert __debug__
    for s in SOURCES:pin(s)
    fixed=[0,3]
    old=[1,2,24,28]
    axes=[9,30,39,44]
    QRg=json.loads(pin(QR))['generator_Q48']
    G=json.loads(pin(OUT))['outer_generator']
    assert len(G)==len(QRg)==24
    assert all(len(r)==48 and set(r)<={0,1,2} for r in G+QRg)
    assert rank(QRg)==24 and all(sum(x*y for x,y in zip(a,b))%3==0 for a in QRg for b in QRg)
    for x in (old,axes):
        assert len(set(fixed+x))==6 and rank([[r[i] for i in fixed+x] for r in QRg])==6
    payload=''.join(' '.join(map(str,r))+'\n' for r in G)
    H=[]
    with TemporaryDirectory() as job:
        exe=Path(job)/'replay'
        subprocess.run([os.environ.get('CXX','g++'),'-std=c++17','-O3',str(P/'replay.cpp'),'-o',str(exe)],check=True)
        for ax in (old,axes):
            r=subprocess.run([str(exe)],input=payload+' '.join(map(str,ax))+'\n',text=True,capture_output=True,check=True)
            z=r.stdout.splitlines()
            assert len(z)==3
            h=z[-1].split()
            assert h[0]=='histogram' and h[6:9]==['count','194304','sparse']
            H.append(([int(v) for v in h[1:6]],int(h[9])))
    assert H==[([37468,80272,58152,16784,1628],959965380),([37392,80576,57696,17088,1552],960535152)]
    prev=json.loads(pin(INC))[0]['lower_bound']
    dense=2*2**4*3**18
    assert dense==12397455648
    gain=H[1][1]-H[0][1]
    assert prev==13357793920
    out=dict(dimension=88,fixed=fixed,chord_axes=axes,previous=prev,
             dense=dense,previous_sparse=H[0][1],sparse=H[1][1],histogram=H[1][0],
             gain=gain,lower_bound=prev+gain)
    assert out==json.loads((P/'claims.json').read_text())
    print(json.dumps(out,indent=2))
if __name__=='__main__':
    main()
