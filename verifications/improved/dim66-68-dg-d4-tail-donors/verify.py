#!/usr/bin/env python3
"""Exact two-sided D4-tail donor insertion in dimensions 66–68."""
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from tempfile import TemporaryDirectory
from itertools import product,combinations
from subprocess import run,check_output
from functools import reduce
from operator import xor
from contextlib import redirect_stdout
from io import StringIO
import json,runpy,sys
assert __debug__
here=Path(__file__).resolve().parent
repo=here.parent
prior=repo/'dim65-68-dg-diagonal-triples'
base=repo/'dim65-68-dg-retained-atlas'
sys.path[:0]=[str(prior),str(base)]
def check(path,digest): 
    b=path.read_bytes();assert sha256(b).hexdigest()==digest,(path,digest)
    return b
# Pin the predecessor's cited stock and its exact geometry files.
for n,h in json.loads((here/'predecessor.json').read_text()).items():
    check(prior/n,h)
with redirect_stdout(StringIO()):
    ancestor=runpy.run_path(str(prior/'verify.py'))
E=ancestor['E']
old=ancestor['out']['lower_bounds']
assert old=={'65':424177492,'66':434139558,'67':453138696,'68':483338010}
oldnorm=ancestor['verify_diagonal'].N
oldden=ancestor['verify_diagonal'].D
oldv=ancestor['verify_diagonal'].V
profiles=ancestor['verify_diagonal'].P
assert all(p[-1]==p[-2] and abs(p[-1])==40547 for p in profiles)
data=json.loads((here/'selection.json').read_text())
labels=sorted({r[2]for choices in data.values()for r in choices})
assert labels and all(0<=l<2**20 for l in labels)
def mon(u):return sum(1<<i for i in range(64)if i&u==u)
cub=[mon(u)for u in range(64)if u.bit_count()==3]
quad=[mon(u)for u in range(64)if u.bit_count()==2]
lo=[mon(u)for u in range(64)if u.bit_count()<=2]
def word(code,arr):
    assert 0<=code<1<<len(arr)
    return reduce(xor,(x for i,x in enumerate(arr)if code>>i&1),0)
bwords={c:word(c,cub)for c in labels}
C,A=E['C'],E['A']
assert len(C)==1024 and len(A)==128
# A literal quadratic translate of each DG stock preserves the minimum
# distance 24. Distinct translates used within a given cubic are disjoint.
assert len(A)==len(set(A))and all(w.bit_count()>=24 for w in A if w)
def gate(d,n,m):
    return d<=0 or 4*d*d<=n*m
t=F(31,10)
N=64+2*t*t
assert N==F(4161,50)
assert gate(36,N,64) and gate(16+2*t*t,N,N) and gate(t,N,1)
assert gate(2+2*t,N,4)
assert all(gate(x,N,16)for x in (16,14+t,12+2*t,14+2*t-2))
sels={}
for dim,entries in data.items():
    dim=int(dim);sel=[];keys=set()
    assert 66<=dim<=68
    for tag,parity,c,a,b in entries:
        assert 0<=tag<1<<(dim-64)and tag.bit_count()==2 and parity in (0,1)
        assert (tag,parity) not in keys;keys.add((tag,parity))
        i,j=[q for q in range(4)if tag>>q&1]
        tail=tuple(int(k==i)+(1-2*parity)*int(k==j)for k in range(4))
        for s,shift in ((1,a),(-1,b)):
            rot=tuple(x*s for x in tail)
            for nr,p,d in zip(oldnorm,profiles,oldv):
                m=F(36)+F(d,oldden)+t*F(p[-1],oldden)*(rot[0]+rot[1])
                assert gate(m,N,F(nr,oldden**2)),(dim,tag,rot,m)
            sel.append((tag,parity,c,word(shift,quad),rot))
    assert len(set(x[-1]for x in sel))==len(sel)
    for c in labels:
        sets=[{a^shift for a in C} for _,_,v,shift,_ in sel if v==c]
        assert len(set().union(*sets))==len(sets)*1024
    sels[dim]=sel
# The predecessor has already certified its angular inequalities. Extract the
# complete masks and exact total parities from the SAME deterministic C++ replay.
with TemporaryDirectory()as d:
    T=Path(d)
    for name,contents in E['sources'].items():(T/name).write_bytes(contents)
    assignments=(T/'assignment.txt').read_text().splitlines()
    assignments[0]='46';heads=sorted(ancestor['heads'])
    assignments[-1]=' '.join(map(str,heads))
    (T/'assignment.txt').write_text('\n'.join(assignments)+'\n')
    phase=list((T/'phase.txt').read_text().strip())
    for p in ancestor['positions']:phase[p]='1' if phase[p]=='0'else'0'
    (T/'phase.txt').write_text(''.join(phase)+'\n')
    (T/'recover.txt').write_text(''.join(f'{p}\n'for p in ancestor['recovered']))
    (T/'atlas46.inc').write_bytes((prior/'atlas46.inc').read_bytes())
    src=(prior/'verify_atlas46.cpp').read_text()
    q='for(size_t i=0;i<B.size();i++)if(B[i].t==3){'
    zz="zz<<i<<' '<<B[i].b"
    assert src.count(q)==1 and src.count(zz)==1
    src=src.replace(q,'for(size_t i=0;i<B.size();i++){').replace(zz,"zz<<B[i].t<<' '<<i<<' '<<B[i].b")
    (T/'banks.cpp').write_text(src)
    run(['g++','-O3','-march=native','-mbmi2','-std=c++20',str(T/'banks.cpp'),'-o',str(T/'banks')],check=True)
    rs=check_output([str(T/'banks')],cwd=T,text=True)
    assert [(int(z.split()[1]),int(z.split()[2]))for z in rs.splitlines()if z.startswith('ATLAS')]==[(65,40414),(66,45181),(67,54363),(68,68275)]
    banks=[]
    for l in (T/'tag3.txt').read_text().splitlines():
        tag,i,b,r,phase,*H=map(int,l.split())
        assert i==len(banks) and len(H)==r and H[:2]==[0,65535]
        banks.append((tag,b,(phase>>1)&1))
    assert len(banks)==68275
    run(['g++','-O3','-march=native','-std=c++20',str(here/'clearance.cpp'),'-o',str(T/'clearance')],check=True)
    proof=check_output([str(T/'clearance'),*map(str,labels)],cwd=T,text=True)
parents={};atlas={};pairs={}
for ln in proof.splitlines():
    k,*l=ln.split();r=list(map(int,l))
    if k=='PAIR':a,b,v=r;pairs[tuple(sorted((a,b)))]=v
    elif k=='PARENT':a,b,ct=r;assert ct==335544320;parents[a]=b
    else:
        assert k=='HEAD'
        a,h,v=r;atlas[a,h]=v
assert set(parents)==set(labels) and min(parents.values())>=14
assert len(atlas)==126*len(labels)and len(pairs)==len(labels)*(len(labels)-1)//2
out={}
for dim,sel in sels.items():
    for a,b in combinations(sel,2):
        d=16 if a[2]==b[2] else pairs[tuple(sorted((a[2],b[2])))]
        h=sum(x*y for x,y in zip(a[4],b[4]))
        assert gate(64-2*d+t*t*h,N,N),(dim,a,b,d,h)
    bad_banks=set()
    for tag,pa,c,_,_ in sel:
        for i,(bt,w,s)in enumerate(banks):
            if bt!=tag:continue
            assert not any((w&v).bit_count()%2 for v in lo), (tag,i)
            if ((w&bwords[c]).bit_count()%2)^pa==s:bad_banks.add(i)
    tails=[]
    for j in range(dim-64):
        for s in (-2,2):tails.append(tuple(s if j==k else 0 for k in range(4)))
    if dim==68:tails.extend(product((-1,1),repeat=4))
    bad_atlas=set()
    for s,h,u in product(sel,heads,tails):
        d=atlas[s[2],h]
        p=32-2*d+t*sum(a*b for a,b in zip(s[4],u))
        if not gate(p,N,36):bad_atlas.add((h,u))
    # Remaining sparse supports use 0/2/4 tail marks, and the magnitude
    # tests above cover every tag except the exact donor tag. Each such
    # retained tag is paid a forced one-coordinate disagreement of ≥1.
    assert all(tag==0 or tag.bit_count()==2 or tag==15 for tag,_,_ in banks)
    restored=old[str(dim)]
    gain=len(sel)*len(C)*len(A)-(len(bad_banks)+len(bad_atlas))*len(A)*16
    assert gain>0
    out[str(dim)]=dict(lower_bound=restored+gain,gain=gain,layers=len(sel),
        dropped_banks=len(bad_banks),dropped_atlas=len(bad_atlas),
        unique_cubic_labels=len({r[2]for r in sel}))
assert sys.argv[1:] in ([],['--write'])
if sys.argv[1:]:
    (here/'result.json').write_text(json.dumps(out,indent=2)+'\n')
else:assert out==json.loads((here/'result.json').read_text())
summary={'status':'PASS','lower_bounds':{k:v['lower_bound'] for k,v in out.items()}}
assert summary==json.loads((here/'claims.json').read_text())
print(json.dumps(out,sort_keys=True))
print(json.dumps(summary,sort_keys=True))
