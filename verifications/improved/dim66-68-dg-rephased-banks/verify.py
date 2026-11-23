#!/usr/bin/env python3
"""Exact affine sparse-bank sign trades and off-axis DG donors."""
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
labels=sorted({r[0]for choices in data.values()for r in choices})
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
translation=json.loads((here/'translates.json').read_text())
pre=repo/'dim66-68-dg-d4-tail-donors'/'claims.json'
last=json.loads(check(pre,'14f576f2af0db3e1f73a1b8978fd7dd88d72cdac3be8617ff40762d8e632bd5c'))
assert last['lower_bounds']=={'66':434387366,'67':454068488,'68':485390106}
sels={}
for dim,entries in data.items():
    dim=int(dim);sel=[]
    assert 66<=dim<=68
    for c,q,*v in entries:
        assert len(v)==4 and all(not v[i]for i in range(dim-64,4))
        t=tuple(F(a,10)for a in v)
        norm=64+sum(a*a for a in t)
        assert norm>64
        assert gate(36,norm,64) and gate(norm-48,norm,norm)
        assert gate(sum(sorted([1]*64+[abs(a)for a in t],reverse=True)[:4]),norm,4)
        for nr,p,d in zip(oldnorm,profiles,oldv):
            m=F(36)+F(d,oldden)+F(p[-1],oldden)*(t[0]+t[1])
            assert gate(m,norm,F(nr,oldden**2)),(dim,c,t,m)
        sel.append((c,word(q,quad),t,norm))
    assert len(set(x[2]for x in sel))==len(sel)
    for c in labels:
        sets=[{a^shift for a in C} for v,shift,_,_ in sel if v==c]
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
    marker=' for(auto& z:B)z.h.complete();'
    assert src.count(marker)==1
    extra=r'''
  {ofstream zz2("guard_masks.txt");
  for(auto [ij,ab]:guards)zz2<<"E "<<ij.first<<" "<<ij.second<<" "<<ab.first<<" "<<ab.second<<"\n";
  for(size_t i=0;i<B.size();i++)for(auto [v,b]:HG[i])zz2<<"H "<<i<<" "<<v<<"\n";}
'''
    src=src.replace(marker,extra+marker)
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
    rows=[]
    for line in (T/'guard_masks.txt').read_text().splitlines():
        ty,*v=line.split()
        v=list(map(int,v))
        assert len(v)==(4 if ty=='E' else 2) and ty in ('E','H')
        rows.append((ty,*v))
    assert sum(x[0]=='E' for x in rows)==31512
    triple=[]
    P=[int(v)for v in (T/'cubic80.txt').read_text().split()]
    new_lines=(prior/'triple_assignment.txt').read_text().splitlines()
    assert new_lines[0]=='21'
    for l in new_lines[1:22]:
        c,*I=map(int,l.split())
        triple.append((word(P[c],cub),sum(1<<i for i in I)))
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
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def par(a,b):return (a&b).bit_count()%2
for dim,sel in sels.items():
    for a,b in combinations(sel,2):
        d=16 if a[0]==b[0] else pairs[tuple(sorted((a[0],b[0])))]
        assert gate(64-2*d+dot(a[2],b[2]),a[3],b[3]),(dim,a,b,d)
    ks=translation[str(dim)]
    X=dict(ks)
    assert len(ks)==len(X) and all(0<=i<len(banks)and 0<v<65536 for i,v in ks)
    for ty,i,j,*x in rows:
        if ty=='E':
            a,b=x
            assert par(X.get(i,0),a)==par(X.get(j,0),b)
        else:assert not par(X.get(i,0),j)
    new_banks=[]
    for i,(tag,w,s)in enumerate(banks):
        q=s^par(X.get(i,0),65535)
        assert tag==0 or tag.bit_count()==2 or tag==15
        if tag==0:assert q==0
        if tag==3:assert all((w&a).bit_count()<2 or par(w,b)==q for b,a in triple)
        new_banks.append((tag,w,q))
    bad_banks=set()
    for c,_,t,N in sel:
        for i,(tag,w,s) in enumerate(new_banks):
            if tag>=1<<(dim-64):continue
            T=[a for j,a in enumerate(t)if tag>>j&1]
            v=16-len(T)+sum(map(abs,T))
            if gate(v,N,16):continue
            # Orthogonality of the body restriction to RM(2,6) fixes
            # its sign parity on the entire donor quadratic translate.
            assert not any(par(w,b) for b in lo), (dim,i,tag)
            p=par(w,bwords[c])^sum(a<0 for a in T)%2
            if s==p or not gate(v-2*min([F(1),*map(abs,T)]),N,16):
                bad_banks.add(i)
    tails=[]
    poles=[]
    for j in range(dim-64):
        for s in (-2,2):
            u=tuple(s if j==k else 0 for k in range(4))
            tails.append(u);poles.append(tuple(x//2 for x in u))
    if dim==68:tails.extend(product((-1,1),repeat=4))
    bad_atlas=set()
    bad_poles=set()
    for (c,_,t,N),h,u in product(sel,heads,tails):
        d=atlas[c,h];p=32-2*d+dot(t,u)
        if not gate(p,N,36):bad_atlas.add((h,u))
    for (_,_,t,N),u in product(sel,poles):
        if not gate(dot(t,u),N,1):bad_poles.add(u)
    count=old[str(dim)]+len(sel)*len(C)*len(A)-(len(bad_banks)+len(bad_atlas))*len(A)*16-len(bad_poles)
    assert count>last['lower_bounds'][str(dim)]
    out[str(dim)]=count
    print('PASS',dim,count,'additional',count-last['lower_bounds'][str(dim)],
          'sign translations',len(X),'lost banks/atlas/poles',len(bad_banks),len(bad_atlas),len(bad_poles),flush=True)
summary={'status':'PASS','lower_bounds':out}
assert summary==json.loads((here/'claims.json').read_text())
print(json.dumps(summary,sort_keys=True))
