"""Exploratory Gamma72 profile systems and one-coordinate head compression.

Numerical discovery only. No output here is a kissing certificate.
"""
from pathlib import Path
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import prod
import json
import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "verifications/improved/dim68-69-gamma72-cross-sections"
N = 6218175600
M = {0:2603658750,1:1512243200,2:280928256,3:13959168,4:127800}

def build(k):
    if k == 2:
        c=json.loads((ROOT/'verifications/improved/dim70-71-gamma72-cross-sections/affine70/certificate.json').read_text())
        cells=np.array(c['cells'],dtype=np.int64)
        A=np.zeros((len(c['rows']),len(cells)))
        for i,row in enumerate(c['rows']):
            for j,v in row:A[i,j]=v
        scale=np.max(np.abs(A),axis=1);A/=scale[:,None]
        b=np.array(c['rhs'],dtype=float)/scale/N
        r=linprog(np.zeros(len(cells)),A_eq=A,b_eq=b,bounds=(0,None),method='highs')
        assert r.success,r.message
        return dict(A=A,b=b,weights=r.x*N,gram=np.array(c['gram']),cells=cells,groups=np.arange(len(cells)))
    w = json.loads((BASE/f"weyl_averaged/certificate{k}.json").read_text())
    c = json.loads((BASE/f"affine68_69/certificate{k}.json").read_text())
    G = np.array(w["gram"], dtype=np.int64)
    cells = np.array(w["cells"], dtype=np.int64)
    groupid = np.empty(len(cells),dtype=int)
    for i,gg in enumerate(w["groups"]):
        groupid[gg] = i
    roots = [np.array(a)@G for a in product(range(-2,3),repeat=k) if np.array(a)@G@np.array(a)==8]

    @lru_cache(None)
    def gaussian(e):
        if not sum(e): return 1
        i = next(i for i,v in enumerate(e) if v)
        f = list(e); f[i] -= 1
        total = 0
        for j,v in enumerate(f):
            if v:
                f[j] -= 1
                total += int(G[i,j])*v*gaussian(tuple(f))
                f[j] += 1
        return total
    rows=[];rhs=[]
    for tag in c["tags"]:
        if tag[0]=="mom":
            powers=tag[1];m=sum(powers)//2
            vals=np.prod(cells**np.array(powers),axis=1)
            b=Fraction(N*8**m*gaussian(tuple(powers)),prod(72+2*j for j in range(m)))
            b-=sum(prod(int(x)**a for x,a in zip(root,powers)) for root in roots)
        elif tag[0]=="marg":
            i,v=tag[1:];vals=(cells[:,i]==v).astype(int)
            b=Fraction(M[abs(v)]-sum(root[i]==v for root in roots))
        else: raise ValueError(tag)
        row=np.bincount(groupid,weights=vals,minlength=len(w["groups"]))
        scale=max(np.max(np.abs(row)),1)
        rows.append(row/scale);rhs.append(float(b)/scale/N)
    A=np.array(rows);b=np.array(rhs)
    r=linprog(np.zeros(A.shape[1]),A_eq=A,b_eq=b,bounds=(0,None),method="highs")
    assert r.success, r.message
    return dict(A=A,b=b,weights=r.x*N,gram=G,cells=cells,groups=groupid)

def bound(system,chosen):
    objective=np.bincount(system["groups"][chosen],minlength=system["A"].shape[1])
    r=linprog(objective,A_eq=system["A"],b_eq=system["b"],bounds=(0,None),method="highs")
    return r.fun*N if r.success else None

if __name__=="__main__":
    for k in (3,4):
        S=build(k);np.savez(Path(__file__).with_name(f"full_k{k}.npz"),**S)
        print("K",k,"rank",np.linalg.matrix_rank(S["A"]),S["A"].shape)
        H=np.linalg.inv(S["gram"])
        for i,x in enumerate(S["weights"]):
            gg=np.flatnonzero(S["groups"]==i);a=S["cells"][gg[0]]
            e=np.zeros(len(S["weights"]));e[i]=1
            lo=linprog(e,A_eq=S["A"],b_eq=S["b"],bounds=(0,None),method="highs").fun*N
            hi=-linprog(-e,A_eq=S["A"],b_eq=S["b"],bounds=(0,None),method="highs").fun*N
            print(a.tolist(),"orbit",len(gg),"each",round(x),"range",round(lo),round(hi),"q",a@H@a)
