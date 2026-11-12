"""Different Gamma72 rank-two Gram families: affine unions and scalar heads.

Numerical geometry searches are followed by exact rational dual and pair checks.
The program does not update any maintained record.
"""
from pathlib import Path
from fractions import Fraction as F
from math import gcd,prod,isqrt,sqrt
from itertools import combinations_with_replacement
import json
import subprocess
import sys
import numpy as np
import sympy as sp
from scipy.optimize import linprog

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]
sys.path.insert(0,str(REPO/'common'))
from kpoint_lp import build_system,_integralise
N=6218175600
M1={0:2603658750,1:1512243200,2:280928256,3:13959168,4:127800}
OPTIONS=dict(threads=1,presolve=False,primal_feasibility_tolerance=1e-9,dual_feasibility_tolerance=1e-9)

def build(ip):
    G=[[8,ip],[ip,8]]
    raw=json.loads((REPO/'verifications/improved/dim68-69-gamma72-cross-sections/weyl_averaged/embedding.json').read_text())
    B=raw['catalogue_gram'];j=next(j for j in range(1,72) if abs(B[0][j])==ip and B[j][j]==8)
    u=[0]*72;u[0]=1;v=[0]*72;v[j]=1 if B[0][j]>=0 else -1
    (HERE/f'pair{ip}_embedding.json').write_text(json.dumps(dict(normals=[u,v],gram=G),indent=2)+'\n')
    S=build_system(72,N,8,[1,2,3,4],M1,G)
    cells=S['cells'];idx=S['idx'];parent=list(range(len(cells)))
    def root(i):
        if parent[i]!=i:parent[i]=root(parent[i])
        return parent[i]
    def join(a,b):
        a,b=root(a),root(b)
        if a!=b:parent[max(a,b)]=min(a,b)
    # Antipodality and root translations are exact population equalities.
    # Swapping the two normals (and, for ip=0, reversing one normal) is
    # a finite orientation average; it does not assume a shell automorphism.
    for i,c in enumerate(cells):
        for d in [(-c[0],-c[1]),(c[1],c[0])]+([(-c[0],c[1])] if ip==0 else []):
            if d in idx:join(i,idx[d])
        for a,profile in S['spanmin']:
            if sum(x*y for x,y in zip(a,c))==4:
                d=tuple(x-y for x,y in zip(profile,c))
                if d in idx:join(i,idx[d])
    labels=sorted({root(i) for i in range(len(cells))})
    groups=[[i for i in range(len(cells)) if root(i)==label] for label in labels]
    gid={i:j for j,g in enumerate(groups) for i in g}
    rows=[[sum(row[i] for i in g) for g in groups] for row in S['rows']]
    sparse,rhs,_=_integralise(rows,S['rhs'])
    A=[[0]*len(groups) for row in sparse]
    for row,ent in zip(A,sparse):
        for j,x in ent:row[j]=x
    Am=sp.Matrix(A);choose=list(Am.T.rref()[1]);A=[A[i] for i in choose];rhs=[rhs[i] for i in choose]
    sizes=np.asarray(list(map(len,groups)),float);An=np.asarray(A,float)/sizes[None,:];bn=np.asarray(rhs,float)/N
    scales=np.maximum(np.max(abs(An),axis=1),abs(bn));scales=np.maximum(scales,1)
    aa=An/scales[:,None];bb=bn/scales
    moments=[]
    for j in range(len(groups)):
        cc=np.zeros(len(groups));cc[j]=1/sizes[j]
        for sign in (-1,1):
            lp=linprog(sign*cc,A_eq=aa,b_eq=bb,bounds=(0,None),method='highs',options=OPTIONS)
            assert lp.success,lp.message
            moments.append(lp.x/sizes*N)
    weights=np.mean(moments,axis=0)
    H=[[F(8,64-ip*ip),F(-ip,64-ip*ip)],[F(-ip,64-ip*ip),F(8,64-ip*ip)]]
    inner=lambda c,d:sum(F(c[i])*H[i][j]*d[j] for i in range(2) for j in range(2))
    return dict(ip=ip,G=G,H=H,cells=cells,groups=groups,gid=gid,A=A,b=rhs,aa=aa,bb=bb,sizes=sizes,weights=weights,
                tags=[S['tags'][i] for i in choose],known=list(S['known']),inner=inner)

def evaluate(M,selected,heights=None,height_scale=F(1),label=None):
    objective=[sum(i in selected for i in g) for g in M['groups']]
    lp=linprog(np.asarray(objective)/M['sizes'],A_eq=M['aa'],b_eq=M['bb'],bounds=(0,None),method='highs',options=OPTIONS)
    assert lp.success,lp.message
    result=dict(gram=M['G'],cells=M['cells'],groups=M['groups'],rows=M['A'],rhs=M['b'],tags=M['tags'],
                objective=objective,selected=selected,profiles=[M['cells'][i] for i in selected],
                heights=heights,height_scale_squared=str(height_scale),value=lp.fun*N)
    if label:
        A=sp.Matrix(M['A']);b=sp.Matrix(M['b']);c=sp.Matrix(objective)
        active=[i for i,x in enumerate(lp.x) if x>1e-10]
        sols=list(sp.linsolve((A[:,active].T,c[active,:])))
        if sols:
            free=set().union(*(x.free_symbols for x in sols[0]));y=sp.Matrix([x.subs({f:0 for f in free}) for x in sols[0]])
            slack=c-A.T*y
            if all(x>=0 for x in slack):
                exact=(b.T*y)[0]
                result.update(exact_value=str(exact),exact_bound=int(sp.ceiling(exact)),dual_exact=[str(x) for x in y],slack_exact=[str(x) for x in slack])
        (HERE/label).write_text(json.dumps(result,indent=2)+'\n')
    return result

def maximal_cliques(adj):
    n=len(adj)
    def rec(r,p,x):
        if not p and not x:
            yield r;return
        union=p|x
        pivot=max((i for i in range(n) if union>>i&1),key=lambda i:(p&adj[i]).bit_count(),default=0)
        choices=p&~adj[pivot]
        while choices:
            bit=choices&-choices;i=bit.bit_length()-1
            yield from rec(r+[i],p&adj[i],x&adj[i])
            p&=~bit;x|=bit;choices&=~bit
    yield from rec([],(1<<n)-1,0)

def affine(M):
    cells=M['cells'];inner=M['inner'];q=[inner(c,c) for c in cells]
    adj=[0]*len(cells)
    for i in range(len(cells)):
        assert q[i]<8
        for j in range(i):
            d=4-inner(cells[i],cells[j])
            if d<=0 or 4*d*d<=(8-q[i])*(8-q[j]):adj[i]|=1<<j;adj[j]|=1<<i
    cliques=list(maximal_cliques(adj))
    ranked=sorted(cliques,key=lambda cc:sum(M['weights'][M['gid'][i]] for i in cc),reverse=True)
    best=None
    for cc in ranked[:250]:
        r=evaluate(M,cc)
        if best is None or r['value']>best['value']:best=r
    result=evaluate(M,best['selected'],['0']*len(best['selected']),label=f'pair{M["ip"]}_affine70.json')
    print('AFFINE70 ip',M['ip'],'cells',len(cells),'groups',len(M['groups']),'maximal cliques',len(cliques),
          'checked',min(250,len(cliques)),'best',result.get('exact_value',result['value']),'slices',len(best['selected']),flush=True)

def heads(M):
    ip=M['ip'];cells=M['cells'];inner=M['inner'];scale=F(8,64-ip*ip)
    L=np.linalg.cholesky(np.asarray(M['H'],float));vertices=[];metadata=[]
    for i,c in enumerate(cells):
        q=inner(c,c);limit=isqrt(int(q/scale*144));z=np.asarray(c)@L
        for j in range(-limit,limit+1):
            u=F(j,12);t2=scale*u*u
            assert t2<=q
            vertices.append((i,sqrt(float((8-t2)/(8-q))),float(u)*sqrt(float(scale)),M['weights'][M['gid'][i]],*z))
            metadata.append((i,u))
    inp=HERE/f'pair{ip}_head_vertices.txt';out=HERE/f'pair{ip}_head_search.json'
    with inp.open('w') as f:
        f.write(f'{len(vertices)} 2\n')
        for row in vertices:f.write(' '.join(map(str,row))+'\n')
    subprocess.run([str(REPO/'research/dim32plus_20261113/theory/head_search'),str(inp),str(out),'1000000'],check=True)
    selected=[metadata[j] for j in json.loads(out.read_text())['vertices']]
    pp=[cells[i] for i,u in selected];uu=[u for i,u in selected];q=[inner(c,c) for c in pp]
    a2=[(8-scale*u*u)/(8-r) for u,r in zip(uu,q)]
    for i,j in combinations_with_replacement(range(len(pp)),2):
        d=4-inner(pp[i],pp[j]);r=4-scale*uu[i]*uu[j]
        assert d>=0 and r>=0 and r*r-a2[i]*a2[j]*d*d>=0,(ip,i,j)
    result=evaluate(M,[i for i,u in selected],[str(u) for i,u in selected],scale,f'pair{ip}_head71.json')
    print('HEAD71 ip',ip,'vertices',len(vertices),'best',result.get('exact_value',result['value']),'slices',len(selected),flush=True)

if __name__=='__main__':
    models=[build(ip) for ip in (0,1,2,3)]
    for M in models:affine(M)
    if '--heads' in sys.argv:
        for M in models:heads(M)
