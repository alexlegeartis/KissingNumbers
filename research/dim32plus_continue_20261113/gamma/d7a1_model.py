"""Mixed D7 x A1 moments, including exact membership of boundary profiles."""
from pathlib import Path
from fractions import Fraction as F
from math import gcd,prod
import json
import numpy as np
import sympy as sp
from scipy.optimize import linprog
import extension_model as m
from verify_extension import inverse,multiply

HERE=Path(__file__).resolve().parent
N=m.N

def build(dim):
    emb=json.loads((HERE/'d7a1_embedding.json').read_text())
    U=emb['normals']
    source=json.loads((m.REPO/'verifications/improved/dim68-69-gamma72-cross-sections/weyl_averaged/embedding.json').read_text())
    G=multiply(multiply(U,source['catalogue_gram']),list(map(list,zip(*U))))
    assert G==emb['gram']
    assert G[-1]==[0]*7+[8] and all(row[-1]==0 for row in G[:-1])
    B=[[1,-1,0,0,0,0,0],[0,1,-1,0,0,0,0],[0,1,1,0,0,0,0],[-1,0,0,1,0,0,0],
       [0,0,0,-1,1,0,0],[0,0,0,0,-1,1,0],[0,0,0,0,0,-1,1]]
    W=multiply([[2*x for x in row] for row in inverse(B)],U[:7])
    assert all(x.denominator==1 for row in W for x in row)
    W=[[int(x) for x in row] for row in W]
    reps=m.representatives(7)
    groups=[];boundary=[]
    for h in reps:
        for c in range(5):
            q16=sum(x*x for x in h)+2*c*c
            if q16<128:groups.append((h,c))
            elif q16==128:
                count=0;examples=[]
                for hh in m.orbit(h):
                    for cc in ((-c,c) if c else (0,)):
                        numerator=[sum(hh[i]*W[i][j] for i in range(7))+2*cc*U[7][j] for j in range(72)]
                        if all(x%16==0 for x in numerator):
                            count+=1
                            if len(examples)<2:examples.append([x//16 for x in numerator])
                boundary.append(dict(h=list(h),c=c,population=count,examples=examples))
    profiles,heights=m.selected_profiles(dim)
    marginal=m.a2_marginal_weights(profiles) if dim==70 else m.d4_marginal_weights(profiles)
    size=len(next(iter(marginal)))
    sizes=[m.orbit_size(h)*(2 if c else 1) for h,c in groups]
    objective=[m.objective(h,marginal)*(2 if c else 1) for h,c in groups]
    known=[(h,0) for h in m.normal_roots(7)]+[((0,)*7,-8),((0,)*7,8)]
    known_obj=sum((marginal.get(tuple(sorted(map(abs,h[:size]))),F(0)) for h,c in known),F(0))
    for z in boundary:known_obj+=F(z['population'],m.orbit_size(z['h']))*m.objective(z['h'],marginal)
    rows=[];rhs=[];tags=[]
    for total in range(6):
        for a in range(total+1):
            for part in m.partitions(total-a,7):
                pp=tuple(2*x for x in part);cp=2*a
                row=[m.orbit_moment(h,pp)*(2 if c else 1)*c**cp for h,c in groups]
                b=F(N*128**(total-a)*64**a*prod(prod(range(1,p,2)) for p in pp)*prod(range(1,cp,2)),
                    prod(72+2*j for j in range(total)))
                b-=sum(prod(x**p for x,p in zip(h,pp))*c**cp for h,c in known)
                for z in boundary:
                    b-=F(z['population'],m.orbit_size(z['h']))*m.orbit_moment(tuple(z['h']),pp)*z['c']**cp
                row=[x*b.denominator for x in row];bn=b.numerator;g=gcd(*row,bn)
                rows.append([x//g for x in row]);rhs.append(bn//g);tags.append([list(pp),cp])
    return dict(dimension=dim,groups=[[list(h),c] for h,c in groups],sizes=sizes,
                rows=rows,rhs=rhs,tags=tags,objective=[str(x) for x in objective],
                known_objective=str(known_obj),boundary=boundary,profiles=profiles,heights=heights)

def solve(dim):
    d=build(dim)
    A=sp.Matrix(d['rows']);b=sp.Matrix(d['rhs']);choose=list(A.T.rref()[1])
    A=A[choose,:];b=b[choose,:]
    for key in ('rows','rhs','tags'):d[key]=[d[key][i] for i in choose]
    sizes=np.asarray(d['sizes'],float);An=np.asarray(A,float)/sizes[None,:];bn=np.asarray(b,float).ravel()/N
    scales=np.maximum(abs(bn),1e-15);cn=np.asarray([float(F(x)) for x in d['objective']])/sizes
    lp=linprog(cn,A_eq=An/scales[:,None],b_eq=bn/scales,bounds=(0,None),method='highs',
               options=dict(threads=1,presolve=False,primal_feasibility_tolerance=1e-9,dual_feasibility_tolerance=1e-9))
    assert lp.success,lp.message
    active=[i for i,x in enumerate(lp.x) if x>1e-11]
    c=sp.Matrix([sp.Rational(x) for x in d['objective']]);sol=next(iter(sp.linsolve((A[:,active].T,c[active,:]))))
    free=set().union(*(x.free_symbols for x in sol));y=sp.Matrix([x.subs({f:0 for f in free}) for x in sol])
    slack=c-A.T*y
    exact=None
    if all(x>=0 for x in slack):
        exact=(b.T*y)[0]+sp.Rational(d['known_objective'])
        d.update(exact_value=str(exact),exact_bound=int(sp.ceiling(exact)),dual_exact=[str(x) for x in y],slack_exact=[str(x) for x in slack])
    d.update(value=lp.fun*N+float(F(d['known_objective'])),primal=(lp.x/sizes).tolist())
    (HERE/f'd7a1_compact{dim}.json').write_text(json.dumps(d,indent=2)+'\n')
    print('D7A1',dim,'groups',len(sizes),'rank',len(choose),'boundary',[(z['h'],z['c'],z['population']) for z in d['boundary']],
          'float',d['value'],'exact',exact,flush=True)

if __name__=='__main__':
    for dim in (68,69,70):solve(dim)
