"""Exact maximum-clique audit of the stated 1296 quaternionic-cell family.

Every scalar-product maximum is (a+b*sqrt(2))/8 with nonnegative integers.
The test against 1/2 is integer arithmetic, including the case a >= 4.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import itertools,json,hashlib,time
from pathlib import Path
import numpy as np
from numba import njit

ROOT=Path(__file__).resolve().parent


def imul(a,b):
    return np.r_[a[0]*b[0]-a[1:]@b[1:],a[0]*b[1:]+b[0]*a[1:]+np.cross(a[1:],b[1:])].astype(np.int64)


def exact_frames():
    powers=np.array([[2,0,0,0],[1,1,1,1],[-1,1,1,1]],dtype=np.int64)
    cs=[];ds=[]
    for p in powers:
        c=[];d=[]
        for e in 2*np.eye(4,dtype=np.int64):
            product=imul(p,e);assert np.all(product%2==0)
            ci=product//2;c.append(ci)
            product=imul(np.array([1,1,0,0]),ci);assert np.all(product%2==0)
            d.append(product//2)
        cs.append(c);ds.append(d)
    c=np.array(cs);d=np.array(ds)
    # c/2 and d/sqrt(2) are exact coordinate descriptions.
    for f in c:assert np.array_equal(f@f.T,4*np.eye(4,dtype=int))
    for f in d:assert np.array_equal(f@f.T,2*np.eye(4,dtype=int))
    return c,d


def scalar_table(c,d):
    # dot = (A+B*sqrt(2))/2; only absolute values enter sign-lift compatibility.
    table=np.zeros((6,6,4,4,2),dtype=np.int64)
    for f,g in itertools.product(range(6),repeat=2):
        if f<3 and g<3:
            dots=abs(c[f]@c[g].T);assert np.all(dots%2==0)
            table[f,g,:,:,0]=dots//2
        elif f>=3 and g>=3:
            table[f,g,:,:,0]=abs(d[f-3]@d[g-3].T)
        else:
            dots=abs(c[f]@d[g-3].T) if f<3 else abs(d[f-3]@c[g].T)
            assert np.all(dots%2==0)
            table[f,g,:,:,1]=dots//2
    return table


@njit(cache=True)
def exact_compatibility(cells,table):
    n=len(cells);adj=np.zeros((n,n),dtype=np.bool_)
    for i in range(n):
        h,fi,gi,ki,di=cells[i]
        for j in range(i):
            hh,fj,gj,kj,dj=cells[j]
            # Weight products = (A+B*sqrt(2))/4.
            wa=np.ones(3,dtype=np.int64);wb=np.zeros(3,dtype=np.int64)
            if h==hh:wa[h]=2
            else:wa[h]=0;wa[hh]=0;wb[h]=1;wb[hh]=1
            good=True
            for a in range(4):
                if not good:break
                for b in range(4):
                    if not good:break
                    cc=a^b^di
                    for aa in range(4):
                        if not good:break
                        for bb in range(4):
                            dd=aa^bb^dj
                            qa=np.array([table[fi,fj,a,aa,0],table[gi,gj,b,bb,0],table[ki,kj,cc,dd,0]])
                            qb=np.array([table[fi,fj,a,aa,1],table[gi,gj,b,bb,1],table[ki,kj,cc,dd,1]])
                            ra=0;rb=0
                            for k in range(3):
                                ra+=wa[k]*qa[k]+2*wb[k]*qb[k]
                                rb+=wa[k]*qb[k]+wb[k]*qa[k]
                            # Maximum (ra+rb*sqrt(2))/8 must be <= 1/2.
                            if (ra>4 or (ra==4 and rb>0) or
                                (ra<4 and 2*rb*rb>(4-ra)*(4-ra))):
                                good=False;break
            if good:adj[i,j]=True;adj[j,i]=True
    return adj


def find_clique(bits,target):
    nodes=0
    def visit(chosen,available):
        nonlocal nodes
        nodes+=1
        if len(chosen)==target:return chosen
        while available.bit_count()+len(chosen)>=target:
            vbit=available&-available
            available^=vbit
            v=vbit.bit_length()-1
            result=visit(chosen+[v],available&bits[v])
            if result is not None:return result
        return None
    answer=visit([],(1<<len(bits))-1)
    return answer,nodes


def main():
    start=time.monotonic();c,d=exact_frames();table=scalar_table(c,d)
    cells=[]
    for h in range(3):
        for fs in itertools.product(range(6),repeat=3):
            if fs[h]<3:continue
            for delta in range(4):cells.append([h,*fs,delta])
    cells=np.array(cells,dtype=np.int64)
    adj=exact_compatibility(cells,table)
    bits=[sum(1<<int(j) for j in np.flatnonzero(row)) for row in adj]
    control=[]
    for r in range(3):
        control.append([0,3+r,r,r,0])
        control.append([1,r,3+r,3+r,[2,1,3][r]])
    indices=[int(np.flatnonzero(np.all(cells==row,axis=1))[0]) for row in control]
    assert all(adj[i,j] for i,j in itertools.combinations(indices,2))
    seven,nodes=find_clique(bits,7)
    print('seven_cell_clique',seven,'search_nodes',nodes,flush=True)
    assert seven is None,'An exact >=968-point candidate requires immediate verification.'
    report=dict(candidate_cells=len(cells),compatible_pairs=int(adj.sum()//2),
                exact_scalar_arithmetic='nonnegative (a+b sqrt(2))/8, integer comparison',
                adjacency_sha256=hashlib.sha256(adj.tobytes()).hexdigest(),
                exhaustive_seven_clique_search=True,seven_clique_exists=False,search_nodes=nodes,
                six_clique_witness=indices,six_cells=control,
                conclusion='This precise union-of-whole-cells ansatz has maximum 840 points including 72 axes.',
                scope='Fixed six frames; xor affine lifts; independent component signs; fixed weight permutations; all 72 axial points.',
                elapsed_seconds=time.monotonic()-start)
    (ROOT/'exact_cells_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
