"""Construct an exactly checkable rational upper bound for a finite exchange."""
import argparse,json,time
from pathlib import Path
from fractions import Fraction
import numpy as np
from scipy.sparse import coo_matrix
from scipy.optimize import linprog

def cliques(adjacency):
    def visit(R,P,X):
        if not P and not X:yield tuple(sorted(R));return
        pivot=max(P|X,key=lambda u:len(P&adjacency[u]),default=None)
        for v in sorted(P-(adjacency[pivot] if pivot is not None else set())):
            yield from visit(R|{v},P&adjacency[v],X&adjacency[v]);P.remove(v);X.add(v)
    yield from visit(set(),set(range(len(adjacency))),set())

def main():
    ap=argparse.ArgumentParser();ap.add_argument('pool',type=Path);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--threads',type=int,default=1);a=ap.parse_args();start=time.monotonic();z=np.load(a.pool);C=z['candidates'].astype(np.float32);ptr=z['blocker_indptr'];idx=z['blocker_indices'];old,inv=np.unique(idx,return_inverse=True);n=len(C);m=len(old);ownerpoints=[[] for _ in old];G=C@C.T>64
    for i in range(n):
        for j in inv[ptr[i]:ptr[i+1]]:ownerpoints[j].append(i)
    rr=[];cc=[];vv=[];upper=[];descriptions=[]
    for j,points in enumerate(ownerpoints):
        points=np.asarray(points,int);local=G[np.ix_(points,points)];np.fill_diagonal(local,False);adj=[set(map(int,np.flatnonzero(row))) for row in local]
        for group in cliques(adj):
            members=points[list(group)];row=len(upper);rr.extend([row]*(len(members)+1));cc.extend(members.tolist()+[n+j]);vv.extend([1]*len(members)+[-1]);upper.append(0);descriptions.append((int(old[j]),members.tolist()))
    # Candidate-only edges can improve a fractional bound even when owner
    # cliques alone do not. They are ordinary exact conflict inequalities.
    ii,jj=np.where(np.triu(G,1))
    for i,j in zip(ii,jj):
        row=len(upper);rr.extend([row,row]);cc.extend([int(i),int(j)]);vv.extend([1,1]);upper.append(1);descriptions.append((-1,[int(i),int(j)]))
    M=coo_matrix((vv,(rr,cc)),shape=(len(upper),n+m)).tocsr();print('LP',M.shape,flush=True);res=linprog(np.r_[-np.ones(n),np.ones(m)],A_ub=M,b_ub=upper,bounds=(0,1),method='highs',options={'threads':a.threads});assert res.success,res.message;print('numerical_bound',-res.fun,res.message,flush=True)
    if -res.fun>=1-1e-8:raise SystemExit('Insufficient fractional bound')
    lhs=[Fraction(0) for _ in range(n+m)];rhs=Fraction(0);certificate=[]
    for row,value in enumerate(-res.ineqlin.marginals):
        if value<1e-9:continue
        w=Fraction(round(float(value)*(1<<24)),1<<24);owner,members=descriptions[row];certificate.append(dict(owner=owner,members=members,weight=[w.numerator,w.denominator]));rhs+=w*upper[row]
        for i in members:lhs[i]+=w
        if owner>=0:lhs[n+int(np.searchsorted(old,owner))]-=w
    bound_weights=[]
    for i,value in enumerate(-res.upper.marginals):
        if value<1e-9:continue
        w=Fraction(round(float(value)*(1<<24)),1<<24);lhs[i]+=w;rhs+=w;bound_weights.append(dict(kind='candidate' if i<n else 'deletion',index=int(i if i<n else old[i-n]),weight=[w.numerator,w.denominator]))
    for i in range(n+m):
        target=1 if i<n else -1
        if lhs[i]<target:
            w=target-lhs[i];lhs[i]+=w;rhs+=w;bound_weights.append(dict(kind='candidate' if i<n else 'deletion',index=int(i if i<n else old[i-n]),weight=[w.numerator,w.denominator]))
    assert all(v>= (1 if i<n else -1) for i,v in enumerate(lhs));assert rhs<1,rhs
    report=dict(schema_version=1,pool_file=a.pool.name,candidates=n,source_deletions=m,cliques=certificate,variable_bounds=bound_weights,exact_gain_upper_bound=[rhs.numerator,rhs.denominator],integer_gain_upper_bound=rhs.numerator//rhs.denominator,proof='All variables are nonnegative and at most one. A source-owner clique gives sum(x_i)-y_owner <= 0; a candidate-only clique gives sum(x_i) <= 1. The nonnegative rational combination dominates sum(x_i)-sum(y_j).',seconds=time.monotonic()-start)
    a.out.write_text(json.dumps(report,indent=2)+'\n');print('EXACT GAIN BOUND',str(rhs),'rows',len(certificate),'bounds',len(bound_weights),'seconds',time.monotonic()-start,flush=True)

if __name__=='__main__':main()
