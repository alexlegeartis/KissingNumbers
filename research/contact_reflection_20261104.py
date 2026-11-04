"""Exact contact-reflection candidate generation and coupled exchange search.

A finite search, not a global kissing-number claim. The input is a valid,
common-norm integer code in an integer metric whose positive definiteness is independently verified. For a pair
with <x,y>=N/2, the root reflection x-y has norm N; for <x,y>=-N/2,
x+y has norm N. The generator also considers their negatives. A bounded seed
set, pair limit, candidate limit, and optional blocker cutoff define precisely
which part of this search is covered. Screening and solver outputs are checked
with exact integer products before a witness is returned.
"""
from __future__ import annotations
import itertools
import numpy as np
from scipy.optimize import milp, Bounds, LinearConstraint
from scipy.sparse import coo_matrix


def _dot(a, g, b):
    return sum(int(a[i])*int(g[i][j])*int(b[j]) for i in range(len(a)) for j in range(len(b)))


def search(points, metric, seeds=None, pair_limit=100000, candidate_limit=5000,
           blocker_limit=None, seconds=60):
    def as_integer(v):
        if isinstance(v,(bool,np.bool_)) or not isinstance(v,(int,np.integer)):
            raise ValueError('coordinates and metric must be exact integers')
        return int(v)
    x=[tuple(map(as_integer,row)) for row in points];g=[tuple(map(as_integer,row)) for row in metric]
    if not x or not g or any(len(row)!=len(g) for row in g) or any(len(v)!=len(g) for v in x):
        raise ValueError('nonempty points and square dimension-matched metric required')
    if any(g[i][j]!=g[j][i] for i in range(len(g)) for j in range(len(g))):
        raise ValueError('metric must be symmetric')
    if len(set(x))!=len(x):raise ValueError('duplicate source points')
    N=_dot(x[0],g,x[0])
    if N<=0:raise ValueError('positive common norm required')
    if pair_limit<1 or candidate_limit<1 or seconds<=0:raise ValueError('positive resource limits required')
    if blocker_limit is not None and blocker_limit<0:raise ValueError('negative blocker limit')
    if seeds is None: seeds=list(range(len(x)))
    if len(set(seeds))!=len(seeds) or any(type(i)!=int or not 0<=i<len(x) for i in seeds):
        raise ValueError('invalid seed indices')
    old=set(x);candidates=[];seen=set();pairs=0;contacts=0;truncated=None
    dim=len(g);maxx=max(abs(v) for row in x for v in row);maxg=max(abs(v) for row in g for v in row)
    # Candidate entries can double; this conservative bound also covers the
    # final factor-two threshold test. Fall back to Python integers otherwise.
    fast=4*dim*dim*maxx*maxx*maxg < 2**62
    if fast:
        xx=np.asarray(x,dtype=np.int64);gg=np.asarray(g,dtype=np.int64);xg=xx@gg
        if np.any(np.einsum('ij,ij->i',xg,xx)!=N):raise ValueError('positive common norm required')
    elif any(_dot(v,g,v)!=N for v in x):raise ValueError('positive common norm required')
    for i in seeds:
        all_dots=xg@xx[i] if fast else None
        for j in range(len(x)):
            if i==j:continue
            if pairs>=pair_limit:truncated='pair_limit';break
            pairs+=1;dot=int(all_dots[j]) if fast else _dot(x[i],g,x[j])
            if 2*abs(dot)!=N:continue
            contacts+=1;sgn=1 if dot>0 else -1
            v=tuple(a-sgn*b for a,b in zip(x[i],x[j]))
            for y in (v,tuple(-a for a in v)):
                if y in old or y in seen:continue
                if _dot(y,g,y)!=N:raise ArithmeticError('reflection norm invariant failed')
                seen.add(y);candidates.append(y)
                if len(candidates)>=candidate_limit:truncated='candidate_limit';break
            if truncated:break
        if truncated:break
    # The source's validity is an explicit input assumption; every returned
    # modified output is checked below rather than trusting the solver.
    blockers=[];keep=[]
    for v in candidates:
        bad=tuple(map(int,np.flatnonzero(2*(xg@np.asarray(v,dtype=np.int64))>N))) if fast else tuple(i for i,u in enumerate(x) if 2*_dot(v,g,u)>N)
        if blocker_limit is None or len(bad)<=blocker_limit:
            keep.append(v);blockers.append(bad)
    candidates=keep;m=len(candidates);n=len(x)
    receipt={'source_points':n,'norm':N,'seed_count':len(seeds),'pairs_examined':pairs,
             'contact_pairs':contacts,'candidates_after_filter':m,'truncated':truncated,
             'blocker_limit':blocker_limit,'arithmetic':'bounded int64' if fast else 'Python integers','source_validity':'assumed; independently verify input','metric_positive_definiteness':'assumed; independently verify input'}
    if not m:return {'gain':0,'removed':[],'added':[],'receipt':receipt}
    blocked_ids=sorted({j for bad in blockers for j in bad}); pos={j:k for k,j in enumerate(blocked_ids)}; nb=len(blocked_ids)
    rows=[];cols=[];values=[];upper=[]
    def edge(a,b,old):
        r=len(upper);rows.extend((r,r));cols.extend((a,b));values.extend((1,-1 if old else 1));upper.append(0 if old else 1)
    for i,bad in enumerate(blockers):
        for j in bad:edge(i,m+pos[j],True)
    for i in range(m):
        for j in range(i):
            if 2*_dot(candidates[i],g,candidates[j])>N:edge(i,j,False)
    A=coo_matrix((values,(rows,cols)),shape=(len(upper),m+nb)).tocsr()
    result=milp(np.r_[-np.ones(m),np.ones(nb)],integrality=np.ones(m+nb),
                bounds=Bounds(0,1),constraints=LinearConstraint(A,-np.inf,upper),
                options={'time_limit':float(seconds)})
    receipt['solver_status']=str(result.message)
    if result.x is None:return {'gain':0,'removed':[],'added':[],'receipt':receipt}
    chosen=[i for i,v in enumerate(result.x[:m]) if v>.5]
    removed=sorted({j for i in chosen for j in blockers[i]})
    if len(chosen)<=len(removed):return {'gain':0,'removed':[],'added':[],'receipt':receipt}
    removed_set=set(removed)
    for i in chosen:
        if any(j not in removed_set for j in blockers[i]):raise ArithmeticError('unremoved blocker')
    if any(2*_dot(candidates[i],g,candidates[j])>N for i,j in itertools.combinations(chosen,2)):
        raise ArithmeticError('solver proposal failed exact candidate check')
    receipt['output_scope']='all changed pairs checked exactly; old-old pairs require verified source'
    return {'gain':len(chosen)-len(removed),'removed':removed,
            'added':[list(candidates[i]) for i in chosen],'receipt':receipt}


# Standalone correctness controls. Run this file to execute them.
import unittest

class ContactSearch(unittest.TestCase):
 def test_recovers_hexagon(self):
  G=[[2,-1],[-1,2]]; full=[[1,0],[0,1],[1,1],[-1,0],[0,-1],[-1,-1]]
  out=search(full[:-1],G);self.assertEqual(out['gain'],1);self.assertEqual(out['added'],[[-1,-1]])
 def test_shared_deletion_enables_gain(self):
  X=[[-1,2,1],[2,1,1],[1,-2,1],[-2,-1,-1],[-1,-1,2],[1,1,-2]]
  old=[tuple(a+b for a,b in zip(x,y)) for i,x in enumerate(X) for y in X[:i]
       if sum(a*b for a,b in zip(x,y))==-3]
  self.assertEqual(set(old),{(2,-1,-1)})
  out=search(X,[[1,0,0],[0,1,0],[0,0,1]])
  self.assertEqual(out['gain'],1);self.assertEqual(len(out['removed']),1);self.assertEqual(len(out['added']),2)
 def test_python_integer_fallback(self):
  K=10**30;G=[[2*K,-K],[-K,2*K]];X=[[1,0],[0,1],[1,1],[-1,0],[0,-1]]
  out=search(X,G);self.assertEqual(out['gain'],1);self.assertEqual(out['receipt']['arithmetic'],'Python integers')
 def test_rejects_inexact_input(self):
  with self.assertRaises(ValueError):search([[1.1,0]],[[1,0],[0,1]])
 def test_seed_and_resource_limits(self):
  G=[[2,-1],[-1,2]]; x=[[1,0],[0,1],[1,1],[-1,0],[0,-1]]
  out=search(x,G,seeds=[0],pair_limit=1);self.assertEqual(out['receipt']['pairs_examined'],1);self.assertEqual(out['receipt']['truncated'],'pair_limit')
 def test_rejects_bad_source(self):
  with self.assertRaises(ValueError):search([[1,0],[1,0]],[[1,0],[0,1]])
  with self.assertRaises(ValueError):search([[1,0],[0,2]],[[1,0],[0,1]])


if __name__=='__main__': unittest.main()
