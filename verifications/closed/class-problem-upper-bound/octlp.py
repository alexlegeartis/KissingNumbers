import numpy as np, collections
from scipy.optimize import linprog
from ortools.sat.python import cp_model
import os, sys, time
O=np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'octads.npy')).astype(np.int64); n=len(O)
M=O@O.T
rel=np.zeros((n,n),dtype=np.int8)          # 0 self, 1 meet4, 2 meet2, 3 meet0
rel[M==4]=1; rel[M==2]=2; rel[M==0]=3; rel[M==8]=0
k=[int((rel[0]==i).sum()) for i in range(4)]
Pint=np.zeros((4,4,4),dtype=np.int64)
rng=np.random.default_rng(0); ok=True
for h in range(4):
    idx=np.nonzero(rel[0]==h)[0]; tab=None
    for _ in range(4):
        v=int(rng.choice(idx))
        t=np.array([[int(((rel[0]==i)&(rel[v]==j)).sum()) for j in range(4)] for i in range(4)])
        if tab is None: tab=t
        elif not (tab==t).all(): ok=False
    Pint[h]=tab
print("octad scheme valencies",k,"association scheme:",ok)
B=[Pint[:,i,:].astype(float) for i in range(4)]
w,V=np.linalg.eig(B[1]); o=np.argsort(-w.real); V=V.real[:,o]
Pm=np.zeros((4,4))
for c in range(4):
    v=V[:,c]/V[0,c]
    for i in range(4): Pm[c,i]=(B[i]@v)[0]
Pm=np.round(Pm,6)
mult=np.array([n/sum(Pm[c,i]**2/k[i] for i in range(4)) for c in range(4)])
Q=np.array([[mult[j]*Pm[j,i]/k[i] for j in range(4)] for i in range(4)])
print("P=\n",Pm,"\nmult",np.round(mult,4))
free=[2,3]
A=np.array([[-Q[i,j] for i in free] for j in range(1,4)]); b=np.array([Q[0,j] for j in range(1,4)])
r=linprog(-np.ones(2),A_ub=A,b_ub=b,bounds=[(0,None)]*2,method='highs')
print("Delsarte LP (= theta') for the meet-4 graph: %.6f  -> alpha <= %d"%(1-r.fun,int(np.floor(1-r.fun))))
r2=linprog(-np.ones(2),A_ub=A,b_ub=b,bounds=[(None,None)]*2,method='highs')
print("Lovasz theta: %.6f"%(1-r2.fun))
A4=np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'meet4.npy'))
for T in range(10,17):
    mm=cp_model.CpModel(); x=[mm.NewBoolVar('') for _ in range(n)]
    ii,jj=np.nonzero(np.triu(A4,1))
    for a,bb in zip(ii.tolist(),jj.tolist()): mm.Add(x[a]+x[bb]<=1)
    mm.Add(sum(x)>=T)
    s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=60; s.parameters.num_search_workers=4
    st=s.Solve(mm); print("  exists %d pairwise-(meet 0 or 2) octads? %s (%.0fs)"%(T,s.StatusName(st),s.WallTime()),flush=True)
    if s.StatusName(st)=="INFEASIBLE": break
    if s.StatusName(st)=="UNKNOWN": break
