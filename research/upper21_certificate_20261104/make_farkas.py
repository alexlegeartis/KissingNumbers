"""Discover a rational Farkas combination for the saved exact-cut LP."""
import json
import numpy as np
import sympy as sy
from scipy.optimize import linprog
from model import DATA

q=np.load(DATA/'outer_lp.npz');A=q['A'];b=q['b']
E=np.r_[q['E'],q['row'][None,:]];e=np.r_[q['e'],972]
scale=np.maximum(1,np.max(abs(A),axis=1))
M=np.c_[A.T/scale,E.T];normalization=np.r_[np.ones(len(A)),np.zeros(len(E))]
r=linprog(-np.r_[b/scale,e],A_eq=np.r_[M,normalization[None,:]],
          b_eq=np.r_[np.zeros(A.shape[1]),1],
          bounds=[(0,None)]*len(A)+[(None,None)]*len(E),method='highs')
assert r.success and r.fun<0
support=np.flatnonzero(abs(r.x)>1e-8)
integer_matrix=np.c_[A.T,E.T]
anchor=next(int(i) for i in support if i<len(A))
extra=np.zeros(len(r.x),dtype=np.int64);extra[anchor]=1
exact=sy.Matrix(np.r_[integer_matrix,extra[None,:]][:,support])
rhs=sy.Matrix([0]*A.shape[1]+[1]);v=exact.gauss_jordan_solve(rhs)[0]
assert not any(a.free_symbols for a in v)
full=[sy.Rational(0)]*len(r.x)
for i,a in zip(support,v):full[i]=a
assert all(a>=0 for a in full[:len(A)])
assert sy.Matrix(integer_matrix)*sy.Matrix(full)==sy.zeros(A.shape[1],1)
gap=(sy.Matrix([np.r_[b,e].tolist()])*sy.Matrix(full))[0];assert gap>0
report=dict(target=972,upper_component=971,
            inequality_terms=[(int(i),str(a)) for i,a in enumerate(full[:len(A)]) if a],
            equality_terms=[(int(i),str(a)) for i,a in enumerate(full[len(A):]) if a],
            contradiction=str(gap))
(DATA/'label_farkas972.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(dict(status='exact combination found',terms=len(support),contradiction=str(gap))),flush=True)
