"""Reconstruct and exactly verify the discovered 5D control's rational Gram.
Rounding proposes data; every mathematical condition is then checked exactly.
"""
from pathlib import Path
from fractions import Fraction
import json,time
import numpy as np
from scipy.linalg import qr
if not __debug__:raise SystemExit('Exact certification refuses optimized Python mode.')
P=Path(__file__).resolve().parent

def det(A):
    A=[[Fraction(v) for v in row] for row in A];answer=Fraction(1)
    for j in range(len(A)):
        p=next((i for i in range(j,len(A)) if A[i][j]),None)
        if p is None:return Fraction(0)
        if p!=j:A[p],A[j]=A[j],A[p];answer=-answer
        pivot=A[j][j];answer*=pivot
        for i in range(j+1,len(A)):
            a=A[i][j]/pivot
            for k in range(j+1,len(A)):A[i][k]-=a*A[j][k]
    return answer

def inverse(A):
    n=len(A);A=[[Fraction(x) for x in row]+[Fraction(i==j) for j in range(n)] for i,row in enumerate(A)]
    for j in range(n):
        p=next(i for i in range(j,n) if A[i][j]);A[p],A[j]=A[j],A[p]
        z=A[j][j];A[j]=[v/z for v in A[j]]
        for i in range(n):
            if i!=j:
                z=A[i][j];A[i]=[a-z*b for a,b in zip(A[i],A[j])]
    return [row[n:] for row in A]

def verify(data):
    H=data['gram_twice'];ids=data['basis_indices'];n=len(H);r=len(ids)
    if n!=40 or r!=5 or any(len(row)!=n for row in H):raise ValueError('shape')
    if any(H[i][i]!=2 for i in range(n)):raise ValueError('norm')
    if any(H[i][j]!=H[j][i] or H[i][j]>1 for i in range(n) for j in range(i)):raise ValueError('angle/symmetry')
    B=[[H[i][j] for j in ids] for i in ids]
    minors=[det([row[:k] for row in B[:k]]) for k in range(1,r+1)]
    if any(v<=0 for v in minors):raise ValueError('not positive definite')
    inv=inverse(B);A=[[H[i][j] for j in ids] for i in range(n)]
    C=[[sum(Fraction(A[i][a])*inv[a][b] for a in range(r)) for b in range(r)] for i in range(n)]
    if any(sum(C[i][b]*A[j][b] for b in range(r))!=H[i][j] for i in range(n) for j in range(n)):raise ValueError('Gram does not reconstruct')
    return {'points':n,'exact_rank':r,'norm_squared':2,'maximum_distinct_inner_product':max(H[i][j] for i in range(n) for j in range(i)),'pairs_checked':n*(n-1)//2,'basis_leading_principal_minors':[str(v) for v in minors],'positive_semidefinite':'H=A*B^-1*A^T with B positive definite (exact Fraction arithmetic)','outcome':'Known 40-point 5D control independently discovered numerically then exactly reconstructed; no new SOTA.'}

def main():
    started=time.perf_counter();X=np.load(P/'d5_c5_40_polished_control.npz')['vectors']
    H=np.rint(2*(X@X.T)).astype(int);_,_,piv=qr(X.T,pivoting=True)
    data={'gram_twice':H.tolist(),'basis_indices':piv[:5].astype(int).tolist(),'origin':'Rounded numerical Gram proposal; exact verification supplies all proof claims.'}
    report=verify(data);report['seconds']=time.perf_counter()-started
    (P/'d5_40_exact_gram.json').write_text(json.dumps(data,indent=2)+'\n')
    (P/'d5_40_exact_receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
    bad=json.loads(json.dumps(data));bad['gram_twice'][0][0]=3
    try:verify(bad)
    except ValueError:pass
    else:raise RuntimeError('norm mutation not rejected')
    bad=json.loads(json.dumps(data));bad['gram_twice'][0][1]=bad['gram_twice'][1][0]=2
    try:verify(bad)
    except ValueError:pass
    else:raise RuntimeError('angle mutation not rejected')
if __name__=='__main__':main()
