#!/usr/bin/env python3
"""Pinned integer-leaf replay of the alternate A2 normal fibers."""
from fractions import Fraction as Q
from pathlib import Path
from itertools import combinations, combinations_with_replacement
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
IMP = HERE.parent
INPUTS = {
    'count128.cpp': ('dim70-gamma72-d4-planar-heads/count128.cpp', 'b6f8cbf2cbb6b624f957e8d009c05dc57ee3588726294e8857d0bf1c2cbd8381'),
    'leech_structure.npz': ('dim68-91-lattice-code-continuation/mixed_shell71_20261117/leech_structure.npz', '3253520e074ea6866bb39c22bf05b595945b6183b48d9ae145758c65db79752e'),
    'min4keys.bin': ('dim68-91-lattice-code-continuation/mixed_shell71_20261117/min4keys.bin', '015256fde6d37be571af8e34fe2d5f180e4cd488e0213bd5eb4ffaed414eb45d'),
    'build_table.py': ('dim68-91-lattice-code-continuation/mixed_shell71_20261117/build_table.py', '382ef81bc43e015c0fce8befb008b9e542bcfe2c5669fe2ddb58e5111e2fc357'),
}
PROFILES = [(-1, -1), (-1, 0), (-2, -1)]
TABLE_SHA = '356cce48ef46e88bf114740aeb1a3afeddc092d1afa326f5225044341bee2037'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def geometry():
    assert len(set(PROFILES)) == 3
    def prod(a, b):
        return Q(8*a[0]*b[0] - 4*(a[0]*b[1]+a[1]*b[0]) + 8*a[1]*b[1], 48)
    assert [8-prod(p,p) for p in PROFILES] == [Q(47,6), Q(47,6), Q(15,2)]
    for a,b in combinations_with_replacement(PROFILES, 2):
        x,y = 8-prod(a,a), 8-prod(b,b)
        upper = 4-prod(a,b)
        assert x>0 and y>0 and (upper<=0 or 4*upper*upper<=x*y)


def parity(v):
    return sum((int(x)&1)<<i for i,x in enumerate(v))


def reduced_rows(K,F):
    # Any LLL decision is only a speed heuristic; all changes are integral, and
    # the exact determinant and matrix identity certify unchanged lattice span.
    U=np.eye(len(K),dtype='i8')
    original=K.copy()
    def gs():
        X=np.linalg.cholesky(K@F@K.T)
        return X/X.diagonal()[None,:], X.diagonal()**2
    mu,d=gs()
    k=1
    while k<len(K):
        for j in range(k-1,-1,-1):
            z=round(mu[k,j])
            if z:
                K[k]-=z*K[j]
                U[k]-=z*U[j]
                mu,d=gs()
        if d[k]>=(.99-mu[k,k-1]**2)*d[k-1]:
            k+=1
        else:
            K[[k,k-1]]=K[[k-1,k]]
            U[[k,k-1]]=U[[k-1,k]]
            mu,d=gs()
            k=max(1,k-1)
    assert np.array_equal(K,U@original) and abs(sp.Matrix(U).det())==1
    return K


def replay(work):
    for name, (file, h) in INPUTS.items():
        source = IMP / file
        assert sha(source) == h, file
        shutil.copyfile(source, work / name)
    subprocess.run([sys.executable,str(work/'build_table.py')],check=True)
    table=work/'fiber8_table.bin'
    assert sha(table)==TABLE_SHA
    binary=work/'count128'
    subprocess.run(['g++','-O3','-std=c++17',str(work/'count128.cpp'),'-o',str(binary)],check=True)
    with np.load(work/'leech_structure.npz',allow_pickle=False) as L:
        F=L['gram']
    witness=json.loads((HERE/'normals.json').read_text())
    assert len(witness)==2 and all(len(r)==24 for r in witness)
    assert all(type(x)==int and -64<=x<=64 for r in witness for x in r)
    N=np.array(witness,dtype='i8')
    assert np.array_equal(2*N@F@N.T,np.array([[8,4],[4,8]]))
    M=F@N.T
    rows=next([a,b] for a,b in combinations(range(24),2) if abs(int(M[a,0]*M[b,1]-M[a,1]*M[b,0]))==1)
    tail=[i for i in range(24) if i not in rows]
    inv=np.array(sp.Matrix(M[rows]).inv(),np.int64)
    K=np.eye(24,dtype=np.int64)[tail]
    K[:,rows]=-M[tail]@inv
    assert np.array_equal(K@M,np.zeros((22,2),dtype='i8'))
    assert np.array_equal(K[:,tail],np.eye(22,dtype='i8'))
    K=reduced_rows(K,F)
    assert np.array_equal(K@M,np.zeros((22,2),dtype='i8'))
    H=K@F@K.T
    Hi=sp.Matrix(H).inv()
    L,D=sp.Matrix(H).LDLdecomposition(hermitian=False)
    assert L*D*L.T==sp.Matrix(H) and all(D[i,i]>0 for i in range(22))
    d=np.diag(np.array(D,float));keys=[parity(row) for row in K];L=np.array(L,float)
    source=work/'affine.input'
    values=[]
    for j in PROFILES:
        v=np.zeros(24,np.int64);v[rows]=np.array(j)@inv
        assert np.array_equal(v@M,j) and not v[tail].any() and any(t%2 for t in j)
        h=K@F@v;q=int(v@F@v)
        sh=Hi*sp.Matrix(h)
        c=q-(sp.Matrix(h).T*sh)[0]
        assert c>0
        sh=np.array(sh,float).reshape(-1)
        with source.open('w') as f:
            f.write(f'22 8 {q} {float(c):.17g} {parity(v)}\n')
            for row in H: f.write(' '.join(map(str,row))+'\n')
            for arr in [h,keys,d,sh]: f.write(' '.join(map(str,arr))+'\n')
            for row in L: f.write(' '.join(map(str,row))+'\n')
        result=subprocess.run([str(binary),str(source),str(table),'unused','8'],capture_output=True,text=True,check=True)
        data=json.loads(result.stdout)
        values.append(data['weighted'])
        print(j,data,flush=True)
    return values


def main():
    assert __debug__, 'Run without -O'
    geometry()
    info=json.loads((HERE/'claims.json').read_text())
    with tempfile.TemporaryDirectory(prefix='gamma72-a2-') as directory:
        values=replay(Path(directory))
    assert info['dimension']==70 and info['previous']==1361137216
    assert values==info['populations'] and sum(values)==info['lower_bound']==1361142400
    print('PASS',info['lower_bound'],flush=True)


if __name__=='__main__':
    main()
