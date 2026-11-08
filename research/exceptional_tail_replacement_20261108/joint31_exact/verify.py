"""Replay the finite 31D no-gain certificate without any optimization call."""
import os
for k in('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='2'
from pathlib import Path
import numpy as np,json,hashlib,sympy as sp,time
from scipy.sparse import coo_matrix,csc_matrix
out=Path(__file__).resolve().parent/'data';meta=json.loads((out/'lp_constraints.json').read_text());source=out/'incumbent251200.npz';assert hashlib.sha256(source.read_bytes()).hexdigest()==meta['source_sha256'];pool=out/'repriced_pool.npz';assert hashlib.sha256(pool.read_bytes()).hexdigest()==meta['pool_sha256'];X=np.load(source)['points'].astype(np.int16);Z=np.load(pool)['candidates'].astype(np.int16);n=len(Z);cert=np.load(out/'lp_dual_certificate.npz');old=cert['old_source_indices'];oldmap={int(j):i for i,j in enumerate(old)};assert len(oldmap)==len(old);start=time.time()
assert X.shape==(251200,31) and Z.shape==(3762,31)
assert np.all((X.astype(np.int64)**2).sum(1)==128) and np.all((Z.astype(np.int64)**2).sum(1)==128)
assert len(np.unique(X,axis=0))==len(X) and len(np.unique(Z,axis=0))==len(Z)
assert not ({r.tobytes() for r in X}&{r.tobytes() for r in Z})
assert max(int(abs(X).max()),int(abs(Z).max()))<=8 and 31*8*8<2**24
# Exact bounded-integer binary32 products: all non-variable old points clear all new points.
keep=np.ones(len(X),bool);keep[old]=False;K=X[keep].astype(np.float32)
for at in range(0,n,128):assert np.max(Z[at:at+128].astype(np.float32)@K.T)<=64
rr=[];cc=[];vv=[];rhs=[];Zi=Z.astype(np.int32)
for o,C in meta['owner_cliques']:
 assert 0<=o<len(old) and len(set(C))==len(C)
 U=Zi[C];assert np.all(U@X[old[o]].astype(np.int32)>64);G=U@U.T;assert np.all(G[np.triu_indices(len(C),1)]>64)
 k=len(rhs);rr.extend([k]*(len(C)+1));cc.extend(C+[n+o]);vv.extend([1]*len(C)+[-1]);rhs.append(0)
for C in meta['global_cliques']:
 assert len(set(C))==len(C);U=Zi[C];G=U@U.T;assert np.all(G[np.triu_indices(len(C),1)]>64)
 k=len(rhs);rr.extend([k]*len(C));cc.extend(C);vv.extend([1]*len(C));rhs.append(1)
rank_rows={tuple(r['tail']):r for r in json.loads((out/'rank_certificates.json').read_text())}
for f in meta['fibers']:
 new=f['candidates'];oldf=f['incumbent_source'];cap=f['capacity'];q=Z[new[0],24:];assert int(q.astype(np.int64)@q)==64
 assert np.array_equal(np.flatnonzero((Z[:,24:]==q).all(1)),new);assert np.array_equal(np.flatnonzero((X[:,24:]==q).all(1)),oldf)
 U=np.r_[X[oldf,:24],Z[new,:24]].astype(np.int64);assert np.all((U*U).sum(1)==64)
 if cap<48:
  row=rank_rows[tuple(q)];assert row['capacity']==cap;N=sp.Matrix(row['null_basis']);assert N.shape==(24-cap//2,24);assert cap%2==0
  assert sp.Matrix(U.tolist())*N.T==sp.zeros(len(U),N.rows);assert (N*N.T).det()>0
 else:assert cap==48
 sub=[j for j in oldf if j in oldmap];k=len(rhs);rr.extend([k]*(len(new)+len(sub)));cc.extend(new+[n+oldmap[j] for j in sub]);vv.extend([1]*len(new)+[-1]*len(sub));rhs.append(cap-len(oldf))
A=coo_matrix((vv,(rr,cc)),shape=(len(rhs),n+len(old)),dtype=np.int64).tocsc();stored=csc_matrix((cert['matrix_data'],cert['matrix_indices'],cert['matrix_indptr']),shape=tuple(cert['matrix_shape']));assert (A-stored).nnz==0;assert np.array_equal(rhs,cert['rhs'])
y=cert['dual_inequalities'];u=cert['dual_upper_bounds'];scale=int(cert['scale']);assert scale>0 and np.all(y<=0) and np.all(u<=0);c=np.r_[-np.ones(n,np.int64),np.ones(len(old),np.int64)];assert int(abs(y).max(initial=0))*int(abs(A).sum(axis=0).max())<2**62
residual=c*scale-A.T@y-u
terms=[int(b)*int(v) for b,v in zip(rhs,y)]+list(map(int,u))+[min(0,int(v)) for v in residual]
num=sum(terms);assert num>-scale
report={'status':'PASS','candidate_count':n,'incumbent_count':len(X),'old_deletion_variables':len(old),'validated_constraints':len(rhs),'exact_lower_bound_for_negative_gain':f'{num}/{scale}','integer_gain_upper_bound':0,'source_sha256':meta['source_sha256'],'pool_sha256':meta['pool_sha256'],'seconds':time.time()-start,'scope':'All exchanges using this explicit 3762-point pool; includes every old-point deletion. No bound outside this finite pool.'};print(json.dumps(report),flush=True)
