from pathlib import Path
import numpy as np,sympy as s,json,hashlib,math
ROOT=Path(__file__).resolve().parents[2];pool=ROOT/'research/exceptional_tail_replacement_20261108/tetrads29/pool.npz';source=ROOT/'verifications/improved/dim27-31-joint-refills/data/coordinates29_231697.npz';C=np.load(pool)['candidates'];X=np.load(source)['points'];H=X[np.sum(X[:,24:].astype(np.int64)**2,axis=1)==64];out=[]
for t in np.unique(C[:,24:],axis=0):
 A=np.r_[H[np.all(H[:,24:]==t,axis=1),:24],C[np.all(C[:,24:]==t,axis=1),:24]].astype(np.int64)
 if np.linalg.matrix_rank(A.astype(float))==24:continue
 ns=s.Matrix(A.tolist()).nullspace();V=[]
 for v in ns:
  den=math.lcm(*(int(x.q) for x in v));w=[int(x*den) for x in v];assert all(sum(int(a)*b for a,b in zip(row,w))==0 for row in A);V.append(w)
 W=s.Matrix(V);det=(W*W.T).det();assert det!=0;rank=24-len(V)
 out.append({'tail':t.tolist(),'rank_upper':rank,'capacity':2*rank,'rows':len(A),'null_basis_rows':V,'null_gram_determinant':str(det)})
rep={'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'pool_sha256':hashlib.sha256(pool.read_bytes()).hexdigest(),'certificates':out,'scope':'Exact rank upper certificates for the union of incumbent and candidate bodies on each listed tail; lower rank is not asserted.'};pool.with_name('tail_rank_certificates.json').write_text(json.dumps(rep,indent=2)+'\n');print([(x['tail'],x['capacity']) for x in out])
