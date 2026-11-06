"""Independent exact angle audit. Requires NumPy and SymPy. Reuses the separately supplied geometry replay, then independently reconstructs labels and checks raw-coordinate products."""
import numpy as np, json, hashlib, argparse
from pathlib import Path
P=Path(__file__).resolve().parent; root=P
parser=argparse.ArgumentParser(); parser.add_argument('--output',type=Path,default=P/'independent_replay');args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
import verify_families as vf
_,residues,_,tiny=vf.geometry()
f=np.load(P/'configuration_43.npz');x=np.load(P/'explicit_points_43.npz');g=np.load(P/'class_geometry.npz')
bc=np.load(P/'body_complement.npz');Q=bc['q'].astype(np.int64);B=bc['b'].astype(np.int64)
cg={'body_gram':Q,'body_residues':residues,'tiny':tiny}
H=np.load(P/'tiny_bodies.npz')['H'].astype(np.int64)
assert H.shape==(30240,48) and B.shape==(40,48)
assert len(np.unique(H,axis=0))==len(H)
assert np.all(np.einsum('ij,ij->i',H,H)==2880)
assert np.array_equal(B@B.T,96*Q)
import sympy as sp
# Positive definiteness follows from B B^T and the exact inverse below.
# Reconstruct integer coefficients with an exact rational inverse, not float rounding.
J=np.array((sp.Matrix(Q.tolist()).inv()*3).tolist(),dtype=np.int64)
num=H@B.T@J; assert np.all(num%288==0); coeff=num//288
assert np.array_equal(coeff@B,H)
lookup={tuple(v%3):i for i,v in enumerate(cg['body_residues'])}; assert len(lookup)==6561
labels=np.array([lookup[tuple(v%3)] for v in coeff]); tiny=np.bincount(labels,minlength=6561); assert np.array_equal(tiny,cg['tiny'])
C=f['classes']; R=f['r3']; Z=f['head_numerators'].astype(object); d=int(f['head_denominator']); U=f['pure_directions'].astype(object)
I=x['body_indices']; T=x['head_numerators'].astype(object); D=int(x['head_denominator'])
assert Z.shape==(118,3) and T.shape==(len(I),3) and U.shape==(5,3) and d>0 and D>0
assert np.all((I>=0)&(I<len(H))) and all(a.dtype.kind in 'iu' for a in (C,R,f['weights'],f['head_numerators'],x['body_indices'],x['head_numerators']))
assert set(R.tolist())<={12,14,16,18}
assert len(set(zip(C.tolist(),R.tolist(),map(tuple,Z))))==118
assert len(set(zip(I.tolist(),map(tuple,T))))==len(I)
assert np.array_equal(f['family_residues']%3,cg['body_residues'][C]%3)
for c,r,w in zip(C,R,f['weights']):
 mu=int(g['mu'][c]); assert (r,mu%24) in {(12,0),(14,16),(16,8),(18,0)}
 expected={12:495,14:9180,16:116235-40*int(tiny[c]),18:1092000 if c==0 else 1084590}[int(r)]
 assert int(w)==expected
 assert r!=12 or c!=0
FN=R.astype(object)*d*d+np.sum(Z*Z,axis=1); TN=10*D*D+np.sum(T*T,axis=1); UN=np.sum(U*U,axis=1)
assert all(v>0 for v in np.r_[FN,TN,UN])
def test(p,n,m):
 assert np.all((p<=0)|(p*p<=n*m))
for i,c in enumerate(C):
 diff=((g['digits'][C]-g['digits'][c])%3)@g['powers']; beta=R+R[i]-18+g['mu'][diff]//4
 beta[(C==c)&(R==R[i])]=2*R[i];beta[i]=2*R[i]-18
 test(beta.astype(object)*d*d+2*Z@Z[i],FN,FN[i]);test(2*U@Z[i],UN,FN[i])
p=2*U@U.T; np.fill_diagonal(p,0);test(p,UN[:,None],UN)
refined=0
for i,(bi,t) in enumerate(zip(I,T)):
 diff=((g['digits'][C]-g['digits'][labels[bi]])%3)@g['powers']; mu=g['mu'][diff]//4;beta=R-8+mu
 for j in np.flatnonzero(mu==8):
  # Direct 48-coordinate products establish absence of the tiny difference witness.
  if not np.any(H[labels==diff[j]]@H[bi]==144*(int(R[j])-20)):
   beta[j]-=6;refined+=1
 test(beta.astype(object)*D*d+2*Z@t,FN,TN[i]);test(2*U@t,UN,TN[i])
 if i%500==0:print('cross',i,flush=True)
for i in range(0,len(I),32):
 j=min(i+32,len(I)); dots=H[I[i:j]]@H[I].T; assert np.all(dots%144==0)
 p=(dots//144).astype(object)*D*D+2*T[i:j]@T.T;p[np.arange(j-i),np.arange(i,j)]=0
 test(p,TN[i:j,None],TN)
print('PASS',flush=True)
report={'status':'PASS','total':int(f['weights'].sum())+len(I)+len(U),'family_points':int(f['weights'].sum()),'tiny_points':len(I),'pure_points':len(U),'all_tiny_ordered_pairs':len(I)**2,'refined_cross_bounds':refined,'dimension':43,'tiny_body_uniqueness':len(H),'scope':'Independent exact raw-coordinate all tiny pairs and tiny/family, family/family, pure pairs; complete populations and gluing/lattice minimum remain inherited mathematical premises. No full enumeration of large families.'}
report['sha256']={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [P/'configuration_43.npz',P/'explicit_points_43.npz']}
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
