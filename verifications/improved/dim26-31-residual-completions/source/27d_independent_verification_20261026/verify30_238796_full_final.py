"""Independent exact verification; integer shell construction + proved shell lemma.
Uses Python/NumPy only, no creator code imports. One CPU, <=512MiB address space.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import resource
resource.setrlimit(resource.RLIMIT_AS,(512*1024**2,512*1024**2))
import numpy as np
from pathlib import Path
import itertools,hashlib,json,time,sys
def popcount(n):
 return n.bit_count() if hasattr(n,'bit_count') else bin(n).count('1')
start=time.time();HERE=Path(__file__).resolve().parent;D=HERE.parent/'27d_local_replacements_20261026'
def stage(s):print(round(time.time()-start,2),s,flush=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
base=D/'base30_groups325.npz';extra=D/'extension30_groups325.npz'
assert sha(base)=='a0bc0ef49d0a163e84c7515d9b6787e28bec7675d1ca2af295c1cab92090138d'
assert sha(extra)=='4c48503313b62b3aa5213afdfdc1bc3ac0fc5ad3522c5bd6e37468a0e5e0cb8d'
b=np.load(base,allow_pickle=False);e=np.load(extra,allow_pickle=False)
assert set(b.files)=={'new_points','removed_leech','midpoint_owner_indices'}
assert set(e.files)=={'free_heads','base_certificate_sha256'}
for z in [b['new_points'],b['removed_leech'],b['midpoint_owner_indices'],e['free_heads']]:assert z.dtype.kind in 'iu'
old=b['new_points'].astype(np.int32);removed=b['removed_leech'].astype(np.int32);owners=b['midpoint_owner_indices'].astype(np.int32);new=e['free_heads'].astype(np.int32)
assert old.shape==(83200,30) and removed.shape==(41600,24) and owners.shape==(83200,2) and new.shape==(576,30)
assert str(e['base_certificate_sha256'])==sha(base)
# Independent bitwise encoder, not matrix multiplication or creator imports.
generator=0xC75
basis=[]
for i in range(12):
 v=generator<<i;basis.append(v|((popcount(v)&1)<<23))
words=[0]
for v in basis:words += [w^v for w in words]
assert len(set(words))==4096
weights={k:sum(popcount(w)==k for w in words) for k in range(25)}
assert {k:v for k,v in weights.items() if v}=={0:1,8:759,12:2576,16:759,24:1}
assert all(popcount(a&b)%2==0 for a in basis for b in basis)
# Linear XOR construction + nonzero minweight8 gives all pairwise distances>=8.
# Orthogonal basis implies all word-pairs have even intersection. These exact
# facts prove every shell pair bound via SHELL_PROOF.md, no external theorem.
S=np.empty((196560,24),np.int32);at=0
for i in range(24):
 for j in range(i+1,24):
  for si,sj in itertools.product((-1,1),repeat=2):
   row=np.zeros(24,np.int32);row[i]=4*si;row[j]=4*sj;S[at]=row;at+=1
for w in words:
 if popcount(w)!=8:continue
 support=[i for i in range(24) if (w>>i)&1]
 for mask in range(256):
  if popcount(mask)%2:continue
  row=np.zeros(24,np.int32)
  for k,i in enumerate(support):row[i]=-2 if(mask>>k)&1 else 2
  S[at]=row;at+=1
for w in words:
 signs=np.array([-1 if(w>>i)&1 else 1 for i in range(24)],np.int32)
 for j in range(24):
  row=signs.copy();row[j]*=-3;S[at]=row;at+=1
assert at==196560
assert np.all(np.einsum('ij,ij->i',S,S,dtype=np.int64)==32)
keys={r.tobytes() for r in S};assert len(keys)==196560
assert len({r.tobytes() for r in removed})==41600 and all(r.tobytes() in keys for r in removed)
assert np.all((0<=owners)&(owners<41600))
assert np.array_equal(2*old[:,:24],removed[owners[:,0]]+removed[owners[:,1]])
assert np.all(owners[:,0]!=owners[:,1])
# For retained shell y, each owner is deleted and distinct from y. Therefore
# y dot midpoint <=16 by the proved full-shell bound. Tail dot y is zero.
assert np.all(np.einsum('ij,ij->i',old,old,dtype=np.int64)==32)
assert len({r.tobytes() for r in old})==83200
assert np.all(np.any(old[:,24:]!=0,axis=1))
stage('BASE_SHELL_AND_MIDPOINT_CERTIFICATE_PASS')
pure=[]
for i,j in itertools.combinations(range(6),2):
 for si,sj in itertools.product((-4,4),repeat=2):
  row=np.zeros(30,np.int32);row[24+i]=si;row[24+j]=sj;pure.append(row)
pure=np.array(pure,np.int32)
assert pure.shape==(60,30)
assert np.all(np.einsum('ij,ij->i',pure,pure,dtype=np.int64)==32)
assert len({r.tobytes() for r in pure})==60
# Bound all integer accumulations: coordinate absolute values<=4;30*16=480
# fits int32, as do norms. Do not infer no overflow only from final dot values.
assert all(np.max(np.abs(x))<=4 for x in (S,removed,old,new,pure))
def pairmax(X,Y,self_check=False,chunk=128):
 worst=-100
 for a in range(0,len(X),chunk):
  P=X[a:a+chunk]@Y.T
  if self_check:
   P[np.arange(len(P)),a+np.arange(len(P))]=-100
  m=int(P.max());assert m<=16,('PAIR_FAILURE',a,m)
  worst=max(worst,m)
 return worst
oldold=pairmax(old,old,True);stage('ALL_OLD_OLD_PAIRS_PASS')
oldpure=pairmax(old,pure);purepure=pairmax(pure,pure,True)
assert np.all(np.einsum('ij,ij->i',new,new,dtype=np.int64)==32)
assert len({r.tobytes() for r in new})==576
assert not ({r.tobytes() for r in new}&{r.tobytes() for r in old})
assert not ({r.tobytes() for r in new}&{r.tobytes() for r in pure})
assert not ({r.tobytes() for r in old}&{r.tobytes() for r in pure})
assert np.all(np.any(new[:,24:]!=0,axis=1))
newold=pairmax(new,old);newpure=pairmax(new,pure);newnew=pairmax(new,new,True)
newshell=pairmax(new[:,:24],S,chunk=16);stage('ALL_NEW_PAIRS_AND_FULL_SHELL_PASS')
assert 196560-41600+83200+60==238220
assert 238220+576==238796
result={'status':'PASS','dimension':30,'base_count':238220,'new_count':576,'total':238796,'gain_over_new_base_without_free_heads':576,'comparison_original_repository_baseline':224280,'gain_over_original_repository_baseline':14516,'norm_squared':32,'max_inner_product':16,'base_sha256':sha(base),'extension_sha256':sha(extra),'old_old_max':oldold,'old_pure_max':oldpure,'pure_pure_max':purepure,'new_old_max':newold,'new_pure_max':newpure,'new_new_max':newnew,'new_full_shell_max':newshell,'full_shell_cardinality':196560,'removed_shell_unique_members':41600,'old_vs_retained_shell':'exact owner-midpoint identity + independently proved shell pair bound','shell_pair_bound':'self-contained three-family proof with enumerated linear binary code weight/parity properties; no quadratic shell scan','arithmetic':'int32 products/accumulations bounded in absolute value by480; int64 norms','scope':'Full construction validated with stated algebraic shell lemma; not a universal current-record survey','elapsed_seconds':time.time()-start,'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(HERE/'verification30_238796_full_replay.json').write_text(json.dumps(result,indent=2)+'\n');stage(json.dumps(result))

# Independent exceptional-tail replay in the same bounded process.
os.execv(sys.executable,[sys.executable,str(HERE/"apply_exact_tail_upgrade_full_final.py"),"30"])
