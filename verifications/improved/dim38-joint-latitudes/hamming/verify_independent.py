from pathlib import Path
from itertools import combinations,product
from collections import Counter
from fractions import Fraction as F
import numpy as np,json,hashlib,time
P=Path(__file__).parent;T=P;S=P.parent;src=S/'source/source';c=json.loads((T/'certificate.json').read_text());old=json.loads((S/'refill/certificate.json').read_text());start=time.time();assert __debug__;assert c['dimension']==38 and c['format']=='orthogonal_hamming_38_v1' and c['new_point_count']==512 and c['previous_size']==old['claimed_size']==611522
assert hashlib.sha256((S/'refill/certificate.json').read_bytes()).hexdigest()==c['baseline_certificate_sha256']
for name,sha in c['source_sha256'].items():assert hashlib.sha256((src/name).read_bytes()).hexdigest()==sha
assert hashlib.sha256((T/'new_points.npz').read_bytes()).hexdigest()==c['new_points_sha256']
# Fresh Golay/Leech reconstruction, without importing either supplied verifier.
basis=[[int(v) for v in line.strip()] for line in (src/'golay_basis.txt').read_text().splitlines()];G=np.array(basis,dtype=np.int64);assert G.shape==(12,24);code=np.array([np.bitwise_xor.reduce(G[[i for i in range(12) if mask>>i&1]],axis=0) if mask else np.zeros(24,dtype=np.int64) for mask in range(4096)]);words={tuple(row) for row in code};assert len(words)==4096 and Counter(code.sum(1))=={0:1,8:759,12:2576,16:759,24:1};assert np.all((G@G.T)%2==0)
lines=[]
for i,j in combinations(range(24),2):
 for sign in [-1,1]:
  v=[0]*24;v[i]=4;v[j]=4*sign;lines.append(v)
for w in code:
 if int(w.sum())==8:
  support=np.flatnonzero(w)
  for mask in range(128):
   if bin(mask).count('1')%2:continue
   v=np.zeros(24,dtype=int);v[support[0]]=2
   for i,k in enumerate(support[1:]):v[k]=-2 if mask>>i&1 else 2
   lines.append(v.tolist())
 for i in range(24):
  v=1-2*w;v=v.copy();v[i]*=-3
  if v[0]>0:lines.append(v.tolist())
A=np.array(sorted(set(map(tuple,lines))),dtype=np.int64);assert A.shape==(98280,24) and np.all(np.sum(A*A,axis=1)==32)
def member(X):
 for row in X:
  parity=int(row[0]%2);assert np.all(row%2==parity);assert tuple(((row-parity)//2)%2) in words;assert int(row.sum())%8==4*parity
member(A)
W=np.loadtxt(src/'dim14_directions.txt',delimiter=',',dtype=np.int64);labels=np.loadtxt(src/'class_labels.txt',dtype=int);tri=np.loadtxt(src/'triples.txt',dtype=int);E=np.loadtxt(src/'equator_lines_6855.txt',dtype=np.int64);assert W.shape==(1932,14) and labels.shape==(98280,) and tri.shape==(644,3);assert np.all((labels>=0)&(labels<644));assert np.all(np.sum(W*W,axis=1)==8);assert len(np.unique(W,axis=0))==1932 and np.array_equal(np.sort(tri.ravel()),np.arange(1932));assert np.all(W[tri].sum(1)==0)
def pairmax(X,absolute=False):
 # Cauchy bounds every absolute partial sum by max squared norm <32767.
 norm=np.sum(X.astype(np.int64)**2,axis=1);assert norm.max()<32767;Y=X.astype(np.int16);best=-99999
 for lo in range(0,len(Y),128):
  block=Y[lo:lo+128]@Y.T
  if absolute:block=np.abs(block)
  block[np.arange(len(block)),np.arange(lo,lo+len(block))]=-32000;best=max(best,int(block.max()))
 return best
assert pairmax(W)<=4
for lab in range(644):
 B=A[labels==lab];assert len(B)>0 and (len(B)<2 or pairmax(B,True)<=8)
assert E.shape==(6855,24) and np.all(np.sum(E*E,axis=1)==48);member(E);canon={tuple(row if row[np.flatnonzero(row)[0]]>0 else -row) for row in E};assert len(canon)==6855;assert pairmax(E,True)<=24
# Retained axis and original-cap compatibility, including restored caps.
Q=np.array(old['rotation_numerator'],dtype=np.int64);assert old['rotation_denominator']==1 and np.array_equal(Q@Q.T,np.eye(14,dtype=int));axis=W[old['axis_keep']]@Q.T;assert len(axis)==1568 and len(np.unique(axis,axis=0))==1568;assert np.max(np.abs(W@axis.T))**2<=48;assert old['additional_coordinate_axes']==list(range(14));assert np.max(W*W)<=6 and np.max(axis*axis)<=2
# Expand the entire retained grid, checking its point set and all pairs.
grid_ids=np.array(old['body_line_indices']);grid_z=np.array(old['doubled_head_lines'],dtype=np.int64);grid=np.array([np.r_[sb*A[i],2*sh*z] for i,z in zip(grid_ids,grid_z) for sb,sh in product([-1,1],repeat=2)]);assert grid.shape==(6272,38);assert len(np.unique(grid,axis=0))==6272;assert np.all(np.sum(grid*grid,axis=1)==128) and np.all(np.sum(grid[:,:24]**2,axis=1)==32);grid_pair_max=pairmax(grid);assert grid_pair_max<=64;assert np.max(np.abs(grid[:,24:]@W.T))<=16 and np.max(np.abs(grid[:,24:]@axis.T))<=16 and np.max(grid[:,24:]**2)<=32
# Original cap/grid bound: (1+sqrt(2))^2 <=6 follows from8<=9.
assert 8<=9
# Fresh construction of the Hamming layer, independent of metadata labels.
z=np.load(T/'new_points.npz');B=z['body_integer'].astype(np.int64);Y=z['head_integer'].astype(np.int64);roots=[]
for a,b in [(14,15),(16,17),(18,19),(20,21)]:
 for sign in [-1,1]:
  u=np.zeros(24,dtype=np.int64);u[a]=4;u[b]=4*sign;roots.append(u)
roots=np.array(roots);assert np.array_equal(roots,np.array(c['body_roots']));assert np.array_equal(roots@roots.T,32*np.eye(8,dtype=int));ids=np.array(c['body_line_indices']);assert np.array_equal(A[ids],roots)
expected=set()
for half in [0,1]:
 for word in range(128):
  syndrome=0
  for k in range(7):
   if word&(1<<k):syndrome^=k+1
  y=np.zeros(14,dtype=int);y[half*7:half*7+7]=[(-1 if word&(1<<k) else 1) for k in range(7)]
  for sg in [-1,1]:expected.add(tuple(np.r_[sg*roots[syndrome],y]))
assert B.shape==(512,24) and Y.shape==(512,14);assert len(expected)==512 and set(map(tuple,np.c_[B,Y]))==expected;assert np.all(np.sum(B*B,axis=1)==32) and np.all(np.sum(Y*Y,axis=1)==7)
# All new points have raw norm128: (4/3)*32+(256/21)*7=128.
assert F(4,3)*32+F(256,21)*7==128
J=28*(B@B.T)+256*(Y@Y.T);assert np.all(np.diag(J)==2688);np.fill_diagonal(J,-1);new_pair_max=int(J.max());assert new_pair_max<=1344
removed=set(map(int,c['removed_cap_body_lines']));gridowners=set(map(int,grid_ids));assert len(gridowners)==14 and removed==gridowners|set(map(int,ids)) and len(removed)==22
retained=np.array([i for i in range(len(A)) if i not in removed]);restored=np.array(sorted(set(old['removed_cap_body_lines'])-removed));assert len(restored)==4
bodycap=int(np.abs(roots@A[retained].T).max());headcap=int(np.abs(Y@W.T).max());assert bodycap<=16 and headcap<=4
# Independent rational proof of 2/sqrt3+4sqrt(2/21)<=sqrt6.
assert F(17,10)**2<=3 and F(31,100)**2>=F(2,21) and F(17,7)**2<=6;assert F(20,17)+F(31,25)<=F(17,7)
assert np.max(np.abs(roots@A[restored].T))<=16;assert np.max(np.abs(A[np.unique(grid_ids)]@A[restored].T))==0
assert np.max(np.abs(B@grid[:,:24].T))==0;grid_tail_max=int(np.abs(Y@grid[:,24:].T).max());assert (16*grid_tail_max)**2<=64**2*21
axismax=int(np.abs(Y@axis.T).max());assert axismax**2<=21;assert F(2,21)<=F(1,4)
eqmax=int(np.abs(B@E.T).max());assert eqmax**2<=1152
# Cross-component distinctness follows from body squared norm fractions.
fractions=[F(2,3),F(1),F(0),F(1,4),F(1,3)];assert len(set(fractions))==5
counts={'caps':6*len(retained),'equator':2*len(E),'axes':len(axis)+28,'grid':len(grid),'new':len(B)};total=sum(counts.values());assert total==611638==c['claimed_size'];report={'status':'PASS','dimension':38,'total':total,'gain_over_611522':total-611522,'counts':counts,'restored_cap_lines':restored.tolist(),'new_new_pairs':512*511//2,'max_new_pair_numerator_denominator':[new_pair_max,21],'max_grid_pair_numerator':grid_pair_max,'new_cap_maxima':[bodycap,headcap],'new_grid_tail_max':grid_tail_max,'new_axis_max':axismax,'new_equator_max':eqmax,'certificate_sha256':hashlib.sha256((T/'certificate.json').read_bytes()).hexdigest(),'new_points_sha256':hashlib.sha256((T/'new_points.npz').read_bytes()).hexdigest(),'scope':'Fresh independent Golay/Leech reconstruction, cap class/head checks, equator pairs, retained-grid pairs, all new interactions and restored-cap bounds. Retained cap-cap global minimum bound uses the source elementary Golay lattice minimum32 proof, not a fresh 9.6-billion Leech-pair enumeration. No external optimality/priority claim.','seconds':time.time()-start};(P/'independent_verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
