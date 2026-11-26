#!/usr/bin/env python3
"""Exact 29D three-source 36-exception obstruction. Usage: python verify_29_three.py path/to/k-nic-KissingNumbers path/to/alexlegeartis-KissingNumbers
Check out k-nic 77f0860bf93a745140e6e25481d3f60773b11441 and Alex's 3c59085cc59cc8208532aa3ca23eae3ed7273ca0 or later revision preserving the pinned data. Building conflict graphs is CPU intensive. NumPy/SciPy, preferably single-thread BLAS.
"""
import os,sys,hashlib,zlib,base64
from pathlib import Path
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np
from scipy.sparse import csr_matrix,bmat
from scipy.sparse.csgraph import maximum_bipartite_matching,connected_components,breadth_first_order
u,v=(Path(x).resolve() for x in sys.argv[1:3])
up=u/'version2/certificates';v1=u/'version1/certificates/data'
def checksum(p,sha):
 assert hashlib.sha256(p.read_bytes()).hexdigest()==sha
checksum(up/'common.py','2253a0f70c1e857c993d1f0ad4ff49b3e33e0598f2d6dfb033870a14aa402dc7')
sys.path.insert(0,str(up));from common import leech_shell,lexsort
def access(A):return np.ascontiguousarray(A.astype(np.int8)).view('V29').ravel()
def isin(A,B):return np.isin(access(A),access(B))
def uniq(A):return A[np.unique(access(A),return_index=True)[1]]
def sha(A):return hashlib.sha256(lexsort(A).astype('<i2').tobytes()).hexdigest()
def load(p,s):checksum(p,s);return np.load(p,allow_pickle=False)
L=leech_shell().astype(np.int8)
def fromcert(p,s):
 z=load(p,s);keep=np.ones(len(L),bool);keep[z['deleted_leech_indices']]=False
 return np.r_[np.pad(L[keep],((0,0),(0,5))),z['added_integer']].astype(np.int8)
b2=fromcert(up/'data/d29_certificate.npz','b27ec26a1d186ed6e9b6885921abf3b2e6bbde132206fde2689dfdaa7e39f824')
b1=fromcert(v1/'d29_certificate.npz','a87126bbc8f8b9627a501ddb2cdb73020cf1d4d655b4d7c3652af1b56615b6ce')
alt=fromcert(v1/'d29_alt_certificate.npz','ded2da5a8141af3606d375df464f2ce3bfa0c45bb02204610d94d19ccdcf8cd5')
ex=load(v/'verifications/improved/dim29-takhanov-auxiliary-exchange/data/exchange29.npz','130c0c144b014e1afa2284242718918fa7ccc2a9af8dd026a10d1a6a5aa9bd3e')
tr=load(v/'verifications/improved/dim29-even-grid-exchange/trade.npz','d85b5cdf65eb128ad55104711cfc57621554de1d31005ad181757b46b2f3f9cd')
o=np.r_[alt[~isin(alt,ex['removed'])],ex['added']].astype(np.int8)
a=np.zeros((14,29),np.int16);m=np.sqrt(128/tr['squared_norms']).round().astype(int);assert np.all(m*m*tr['squared_norms']==128)
a[:,tr['coordinates']]=tr['rays']*m[:,None]
old=np.r_[o[np.all(a@o.astype(np.int16).T<=64,axis=0)],a].astype(np.int8)
assert [sha(w) for w in (b2,old,b1)]==['3675c465c219c0e542f881356774f85dc9d67cb7ad0f83cfb805d0e808138233','85d48a37bccd4aa5950c3b5279b9a673511a71ebc1eef7f0df9647fa1d32495d','e2108120a2c72a40406eef1bb11739ddf9b957a0449ca1433ed1afb97114a224']
assert all(np.all((w.astype(np.int16)**2).sum(1)==128) for w in (b2,old,b1))
def grap(A,B):
 X=A.astype(np.float32);Y=B.astype(np.float32).T.copy(order='F');nn=[];pp=np.zeros(len(A)+1,np.int32)
 # Every dot and every possible partial absolute sum is <=128 by Cauchy–Schwarz.
 for i in range(0,len(A),96):
  w=X[i:i+96]@Y>64;pp[i+1:i+1+len(w)]=w.sum(1);nn.append(np.where(w)[1].astype(np.int32))
 np.cumsum(pp,out=pp);jj=np.concatenate(nn) if nn else np.zeros(0,np.int32)
 return csr_matrix((np.ones(len(jj),np.int8),jj,pp),shape=(len(A),len(B)))
def hall(G):
 n,m=G.shape
 if not n:return np.zeros(0,int),np.zeros(0,int)
 hit=maximum_bipartite_matching(G,perm_type='column');bad=np.flatnonzero(hit<0)
 if not len(bad):return np.zeros(0,int),np.zeros(0,int)
 used=np.flatnonzero(hit>=0)
 back=csr_matrix((np.ones(len(used),np.int8),(hit[used],used)),shape=(m,n))
 from scipy.sparse import bmat
 H=bmat([[None,G,csr_matrix((n,1),dtype=np.int8)],[back,None,csr_matrix((m,1),dtype=np.int8)],[csr_matrix((np.ones(len(bad),np.int8),(np.zeros(len(bad),int),bad)),shape=(1,n)),csr_matrix((1,m),dtype=np.int8),csr_matrix((1,1),dtype=np.int8)]],format='csr')
 reach=breadth_first_order(H,n+m,directed=True,return_predecessors=False)
 l=np.sort(reach[reach<n]);r=np.sort(reach[(reach>=n)&(reach<n+m)]-n)
 assert len(l)-len(r)==len(bad)
 return l,r

def exact_no(K,m):
 N=K.shape[0];full=(1<<N)-1;A=[];nodes=0;ordr=np.argsort(np.diff(K.indptr))[::-1];K=K[ordr][:,ordr].tocsr()
 for i in range(N):
  bad=(1<<i)
  for j in K.indices[K.indptr[i]:K.indptr[i+1]]:bad|=1<<int(j)
  A.append(full^bad)
 def bb(P,k):
  nonlocal nodes
  if not k:return True
  nodes+=1;u=P;lab=[];it=[];c=0
  while u:
   c+=1;s=u
   while s:
    bit=s&-s;x=bit.bit_length()-1;u^=bit;s^=bit;s&=~A[x];lab.append(c);it.append(x)
  for x,h in zip(reversed(it),reversed(lab)):
   if h<k:return False
   if bb(P&A[x],k-1):return True
   P^=(1<<x)
  return False
 res=bb(full,m+1);return res,nodes
d=uniq(np.r_[old[~isin(old,b2)],b1[~isin(b1,b2)]])
b=b2[~(isin(b2,old)&isin(b2,b1))]
G=grap(d,b);L,R=hall(G);d,b,G=d[L],b[R],G[L][:,R].tocsr()
assert G.shape==(62907,60306)
sig=isin(d,old).astype(int)+2*isin(d,b1).astype(int)
A=np.flatnonzero(sig==1);C=np.flatnonzero(sig==2);E0=grap(d[A],d[C]).tocoo()
i,j=A[E0.row],C[E0.col];E=csr_matrix((np.ones(2*len(i),np.int8),(np.r_[i,j],np.r_[j,i])),shape=(len(d),len(d)));assert E.nnz==2*12934
deg=np.diff(G.indptr);cheap=np.flatnonzero(deg<=64);hi=np.flatnonzero(deg>64);assert len(hi)==188 and np.all(sig[hi]==3)
g=G[cheap];e=E[cheap][:,cheap]
H=bmat([[e,g],[g.T,None]],format='csr');num,lab=connected_components(H,directed=False);big=np.argmax(np.bincount(lab))
for z in range(num):
 if z==big:continue
 ii=np.flatnonzero(lab[:len(cheap)]==z);jj=np.flatnonzero(lab[len(cheap):]==z)
 if len(ii)==0:continue
 t=g[ii][:,jj];L,R=hall(t)
 if len(L):
  t=t[L][:,R];h=bmat([[e[ii[L]][:,ii[L]],t],[t.T,None]],format='csr')
  assert not exact_no(h,len(R))[0]
di=cheap[np.flatnonzero(lab[:len(cheap)]==big)];cols=np.flatnonzero(lab[len(cheap):]==big)
pair=np.frombuffer(zlib.decompress(base64.b64decode('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')),dtype='<u2').astype(np.int32).reshape(-1,2)
aa,cc=pair.T;assert np.unique(pair).size==pair.size and np.all(E[di[aa],di[cc]].A1)
sing=np.setdiff1d(np.arange(len(di)),pair.ravel())
from scipy.sparse import vstack
g=G[di][:,cols];K=vstack([g[sing],g[aa].multiply(g[cc])],format='csr');A=G[hi][:,cols]
assert max(K.shape)<100000 and np.all(maximum_bipartite_matching(K,perm_type='column')>=0)
RANK=[]
for j in range(len(hi)):
 N=A.indices[A.indptr[j]:A.indptr[j+1]];avail=np.ones(K.shape[1],bool);avail[N]=0
 mt=maximum_bipartite_matching(K[:,avail],perm_type='column');RANK.append(len(N)-int(sum(mt<0)))
assert RANK==[35]*188
N=np.unique(A.indices);avail=np.ones(K.shape[1],bool);avail[N]=0
mt=maximum_bipartite_matching(K[:,avail],perm_type='column')
assert len(N)==4534 and len(N)-int(sum(mt<0))==35
print('PASS: every improving three-source trade needs >=36 exceptional donors (out of 188); no global record proved.')
