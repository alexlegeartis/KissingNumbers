"""Optimize within independent triple components. Check all selected edges on replay."""
from pathlib import Path
import subprocess, struct, tempfile
import numpy as np

def improve(M,Q,G,rows,occ):
 n=len(rows);k=M.N[0];assert k==13824
 r=np.array(rows)
 ray,c,t=r.T
 dots=G[ray[:,None],ray[None,:]]
 edges=np.where(c[:,None]!=c[None,:],0,np.where(t[:,None]==t[None,:],2,dots!=0)).astype('u1')
 np.fill_diagonal(edges,0)
 W=np.empty((n,k),np.dtype([('g','<u4'),('r','<u2'),('b0','<u2'),('b1','<u2'),('b2','<u2'),('i','<u4')]))
 assert W.dtype.itemsize==16
 B=M.B
 hh=np.array([M.bidx[int(h)]for h in M.H],dtype='u2')
 qq=np.min(np.arange(4096,dtype='u2')[:,None]^hh,axis=1)
 cc,qids=np.unique(qq,return_inverse=True);assert len(cc)==32
 for i,(_,c,t) in enumerate(rows):
  B01=np.searchsorted(B,Q[t,:,1:])
  B2=B01[:,0]^B01[:,1]^M.bidx[M.H[t]]
  assert np.array_equal(B[B2],Q[t,:,1]^Q[t,:,2]^M.H[t])
  A=M.key[0][t]>>24
  W[i]['g']=((c*4096+A.astype('u4'))<<10)|(qids[B01[:,0]]<<5)|qids[B01[:,1]]
  W[i]['r']=i
  W[i]['b0']=B01[:,0];W[i]['b1']=B01[:,1];W[i]['b2']=B2
  W[i]['i']=i*k+np.arange(k)
 seed=np.array([v[0]for v in occ],dtype='u1')
 with tempfile.TemporaryDirectory() as tmp:
  exe=Path(tmp)/'components'
  subprocess.run(['g++','-O3','-std=c++17',str(Path(__file__).with_suffix('.cpp')),'-o',str(exe)],check=True)
  inp=struct.pack('<2I',n,n*k)+edges.tobytes()+W.tobytes()+seed.tobytes()
  p=subprocess.run([str(exe)],input=inp,capture_output=True,check=True)
  print(p.stderr.decode().strip(),flush=True)
  actual=np.frombuffer(p.stdout,dtype='u1').reshape(n,k)
  assert np.max(actual)<=1
  import os
  if os.getenv('EXPORT_FINAL'):
   np.savez_compressed(os.environ['EXPORT_FINAL'],triple=actual.astype(bool),row=np.array(rows),full48=np.array([a[1]for a in occ]),full66=np.array([a[2]for a in occ]))
  gain=int(actual.sum())-int(seed.sum())
  assert gain>=0
  return gain
