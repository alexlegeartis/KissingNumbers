"""Joint amplitude optimization of every coset; no retained dense family.

Numerical search only. Every sign difference in the quotient is constrained,
plus the exact linear maxima over the retained root and odd-octad families.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.optimize import minimize
from probe21 import G,D,line_ids,eta
OUT=Path(__file__).resolve().parent

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--seed',type=int,default=20261104);ap.add_argument('--extra',type=int,default=0);ap.add_argument('--restarts',type=int,default=4);ap.add_argument('--iterations',type=int,default=1500);ap.add_argument('--fill-one',action='store_true');args=ap.parse_args()
 rng=np.random.default_rng(args.seed);bits=((np.array(D)[:,None]>>np.arange(21))&1);sign=1.-2*bits
 syndrome=np.array([eta(w) for w in D])
 extra_masks=rng.integers(1,4096,size=args.extra).tolist()
 for j,mask in enumerate(extra_masks):syndrome|=np.array([(i&mask).bit_count()%2 for i in range(4096)])<<(4+j)
 qsize=1<<(4+args.extra)
 # Accept only surjective maps. A rank loss must not count empty families.
 if len(set(syndrome))!=qsize:raise ValueError('rank-deficient random quotient')
 octads=np.array([[(w>>i)&1 for i in range(21)] for w in G if w.bit_count()==8 and w>>21==0],dtype=float)
 A=[]
 for o in octads:
  for j in np.flatnonzero(o):
   a=o.copy();a[j]-=2;A.append(a)
 for i in range(21):
  for j in range(i):
   a=np.zeros(21);a[i]=a[j]=2;A.append(a)
 A=np.array(A);Amax=np.sqrt(2.)
 active=[q for q in range(qsize) if not args.fill_one or q&15 or q==0]
 count=len(active);ii=[];jj=[];sg=[]
 for i,q in enumerate(active):
  for j,r in enumerate(active[:i+1]):
   ids=np.flatnonzero(syndrome==(q^r))
   if q==r:ids=ids[ids!=0]
   ii.extend([i]*len(ids));jj.extend([j]*len(ids));sg.extend(sign[ids])
 ii=np.array(ii);jj=np.array(jj);sg=np.array(sg)
 labels={eta(D[m]):m for m in reversed(line_ids)}
 seed=np.array([[2 if D[labels.get(q&15,line_ids[0])]>>j&1 else 3 for j in range(21)] for q in active],float)
 if args.fill_one:seed[0]=3
 seed/=np.linalg.norm(seed,axis=1)[:,None]
 history=[];best=None;started=time.monotonic()
 for run in range(args.restarts):
  z=seed.copy()
  if run:z+=rng.normal(0,.025+(run%3)*.015,z.shape)
  penalty=1.
  def objective(flat):
   raw=flat.reshape(count,21);norm=np.linalg.norm(raw,axis=1);a=raw/norm[:,None]
   prod=np.einsum('ij,ij,ij->i',a[ii],a[jj],sg)
   residual=np.maximum(prod-.5,0.)
   base=a@A.T;bres=np.maximum(base-Amax,0.)
   # Equal unit scale for the pairwise and base cosine violations.
   loss=float(residual@residual+.5*np.sum(bres*bres))
   grad=np.zeros_like(a)
   np.add.at(grad,ii,2*residual[:,None]*sg*a[jj]);np.add.at(grad,jj,2*residual[:,None]*sg*a[ii])
   grad+=bres@A
   grad=(grad-a*np.sum(a*grad,axis=1)[:,None])/norm[:,None]
   return loss,grad.ravel()
  result=minimize(objective,z.ravel(),jac=True,method='L-BFGS-B',bounds=[(.04,.65)]*z.size,options={'maxiter':args.iterations,'ftol':1e-16,'gtol':1e-10,'maxls':30})
  raw=result.x.reshape(count,21);a=raw/np.linalg.norm(raw,axis=1)[:,None]
  pairmax=float(np.max(np.einsum('ij,ij,ij->i',a[ii],a[jj],sg)));basemax=float(np.max(a@A.T)/np.sqrt(8.))
  row={'restart':run,'loss':float(result.fun),'max_pair_cosine':pairmax,'max_base_cosine':basemax,'iterations':int(result.nit),'success':bool(result.success),'message':str(result.message),'seconds':time.monotonic()-started}
  print(json.dumps(row),flush=True);history.append(row)
  if best is None or max(pairmax,basemax)<best[0]:best=(max(pairmax,basemax),a.copy())
 filename=f'continuous21-q{qsize}-s{args.seed}'+('-fill' if args.fill_one else '')+'.json'
 (OUT/filename).write_text(json.dumps({'seed':args.seed,'extra_masks':extra_masks,'quotient_size':qsize,'active_cosets':active,'points_proposed':count*(4096//qsize),'history':history,'best_shapes':best[1].tolist(),'verification_status':'numerical','new_sota':False,'scope':'Positive amplitudes, fixed odd-octad/root base; complete selected code-coset families all change jointly'})+'\n')
if __name__=='__main__':main()
