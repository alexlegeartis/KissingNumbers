"""Bounded continuous head search; numerical output is never a certificate."""
from pathlib import Path
from time import monotonic
import json,numpy as np
from scipy.optimize import minimize
from head_system import build,bound

p=Path(__file__).parent;S=build(4);H=np.linalg.inv(S['gram']);cells=S['cells'];qall=np.einsum('ij,jk,ik->i',cells,H,cells)
C=json.loads((p/'twoheads_k4_refined_result.json').read_text())
core=np.array([x['profile'] for x in C['slices']],dtype=int);head0=np.array([x['head'] for x in C['slices']]);q0=np.einsum('ij,jk,ik->i',core,H,core)
choices=[]
grid=np.array([(r*np.cos(t),r*np.sin(t)) for r in [0,.125,.25,.375,.5] for t in np.arange(24)*np.pi/12])
for c in cells[np.abs(qall-.25)<1e-8]:
 if any(np.array_equal(c,a) for a in core):continue
 P=4-core@H@c
 a0=np.sqrt((8-np.sum(head0*head0,axis=1))/(8-q0))
 trial=4-a0[:,None]*np.sqrt((8-np.sum(grid*grid,axis=1))/7.75)[None,:]*P[:,None]-head0@grid.T
 ii=np.argmax(np.min(trial,axis=0));g=trial[:,ii]
 choices.append((float(min(g)),c.copy(),grid[ii].copy(),g))
choices.sort(key=lambda x:-x[0]);start=monotonic();attempts=[];best=-10;result=None
for initial,c,newhead,oldgaps in choices:
 for mode in [0,1,2]:
  if monotonic()-start>170:break
  keep=np.ones(len(core),dtype=bool)
  if mode:
   removable=np.flatnonzero(q0>.99)
   bad=sorted(removable,key=lambda i:oldgaps[i])[:mode]
   keep[bad]=False
  profiles=np.vstack([core[keep],c]);init=np.vstack([head0[keep],newhead]);m=len(profiles)
  q=np.einsum('ij,jk,ik->i',profiles,H,profiles);P=4-profiles@H@profiles.T;I,J=np.triu_indices(m,1);nc=len(I)
  def values(z):
   h=z[:-1].reshape(m,2);a=np.sqrt(np.maximum(1e-10,8-np.sum(h*h,axis=1))/(8-q))
   pair=4-a[I]*a[J]*P[I,J]-np.sum(h[I]*h[J],axis=1)-z[-1]
   return np.r_[pair,q-np.sum(h*h,axis=1)]
  def jac(z):
   h=z[:-1].reshape(m,2);a=np.sqrt(np.maximum(1e-10,8-np.sum(h*h,axis=1))/(8-q));K=np.zeros((nc+m,2*m+1))
   di=(a[J]*P[I,J]/((8-q[I])*a[I]))[:,None]*h[I]-h[J]
   dj=(a[I]*P[I,J]/((8-q[J])*a[J]))[:,None]*h[J]-h[I]
   for l in range(2):K[np.arange(nc),2*I+l]=di[:,l];K[np.arange(nc),2*J+l]=dj[:,l];K[nc+np.arange(m),2*np.arange(m)+l]=-2*h[:,l]
   K[:nc,-1]=-1;return K
  x=np.r_[init.ravel(),-.2];x[-1]=min(values(x)[:nc])-.2
  grad=np.zeros(2*m+1);grad[-1]=-1
  opt=minimize(lambda z:-z[-1],x,jac=lambda z:grad,method='SLSQP',bounds=[(-np.sqrt(q[i//2]),np.sqrt(q[i//2])) for i in range(2*m)]+[(-2,.05)],constraints={'type':'ineq','fun':values,'jac':jac},options={'maxiter':180,'ftol':1e-10,'disp':False})
  h=opt.x[:-1].reshape(m,2);z=opt.x.copy();z[-1]=0;g=float(np.min(values(z)));margin=float(opt.x[-1]);ids=[next(i for i,d in enumerate(cells) if np.array_equal(a,d)) for a in profiles]
  pop=bound(S,ids);rec={'added':c.tolist(),'drop_mode':mode,'pairs_min':g,'margin_variable':margin,'population_lower_numeric':pop,'optimizer_success':bool(opt.success),'iterations':int(opt.nit),'seconds':monotonic()-start};attempts.append(rec);print(rec,flush=True)
  if g>best:best=g
  if g>=-1e-8 and pop>1360986112:
   result=dict(rec,profiles=profiles.tolist(),heads=h.tolist(),status='CANDIDATE ONLY; exact certification mandatory');break
 if result is not None or monotonic()-start>170:break
report={'attempts':attempts,'candidate':result,'best_geometry_margin':best,'seconds':monotonic()-start,'scope':'Bounded numerical attempts; failure is not an obstruction theorem.'}
(p/'continuous70_result.json').write_text(json.dumps(report,indent=2));print('FINAL',report['candidate'],report['seconds'],flush=True)
