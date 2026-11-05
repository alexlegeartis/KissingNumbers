"""Independent connected-component search for one-point exchanges."""
import time
import numpy as np
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix,vstack
from scipy.sparse.csgraph import connected_components


def search_components(m,n,source_i,source_j,conflict_i,conflict_j,seconds):
 """Return first positive exchange; absence of one is a numerical result only.

 Candidate-source edges connect a proposed point to an old point it excludes.
 Candidate-candidate edges connect mutually incompatible new points. Since the
 gain is additive, a feasible global gain implies a positive-gain component.
 """
 u=np.r_[source_i,conflict_i];v=np.r_[m+source_j,conflict_j]
 graph=coo_matrix((np.ones(2*len(u),np.int8),(np.r_[u,v],np.r_[v,u])),shape=(m+n,m+n)).tocsr()
 count,labels=connected_components(graph,directed=False);sizes=np.bincount(labels);order=sorted(set(labels[:m]),key=lambda k:sizes[k]);start=time.monotonic();reports=[];searched=0
 for k in order:
  left=seconds-(time.monotonic()-start)
  if left<=0:reports.append('total time limit');break
  ci=np.flatnonzero(labels[:m]==k);sj=np.flatnonzero(labels[m:]==k);lc=np.full(m,-1,np.int64);ls=np.full(n,-1,np.int64);lc[ci]=np.arange(len(ci));ls[sj]=np.arange(len(sj));a=(lc[source_i]>=0);b=(lc[conflict_i]>=0);x=lc[source_i[a]];y=ls[source_j[a]];i=lc[conflict_i[b]];j=lc[conflict_j[b]];ne=len(x);e=len(i)
  rows=np.r_[np.repeat(np.arange(ne),2),np.repeat(np.arange(ne,ne+e),2)];cols=np.r_[np.c_[x,len(ci)+y].ravel(),np.c_[i,j].ravel()];vals=np.r_[np.tile([1.,-1.],ne),np.ones(2*e)];mat=coo_matrix((vals,(rows,cols)),shape=(ne+e,len(ci)+len(sj))).tocsr();obj=np.r_[-np.ones(len(ci)),np.ones(len(sj))];mat=vstack([mat,coo_matrix(obj.reshape(1,-1))],format='csr');up=np.r_[np.zeros(ne),np.ones(e),-1.]
  searched+=1
  r=milp(obj,integrality=np.ones(len(obj)),bounds=Bounds(0,1),constraints=LinearConstraint(mat,-np.inf,up),options={'time_limit':float(left),'presolve':True});reports.append(str(r.message))
  if r.x is not None:
   pick=ci[r.x[:len(ci)]>.5].tolist();rm=sj[r.x[len(ci):]>.5].tolist()
   if len(pick)>len(rm):return pick,rm,'component proposal; '+str(r.message)
 return [],[],f'no proposal; {searched}/{len(order)} candidate components searched; last: '+(reports[-1] if reports else 'none')
