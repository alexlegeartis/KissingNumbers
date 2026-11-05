"""Materialize supplied 12–16D benchmark seeds from tracked recipes."""
import json,pathlib,sys,numpy as np
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]))
from discovery.core import RunContext
from discovery.engines.partial_group_assemblies import _group
from discovery.engines.group_seed_search import _closure
root=pathlib.Path('research/generic_reflection_campaign_20261105/cache');root.mkdir(parents=True,exist_ok=True);ctx=RunContext(pathlib.Path('.'),{'max_orbit_points':20000,'max_gram_entries':10000000,'max_nodes':3000,'max_rank_work':1000000000})
def params(name):return json.load(open('discovery/recipes/'+name))['nodes'][0]['params']
p=params('group13_quadratic_benchmark.json');A=np.array(p['initial_a'],np.int64);B=np.array(p['initial_b'],np.int64);assert A.shape==B.shape==(1154,13) and all(not s['generators'] for s in p['orbit_specs']);np.savez_compressed(root/'seed13.npz',A=A,B=B)
mask=(B==0).all(1)&(A[:,12]==0)&(np.count_nonzero(A[:,:12],axis=1)==4)&np.isin(A[:,:12],[-2,0,2]).all(1);X=np.r_[A[mask,:12],4*np.eye(12,dtype=np.int64),-4*np.eye(12,dtype=np.int64)];assert X.shape==(840,12);np.savez_compressed(root/'seed12.npz',X=X,g=np.ones(12,np.int64))
p=params('group14_mixed_stabilizer_benchmark.json');G,den=_closure(p['generators'],14,ctx);parts=[]
for x in p['initial_seeds']:parts.append(np.unique(G@np.array(x,np.int64),axis=0))
X=np.concatenate(parts);assert len(X)==1932 and len(np.unique(X,axis=0))==len(X);np.savez_compressed(root/'seed14.npz',X=X,g=np.ones(14,np.int64))
for d in (15,16):
 p=params(f'group{d}_'+('weighted_partial' if d==15 else 'partial')+'_benchmark.json');parts=[]
 for x,s in zip(p['initial_seeds'],p['orbit_specs']):
  G,den=_group(s.get('generators',[]),d,ctx);parts.append(np.unique(G@np.array(x,np.int64),axis=0))
 X=np.concatenate(parts);assert len(X)==({15:2564,16:4320}[d]) and len(np.unique(X,axis=0))==len(X);np.savez_compressed(root/f'seed{d}.npz',X=X,g=np.array(p.get('metric',[1]*d),np.int64))
print('materialized 12–16D tracked seeds; search operators independently validate them')
