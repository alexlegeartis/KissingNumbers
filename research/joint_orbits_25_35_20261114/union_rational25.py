"""Try exact unions with rational points and restored shell from older 25D codes."""
import argparse,collections,importlib.util,json,math,time
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
ROOT=Path(__file__).resolve().parents[2]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(exist_ok=False);t0=time.monotonic()
 src=ROOT/'verifications/improved/dim25-critical-cluster/source';spec=importlib.util.spec_from_file_location('geometry',src/'verify_general.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 cert=json.loads((src/'certificate25.json').read_text());shell,_=m.independent_shell();ss=set(shell);co,ds,removed,*_=m.reconstruct(cert,ss);bulk=np.array([v for v in shell if v not in removed],np.int64)
 iv=[m.enclosure(x,d) for x,d in zip(co,ds)];L=np.array([v[0] for v in iv]);U=np.array([v[1] for v in iv]);err=np.max(U-L,axis=1);l1=np.sum(abs(L),axis=1);bl1=np.sum(abs(bulk),axis=1)
 def key(v,d):
  g=math.gcd(d,math.gcd(*map(int,v.ravel())));return tuple(int(x)//g for x in v.ravel()),d//g
 oldkeys={key(v,d) for v,d in zip(co,ds)};reports=[]
 for donor in ['dim25-point-split','dim25-collective-deformation']:
  f=ROOT/f'verifications/improved/{donor}/source/certificate25.json'
  if not f.exists():continue
  dc=json.loads(f.read_text());_,_,dr,*_=m.reconstruct(dc,ss);vs=[];den=[]
  for x in dc.get('small',[]):
   v=np.zeros((4,25),object);v[0]=list(map(int,x['numerator']));d=int(x['denominator'])
   if key(v,d) not in oldkeys:vs.append(v);den.append(d)
  for u in removed-dr:
   v=np.zeros((4,25),object);v[0,:24]=u
   if key(v,1) not in oldkeys:vs.append(v);den.append(1)
  rows=[];cols=[];fallbacks=0
  for i,(v,d) in enumerate(zip(vs,den)):
   assert m.field_dot(v,v)==[32*d*d,0,0,0]
   lo,hi=m.enclosure(v,d);e=int(np.max(hi-lo));base=bulk@lo[:24];possible=np.flatnonzero(base+bl1*e>16*m.SCALE)
   for j in possible:
    fallbacks+=1
    if m.dot(bulk[j],v[0,:24])>16*d:rows.append(i);cols.append(int(j))
   b=L@lo+err*int(np.sum(abs(lo)))+l1*e+25*e*err
   for j in np.flatnonzero(b>16*m.SCALE*m.SCALE):
    fallbacks+=1;val=[m.dot(v[0],c) for c in co[j]];val[0]-=16*d*ds[j]
    if m.field_sign(val)>0:rows.append(i);cols.append(len(bulk)+int(j))
  G=csr_matrix((np.ones(len(rows),np.int8),(rows,cols)),shape=(len(vs),len(bulk)+len(co)));match=maximum_bipartite_matching(G,perm_type='column');rev=np.full(G.shape[1],-1,int)
  for i,j in enumerate(match):
   if j>=0:rev[j]=i
  ci=match<0;oi=np.zeros(G.shape[1],bool);q=collections.deque(np.flatnonzero(ci))
  while q:
   i=q.popleft()
   for j in G.indices[G.indptr[i]:G.indptr[i+1]]:
    if j==match[i] or oi[j]:continue
    oi[j]=True;assert rev[j]>=0
    if not ci[rev[j]]:ci[rev[j]]=True;q.append(int(rev[j]))
  gain=int(ci.sum()-oi.sum());r=dict(donor=donor,new=len(vs),edges=G.nnz,added=int(ci.sum()),removed=int(oi.sum()),gain=gain,exact_fallbacks=fallbacks);reports.append(r);print('EXACT_UNION25',json.dumps(r),flush=True)
  if gain>0:
   proposal=dict(source=str(src/'certificate25.json'),donor=str(f),source_bulk_count=len(bulk),removed_indices=np.flatnonzero(oi).tolist(),added=[dict(numerator=list(map(str,vs[i][0])),denominator=str(den[i])) for i in np.flatnonzero(ci)],count=197586+gain)
   (a.out/f'proposal25_{197586+gain}.json').write_text(json.dumps(proposal,indent=2)+'\n')
 (a.out/'report.json').write_text(json.dumps(dict(reports=reports,seconds=time.monotonic()-t0),indent=2)+'\n')
if __name__=='__main__':main()
