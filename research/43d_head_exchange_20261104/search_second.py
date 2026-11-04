"""Search a finite two-for-three exchange on the archived first-stage 43D seed."""
import argparse,base64,collections,hashlib,io,itertools,json,os,sys,tempfile,time
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];P=ROOT/'verifications/improved/dim41-46-p48-heads/partial43';sys.path.insert(0,str(P));import verify_families as vf

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);parser.add_argument('--stage',type=int,choices=(2,3),default=2);args=parser.parse_args();raw=base64.b64decode((HERE/'baseline_2919782.npz.b64').read_bytes());assert hashlib.sha256(raw).hexdigest()=='ff4042cfe8bc6bc79f0b33aa4259df41a08c013e37511e03b6833a18aff321c4';old=dict(np.load(io.BytesIO(raw),allow_pickle=False));delta=json.loads((HERE/'delta.json').read_text());cursor=0
 for step in delta['steps'][:args.stage-1]:
  deleted=step.get('removed_indices',step.get('removed_indices_relative_to_intermediate'));count=step['added_count'];old['body_indices']=np.r_[np.delete(old['body_indices'],deleted),np.asarray(delta['added_body_indices'][cursor:cursor+count],np.int64)];old['head_numerators']=np.vstack((np.delete(old['head_numerators'],deleted,axis=0),np.asarray(delta['head_numerators'][cursor:cursor+count],np.int64)));cursor+=count
 assert len(old['body_indices'])==3447+args.stage-1
 start=time.monotonic();g,res,mu,_=vf.geometry();X=np.load(P/'tiny_bodies.npz')['H'].astype(np.int64);q=np.load(P/'body_complement.npz');tc=vf.mm(vf.mm(X,q['b'].T),vf.inverse(q['q'],3));assert not np.any(tc%288);tc//=288;lookup={row.tobytes():i for i,row in enumerate(res)};labels=np.array([lookup[row.tobytes()] for row in tc%3],np.int32);by={int(c):np.flatnonzero(labels==c) for c in np.unique(labels)}
 model=dict(np.load(P/'configuration_43.npz'));D=int(old['head_denominator']);rawheads=np.unique(old['head_numerators'],axis=0);heads=rawheads.astype(float)/D;oldheads=old['head_numerators'].astype(float)/D;oldnorm=10+(oldheads*oldheads).sum(1);Z=model['head_numerators'].astype(float)/int(model['head_denominator']);rad=model['r3'];fn_float=rad+(Z*Z).sum(1);pure=model['pure_directions'].astype(float);digits=g['digits'];powers=g['powers'];classes=model['classes'];classbeta=np.empty((6561,len(classes)),np.int16)
 for i in range(6561):classbeta[i]=10+rad-18+mu[((digits[classes]-digits[i])%3)@powers]
 beta=classbeta[labels].copy()
 for c,rows in by.items():
  dc=((digits[classes]-digits[c])%3)@powers
  for j in np.flatnonzero(mu[dc]==8):
   pool=by.get(int(dc[j]),np.array([],np.int64));contact=np.any(vf.mm(X[rows],X[pool].T)==(int(rad[j])-20)*144,axis=1) if len(pool) else np.zeros(len(rows),bool);beta[rows[~contact],j]-=6
 existing={(int(b),tuple(map(int,h))) for b,h in zip(old['body_indices'],old['head_numerators'])};groups=collections.defaultdict(list)
 with tempfile.TemporaryDirectory(prefix='k43-second-') as temp:
  dots=np.memmap(Path(temp)/'dots.bin',mode='w+',dtype=np.int8,shape=(len(X),len(old['body_indices'])));xf=X.astype(np.float32);Y=xf[old['body_indices']];assert 48*int(np.max(np.abs(xf)))**2<2**24
  for a in range(0,len(X),1000):
   q=xf[a:a+1000]@Y.T;assert np.all(q==np.rint(q)) and not np.any(np.remainder(q,144));q/=144;assert q.min()>=-128 and q.max()<=127;dots[a:a+1000]=q.astype(np.int8)
  for i,h in enumerate(heads):
   hn=10+h@h
   if not np.all(2*pure@h<=np.sqrt(hn*(pure*pure).sum(1))+1e-10):continue
   ids=np.flatnonzero(np.all(beta+2*(Z@h)[None,:]<=np.sqrt(hn*fn_float)[None,:]+1e-10,axis=1));threshold=np.floor(np.sqrt(hn*oldnorm)-2*(oldheads@h)+1e-9).astype(np.int16);hr=tuple(map(int,rawheads[i]))
   for a in range(0,len(ids),128):
    bs=ids[a:a+128];bad=dots[bs]>threshold[None,:];scores=bad.sum(1)
    for body,row,count in zip(bs,bad,scores):
     if 1<=count<=2 and (int(body),hr) not in existing:groups[tuple(map(int,np.flatnonzero(row)))].append((int(body),int(i)))
  fz=model['head_numerators'].astype(object);fd=int(model['head_denominator']);fr=rad.astype(object);fn=fr*fd*fd+(fz*fz).sum(1);oz=old['head_numerators'].astype(object);on=10*D*D+(oz*oz).sum(1);pu=model['pure_directions'].astype(object);pn=(pu*pu).sum(1)
  def exact(v,removed):
   b,i=v;z=rawheads[i].astype(object);N=10*D*D+sum(z*z);p=beta[b].astype(object)*D*fd+2*fz@z
   if not all(x<=0 or x*x<=N*n for x,n in zip(p,fn)):return False
   p=dots[b].astype(object)*D*D+2*oz@z
   if not all(j in removed or x<=0 or x*x<=N*n for j,(x,n) in enumerate(zip(p,on))):return False
   return all(x<=0 or x*x<=N*n for x,n in zip(2*pu@z,pn))
  def pair(v,w):
   b,i=v;c,j=w;x=int(X[b]@X[c]);assert x%144==0;z=rawheads[i].astype(object);q=rawheads[j].astype(object);N=10*D*D+sum(z*z);M=10*D*D+sum(q*q);p=(x//144)*D*D+2*sum(z*q);return p<=0 or p*p<=N*M
  blockers={k for k in groups if len(k)==2};tested=0
  for j,k in blockers:
   V=groups.get((j,),[])+groups.get((k,),[])+groups.get((j,k),[])
   if len(V)<3:continue
   ids=np.array([v[0] for v in V]);hs=np.array([heads[v[1]] for v in V]);d=X[ids]@X[ids].T;assert not np.any(d%144);d//=144;n=10+(hs*hs).sum(1);compatible=d<=np.sqrt(n[:,None]*n[None,:])-2*(hs@hs.T)+1e-9;np.fill_diagonal(compatible,False)
   for i in range(len(V)):
    for jj in np.flatnonzero(compatible[i,i+1:])+i+1:
     for kk in np.flatnonzero(compatible[i,jj+1:]&compatible[jj,jj+1:])+jj+1:
      tested+=1;tri=[V[i],V[jj],V[kk]]
      if all(pair(x,y) for x,y in itertools.combinations(tri,2)) and all(exact(x,{j,k}) for x in tri):
       out={'remove':[j,k],'add':[[b,i] for b,i in tri],'head_denominator':D,'head_numerators':[rawheads[i].tolist() for b,i in tri]}
       if args.output:args.output.write_text(json.dumps(out,indent=2)+'\n')
       print(json.dumps({'status':'two-for-three proposal, exact local checks pass','tested_triangles':tested,'elapsed_seconds':time.monotonic()-start,'delta':out},indent=2));return
 print(json.dumps({'status':'no proposal found in the screened two-blocker pool','tested_triangles':tested}))
if __name__=='__main__':main()
