"""Deterministic finite 43D one-for-two discovery; exact checks gate output."""
import argparse,collections,hashlib,itertools,json,os,sys,tempfile,time
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];P=ROOT/'verifications/improved/dim41-46-p48-heads/partial43';sys.path.insert(0,str(P));import verify_families as vf

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args();base=P/'explicit_points_43.npz';assert hashlib.sha256(base.read_bytes()).hexdigest()=='ff4042cfe8bc6bc79f0b33aa4259df41a08c013e37511e03b6833a18aff321c4';start=time.monotonic();g,res,mu,_=vf.geometry();X=np.load(P/'tiny_bodies.npz')['H'].astype(np.int64);q=np.load(P/'body_complement.npz');tc=vf.mm(vf.mm(X,q['b'].T),vf.inverse(q['q'],3));assert not np.any(tc%288);tc//=288;lookup={row.tobytes():i for i,row in enumerate(res)};labels=np.array([lookup[row.tobytes()] for row in tc%3],np.int32)
 model=dict(np.load(P/'configuration_43.npz'));old=dict(np.load(base));D=int(old['head_denominator']);raw=np.unique(old['head_numerators'],axis=0);heads=raw.astype(float)/D;oldhead=old['head_numerators'].astype(float)/D;oldnorm=10+(oldhead*oldhead).sum(1);Z=model['head_numerators'].astype(float)/int(model['head_denominator']);rad=model['r3'];fn_float=rad+(Z*Z).sum(1);pure=model['pure_directions'].astype(float);digits=g['digits'];powers=g['powers'];classes=model['classes'];beta=np.empty((6561,len(classes)),np.int16)
 for i in range(6561):beta[i]=10+rad-18+mu[((digits[classes]-digits[i])%3)@powers]
 existing={(int(b),tuple(map(int,h))) for b,h in zip(old['body_indices'],old['head_numerators'])};groups=collections.defaultdict(list)
 with tempfile.TemporaryDirectory(prefix='k43-search-') as temp:
  dots=np.memmap(Path(temp)/'dots.bin',mode='w+',dtype=np.int8,shape=(len(X),len(old['body_indices'])));xf=X.astype(np.float32);Y=xf[old['body_indices']];assert 48*int(np.max(np.abs(xf)))**2<2**24
  for a in range(0,len(X),1000):
   q=xf[a:a+1000]@Y.T;assert np.all(q==np.rint(q)) and not np.any(np.remainder(q,144));q/=144;assert q.min()>=-128 and q.max()<=127;dots[a:a+1000]=q.astype(np.int8)
  for i,h in enumerate(heads):
   hn=10+h@h
   if not np.all(2*pure@h<=np.sqrt(hn*(pure*pure).sum(1))+1e-10):continue
   allowed=np.all(beta+2*(Z@h)[None,:]<=np.sqrt(hn*fn_float)[None,:]+1e-10,axis=1);ids=np.flatnonzero(allowed[labels]);threshold=np.floor(np.sqrt(hn*oldnorm)-2*(oldhead@h)+1e-9).astype(np.int16);hr=tuple(map(int,raw[i]))
   for a in range(0,len(ids),128):
    bs=ids[a:a+128];bad=dots[bs]>threshold[None,:];scores=bad.sum(1)
    for body,row,count in zip(bs,bad,scores):
     if count==1 and (int(body),hr) not in existing:groups[int(np.flatnonzero(row)[0])].append((int(body),int(i)))
  fz=model['head_numerators'].astype(object);fd=int(model['head_denominator']);fr=model['r3'].astype(object);fn=fr*fd*fd+(fz*fz).sum(1);oz=old['head_numerators'].astype(object);on=10*D*D+(oz*oz).sum(1);pu=model['pure_directions'].astype(object);pn=(pu*pu).sum(1)
  def exact(v,remove):
   b,i=v;z=raw[i].astype(object);N=10*D*D+sum(z*z);p=beta[labels[b]].astype(object)*D*fd+2*fz@z
   if not all(x<=0 or x*x<=N*n for x,n in zip(p,fn)):return False
   p=dots[b].astype(object)*D*D+2*oz@z
   if not all(j==remove or x<=0 or x*x<=N*n for j,(x,n) in enumerate(zip(p,on))):return False
   return all(x<=0 or x*x<=N*n for x,n in zip(2*pu@z,pn))
  def pair(v,w):
   b,i=v;c,j=w;x=int(X[b]@X[c]);assert x%144==0;z=raw[i].astype(object);q=raw[j].astype(object);N=10*D*D+sum(z*z);M=10*D*D+sum(q*q);p=(x//144)*D*D+2*sum(z*q);return p<=0 or p*p<=N*M
  tested=0
  for blocker,vertices in groups.items():
   for a,b in itertools.combinations(vertices,2):
    tested+=1
    if pair(a,b) and exact(a,blocker) and exact(b,blocker):
     result={'remove':blocker,'add':[[a[0],a[1]],[b[0],b[1]]],'head_denominator':D,'head_numerators':[raw[a[1]].tolist(),raw[b[1]].tolist()]}
     if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
     print(json.dumps({'status':'one-for-two proposal, exact local checks pass; run replay.py on the shipped certificate','tested_pairs':tested,'elapsed_seconds':time.monotonic()-start,'delta':result},indent=2));return
 print(json.dumps({'status':'no proposal found in the finite screened pool','tested_pairs':tested}))
if __name__=='__main__':main()
