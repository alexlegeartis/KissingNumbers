"""Exact gate for exceptional E6/E7 replacement of a full tail-axis layer."""
from pathlib import Path
import argparse,json,time,itertools,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
p=argparse.ArgumentParser();p.add_argument('dimension',type=int,choices=(30,31));p.add_argument('--axis',type=int,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();st=time.monotonic()
source=next((ROOT/('verifications/improved/dim27-31-'+('algebraic-completion' if a.dimension==31 else 'joint-refills')+'/data')).glob(f'coordinates{a.dimension}_*.npz'))
X=np.load(source)['points'].astype(np.int16);d=X.shape[1];assert 24<=a.axis<d
T=np.sum(X[:,24:].astype(np.int64)**2,axis=1); B=np.unique(X[T==32,:24],axis=0);old=X[X[:,a.axis]==0];aux=old[np.any(old[:,24:]!=0,axis=1)]
other=[j for j in range(24,d) if j!=a.axis];tails=[];rad=[]
for s in itertools.product((-1,1),repeat=len(other)):
 if d==31 and np.prod(s)!=1:continue
 for r in ([int(np.prod(s))] if d==30 else (-1,1)):
  t=np.zeros(d-24,np.int16);t[np.array(other)-24]=2*np.array(s);tails.append(t);rad.append(r)
tails=np.array(tails);rad=np.array(rad);m=len(tails)
# Candidate coordinate on the erased axis is 2 sqrt(3) * rad in E6,
# and 2 sqrt(2) * rad in E7. Retained points vanish on that axis,
# so ALL comparisons in this gate are exact integers.
bi=np.repeat(np.arange(len(B)),m);ti=np.tile(np.arange(m),len(B));Y=np.zeros((len(bi),d),np.int16);Y[:,:24]=B[bi];Y[:,24:]=tails[ti]
assert np.all(np.sum(Y.astype(np.int64)**2,axis=1)+(12 if d==30 else 8)==128)
print(json.dumps({'dimension':d,'axis':a.axis,'removed':len(X)-len(old),'bodies':len(B),'spinor_directions':m,'raw_candidates':len(Y),'fixed_aux':len(aux)}),flush=True)
rng=np.random.default_rng(20261108);aux=aux[rng.permutation(len(aux))];active=np.arange(len(Y));checks=0;history=[]
for start in range(0,len(aux),256):
 if not len(active):break
 P=Y[active].astype(np.float32)@aux[start:start+256].astype(np.float32).T;checks+=P.size
 active=active[np.max(P,axis=1)<=64]
 if start%8192==0:print(start,len(active),flush=True);history.append([start,int(len(active))])
# Core comparisons are inherited exactly: every body is the core part of an
# existing cap and the complete original core is retained unchanged here.
# Count is an upper bound before new/new compatibility or other new shells.
a.out.mkdir(parents=True,exist_ok=False)
np.savez_compressed(a.out/'admitted.npz',rational=Y[active],radical_sign=rad[ti[active]],radical_axis=a.axis,radicand=3 if d==30 else 2,retained=old)
report={'source':str(source.relative_to(ROOT)),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'dimension':d,'axis':a.axis,'baseline':len(X),'retained':len(old),'deleted':len(X)-len(old),'unique_bodies':len(B),'spinor_directions':m,'candidate_count':len(Y),'admitted_count':len(active),'potential_total_before_new_new':len(old)+len(active),'exact_products':checks,'seconds':time.monotonic()-st,'history':history,'scope':'All norm96 bodies in incumbent paired with every exceptional spinor tail. Fixed-point gate is exact; no conclusion beyond this family. New/new constraints not yet applied.'}
(a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
