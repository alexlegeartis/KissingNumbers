#!/usr/bin/env python3
"""Rebuild the pinned parent positions and append compatible horizontal groups."""
import os,json,runpy,sys,tempfile,subprocess
from pathlib import Path
from hashlib import sha256
import numpy as np
from refill import Model,sizes
from components import improve
assert __debug__
here=Path(__file__).resolve().parent
if len(sys.argv)==1:
 for dim in ('84','85'):
  subprocess.run([sys.executable,str(Path(__file__).resolve()),dim],check=True)
 raise SystemExit(0)
d=sys.argv[1];assert d in ('84','85')
pins={
 '84':('dim81-84-antipodal-tail-codes',{
 'verify.py':'cc5c2f1ffe83d006e1a7737f6c1b3e71a34100f84bb1c71419bd0164eda1a906',
 'metric_profiles.py':'65f50490e00ee84fb1a699663a766c2d29905dc232e46cf99e07f25107a2eba5',
 'claims.json':'0e31ef9dc03a782d4aaa97788f20aa788652b0e8c8fa8fc0ff6f374086b43c85',
 'phases.json':'cda345b9f3534ac87a3f56240095fc7b16ccf96ec67dcbe3ec87865f1ab41f11',
 'P12a_tails.json':'be558e887459b65345b1981f85852064f0a71fcae38f578cfa85bc75231d2171'}),
 '85':('dim82-85-bw16-tail-sections',{
 'verify.py':'28dc2b85d228045e6974e5ede2d940ff09c2beebf23be8ac99a18e3bdbaa7975',
 'phases.json':'050157ff3f8ee6282e4092dd208bbd8ca1f83ef4b4a6efdda9bdc16735b1f21f',
 'claims.json':'a7ecc17a84cea1ca130e50b75f5273a5ac5ac95eb788569e0ac05b394eba5f1a'}),
}
base,checks=pins[d];base=here.parent/base
assert all(sha256((base/n).read_bytes()).hexdigest()==h for n,h in checks.items())
if d=='85':
 helper=here.parent/'dim81-85-gamma72-half-sign/two_profiles.py'
 assert sha256(helper.read_bytes()).hexdigest()=='0159dc533e090920cc3371b4744437e9371d5868ecbe3dcab103e29a88ca70d4'
with tempfile.TemporaryDirectory() as aux:
 os.environ['EXPORT_PROFILE']=aux;os.environ['EXPORT_TRIPLE']=aux
 sys.argv=[str(base/'verify.py')]
 R=runpy.run_path(sys.argv[0])
 X={k:v for k,v in np.load(Path(aux)/'profile_stock.npz').items()}
 O={k:v for k,v in np.load(Path(aux)/f'profile_{d}.npz').items()}
 T=np.load(Path(aux)/'triple384.npy').astype(bool)
 del os.environ['EXPORT_PROFILE'],os.environ['EXPORT_TRIPLE']
 H=list(map(int,X['H']));R0=R['R']
 Ainv=np.array(sorted(R0['RA'],key=R0['RA'].get),dtype='u4')
 Binv=np.array(sorted(R0['RB'],key=R0['RB'].get),dtype='u4')
 Q=[]
 for t in H:
  A,B=R0['T'][R0['RB'][t]]
  Q.append(np.stack((Ainv[A[:,2]],Binv[B[:,0]],Binv[B[:,1]]),axis=1))
 M=Model(*[X[k]for k in ('H','A8','A6','F','V8','V6','U','V')],np.array(Q))
 if d=='84':
  V=np.array(json.loads((base/'P12a_tails.json').read_text()),dtype='i8')
  rows=[tuple(r)for r in R['phases'][d]['assign']]
  assert R['out'][d]==6291169688
  assert len(V)==420 and np.all(np.diag(V@V.T)==4)
  assert int(T.sum())*16+int(O['end'])*8+840+6218175600==R['out'][d]
  original=R['out'][d]; incumbent=original
  earlier=[]
 else:
  V=R['tails'](13)
  rows=[tuple(r)for r in R['phases'][d]]
  assert R['out'][d]==6292777568
  assert len(V)==453 and np.all(np.diag(V@V.T)==8)
  assert int(T.sum())==4028344
  assert int(T.sum())*16+int(O['end'])*8+768+6218175600==R['out'][d]
  original=R['out'][d];incumbent=6294162634
  e=here.parent/'dim85-shared-phase-refill'/'extra.json'
  assert sha256(e.read_bytes()).hexdigest()=='dcbf8cacdec2c3d6e2e28774a3499cd27d7fdcc92ec5e12e9871df557904c210'
  earlier=json.loads(e.read_text()); assert len(earlier)==69
 G=V@V.T; n=len(G);assert np.max(np.abs(G-np.diag(np.diag(G))))<=G[0,0]//2
 assert len(rows)==384 and len({(c,t)for _,c,t in rows})==384
 assert np.array_equal(O['colors'],[r[1]+1 for r in rows])
 assert np.array_equal(O['phase'],[r[2]for r in rows])
 occ=list(zip(T,O['occ48'],O['occ66']))
 added=json.loads((here/f'append_{d}.json').read_text())
 assert (d=='84' and len(added)==36) or (d=='85' and [a[:4] for a in added]==earlier)
 assert len(set(range(n))-{i for i,c,t in rows})==len(added)
 extra=[0,0,0]
 for i,c,t,a,b,e in added:
  assert 0<=i<n and i not in {j for j,_,_ in rows}
  assert (c,t) in {(y,z)for _,y,z in rows[:384]}
  f=M.add(i,c,t,rows,occ,G)
  assert sizes(f)==(a,b,e) and a>0,(i,sizes(f),(a,b,e))
  extra=[u+v for u,v in zip(extra,(a,b,e))]
  rows.append((i,c,t));occ.append(f)
 pure=2*len(added) if d=='85' else 0
 count=original+16*extra[0]+8*(extra[1]+extra[2])+pure
 g=improve(M,np.array(Q),G,rows,occ)
 assert g=={"84":26764,"85":58720}[d]
 count+=16*g
 claims=json.loads((here/'claims.json').read_text())['lower_bounds']
 assert count==claims[d] and count>incumbent
 print('PASS appended profiles',d,count,'additional',count-incumbent,'groups',extra,'pure',pure,'reallocated triples',g)
