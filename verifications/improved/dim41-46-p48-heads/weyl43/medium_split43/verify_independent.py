from pathlib import Path
import numpy as np, sympy as sp, json, hashlib, sys
import argparse,tempfile,subprocess
P=Path(__file__).resolve().parent;W=P.parent;E=P/'medium_families';CAND=P/'candidate.json'
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=P/'independent_replay');args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
with tempfile.TemporaryDirectory(prefix='kissing43_base_') as directory:
 subprocess.run([sys.executable,str(W/'independent_verify.py'),'--output',directory],check=True,capture_output=True,text=True)
sys.path.insert(0,str(W)); import verify_families as vf
g,residues,mu,tiny=vf.geometry()
f=np.load(W/'configuration_43.npz'); t=np.load(W/'explicit_points_43.npz'); bc=np.load(W/'body_complement.npz'); B=bc['b'];Q=bc['q'];J=np.array((3*sp.Matrix(Q.tolist()).inv()).tolist(),dtype=np.int64)
assert np.array_equal(Q@J,3*np.eye(40,dtype=np.int64)); assert np.linalg.matrix_rank(B)==40
c=f['classes'];r=f['r3'];w=f['weights'];z=f['head_numerators'].astype(object);d=int(f['head_denominator']);u=f['pure_directions'].astype(object);I=t['body_indices'];tz=t['head_numerators'].astype(object);td=int(t['head_denominator']);tinyH=np.load(W/'tiny_bodies.npz')['H'].astype(np.int64)
assert int(f['dimension'])==43 and z.shape[1]==u.shape[1]==3
medium_ids=np.flatnonzero(r==12);assert len(medium_ids)==60 and len(set(c[r==12]))==60
allH=[]
for k in medium_ids:
 data=np.load(E/f'family_{c[k]}.npz');H=data['H'].astype(np.int64)
 assert H.shape==(495,48) and len(np.unique(H,axis=0))==495
 assert np.all(np.sum(H*H,axis=1)==3456)
 coeffnum=H@B.T@J;assert np.all(coeffnum%288==0);coeff=coeffnum//288
 assert np.array_equal(coeff@B,H) and np.all(coeff%3==residues[c[k]])
 assert np.array_equal(coeff,data['coefficients']) and int(data['coset'])==c[k]
 allH.append(H)
MH=np.concatenate(allH);assert len(np.unique(MH,axis=0))==29700
spec=json.loads(CAND.read_text());D=spec['head_denominator'];assert isinstance(D,int) and D>0
assert len(spec['delete'])==len(set(spec['delete']))==6
assert len(spec['delete_tiny'])==len(set(spec['delete_tiny']))==3
assert all(0<=i<len(I) for i in spec['delete_tiny'])
assert len(spec['points'])==12 and len(spec['new_tiny'])==1 and spec['gain']==4
assert set(spec['delete'])==set(p['family']*495+p['row'] for p in spec['points'])
assert all(sum(p['family']*495+p['row']==i for p in spec['points'])==2 for i in spec['delete'])
kt=~np.isin(np.arange(len(I)),spec['delete_tiny'])
lookup={tuple(v):i for i,v in enumerate(residues)}
point_ids=np.concatenate([k*495+np.arange(495) for k in medium_ids]);keep=~np.isin(point_ids,spec['delete']);assert np.count_nonzero(~keep)==6;heads=np.repeat(z[medium_ids],495,axis=0);norms=12*d*d+np.sum(heads*heads,axis=1)
counts={}; allnew=[]; slackmin={}
def check(name,v,nn):
 vals=np.asarray(v,dtype=object).reshape(-1);ns=np.broadcast_to(nn,np.shape(v)).reshape(-1)
 assert all(int(x)<=0 or int(x)**2<=int(n) for x,n in zip(vals,ns)),name
 counts[name]=counts.get(name,0)+len(vals)
 positive=[int(n)-int(x)**2 for x,n in zip(vals,ns) if x>0]
 if positive:slackmin[name]=min(slackmin.get(name,min(positive)),min(positive))
for p in spec['points']+spec['new_tiny']:
 if 'family' in p:
  assert p['family'] in medium_ids and 0<=p['row']<495
  h=MH[np.flatnonzero(point_ids==p['family']*495+p['row'])[0]]; rr=12;cl=int(c[p['family']])
 else:
  assert 0<=p['body_index']<len(tinyH) and p['head_denominator']==D
  h=tinyH[p['body_index']];rr=10
  num=h@B.T@J;assert np.all(num%288==0);co=num//288
  assert np.array_equal(co@B,h);cl=lookup[tuple(co%3)]
 assert int(h@h)==rr*288
 a=np.array(p['head_numerators'],dtype=object);assert a.shape==(3,)
 n=rr*D*D+sum(a*a);assert n==p['normal']
 # Physical vectors: (h/sqrt(864), a/(D sqrt(3))). Scale squared norms by 3 D^2.
 dot=MH@h; assert np.all(dot%144==0)
 check('retained_medium',(dot[keep]//144).astype(object)*D*d+2*heads[keep]@a,n*norms[keep])
 dot=tinyH[I[kt]]@h;assert np.all(dot%144==0)
 check('retained_tiny',(dot//144).astype(object)*D*td+2*tz[kt]@a,n*(10*td*td+np.sum(tz[kt]*tz[kt],axis=1)))
 check('pure',2*u@a,n*np.sum(u*u,axis=1))
 for k in np.flatnonzero(r!=12):
  difference=int(((g['digits'][c[k]]-g['digits'][cl])%3)@g['powers'])
  # Inherited minimum-six parent bound: 2<body,body'> <= (r+r'-18+mu)/3.
  beta=rr+int(r[k])-18+int(mu[difference]); v=beta*D*d+2*sum(a*z[k]);N=int(r[k])*d*d+sum(z[k]*z[k])
  check('major_family_bound',v,n*N)
 for hh,aa,nn in allnew:
  assert int(h@hh)%144==0;v=int(h@hh)//144*D*D+2*sum(a*aa); check('new_pair',v,n*nn)
 allnew.append((h,a,n))
base=int(sum(w))+len(I)+len(u);assert base==2919799
assert len(allnew)-sum(~keep)-sum(~kt)==4
report={'status':'PASS','dimension':43,'base_count':base,'medium_enumerated_distinct':len(MH),'deleted_medium':6,'deleted_tiny':3,'inserted':13,'new_count':base+4,'total':base+4,'checks':counts,'total_checks':sum(counts.values()),'positive_dot_min_integer_slack':slackmin,'candidate_sha256':hashlib.sha256(CAND.read_bytes()).hexdigest(),'scope':'Published parent geometry and old configuration replayed separately. All 60 medium families independently checked against published complement and gluing labels. New interactions exact; retained-retained inherited from freshly rerun published verifier. Minimum-six/11-design and family population/completeness identities remain inherited. No worldwide priority claim.'}
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
