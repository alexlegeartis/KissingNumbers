import argparse,hashlib,itertools,json,sys,subprocess,tempfile
from pathlib import Path
import numpy as np
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import coo_matrix

HERE=Path(__file__).resolve().parent
PACKAGE=HERE.parent
parser=argparse.ArgumentParser(description='Search the finite full-star family; verification is separate.')
parser.add_argument('dimension',type=int,choices=(33,34,35))
parser.add_argument('--output',required=True,type=Path)
args=parser.parse_args()
if args.output.resolve()==PACKAGE or PACKAGE in args.output.resolve().parents:
 parser.error('write search proposals outside the verification package')
r=PACKAGE/'source'
for d in (args.dimension,):
 c=json.loads((r/('multi35/certificate35.json' if d==35 else 'classic/certificate%d.json'%d)).read_text())
 E=np.array([[(-2 if s>>j&1 else 2) if m>>j&1 else 0 for j in range(d)] for m,s in c.get('exceptional_antipodal_lines',[])],dtype=np.int64).reshape(-1,d)
 heads=[];old_removed=set()
 if d==34:
  a,b=c['heavy_pair_extension']['coordinates'];old_removed={(a,b,-1)}
  for t in (-1,1):
   q=[1]*d;q[a],q[b]=3*t,-3*t;heads.append(q)
 elif d==35:heads=c['head_antipodal_lines'];old_removed=set(map(tuple,c['removed_root_lines']))
 heads=np.array(heads,dtype=np.int64).reshape(-1,d)
 stars=[];vectors=[];lines=[]
 full=c.get('full_even_supports',[mask for mask,parity in c.get('full_parity_supports',[])])
 with tempfile.TemporaryDirectory(prefix='root_star_search_') as scratch:
  scratch=Path(scratch);source=scratch/'supports.txt';stars_file=scratch/'stars.txt'
  source.write_text(''.join('%d\n'%mask for mask in full))
  subprocess.run(['c++','-O3','-std=c++17',str(HERE/'enumerate_stars.cpp'),'-o',str(scratch/'enumerate')],check=True)
  subprocess.run([str(scratch/'enumerate'),str(d),str(source),str(stars_file)],check=True)
  candidates=stars_file.read_text().splitlines()
 for line in candidates:
  a,m=map(int,line.split());leaves=[j for j in range(d) if m>>j&1];H=[]
  for sa in (-1,1):
   for mask in range(16):
    if bin(mask).count('1')%2:continue
    v=[0]*d;v[a]=4*sa
    for k,j in enumerate(leaves):v[j]=2-4*((mask>>k)&1)
    H.append(v)
  H=np.array(H,dtype=np.int64)
  if len(E) and np.abs(H@E.T).max()>16:continue
  if len(heads) and np.any(4*(H@heads.T)**2>32*np.sum(heads*heads,axis=1)):continue
  stars.append((a,leaves));vectors.append(H)
  lines.append({tuple(sorted((a,b)))+(sign,) for b in leaves for sign in (-1,1)}-old_removed)
 all_lines=sorted(set.union(*lines));lookup={line:len(stars)+i for i,line in enumerate(all_lines)}
 objective=np.array([-16]*len(stars)+[2]*len(all_lines),dtype=float)
 rr=[];cc=[];vv=[];upper=[]
 def constraint(items,bound):
  row=len(upper);upper.append(bound)
  for col,val in items:rr.append(row);cc.append(col);vv.append(val)
 for i,H in enumerate(vectors):
  for j,K in enumerate(vectors[:i]):
   if np.abs(H@K.T).max()>16:constraint([(i,1),(j,1)],1)
 for i,required in enumerate(lines):
  for line in required:constraint([(i,1),(lookup[line],-1)],0)
 matrix=coo_matrix((vv,(rr,cc)),shape=(len(upper),len(objective))).tocsc()
 result=milp(objective,integrality=np.ones(len(objective)),bounds=Bounds(0,1),constraints=LinearConstraint(matrix,-np.inf,np.array(upper)),options={'time_limit':30})
 if result.x is None:raise RuntimeError('The search returned no feasible witness')
 selected=[stars[i] for i in range(len(stars)) if result.x[i]>.5]
 removed=set.union(*(lines[i] for i in range(len(stars)) if result.x[i]>.5))
 gain=16*len(selected)-2*len(removed)
 print('dimension',d,'stars',len(selected),'gain',gain,'objective optimal',result.success,selected,flush=True)
 baseline='source/multi35/certificate35.json' if d==35 else 'source/classic/certificate%d.json'%d
 certificate=dict(format='root_star_trade_v1',dimension=d,baseline=baseline,
                  baseline_sha256=hashlib.sha256((PACKAGE/baseline).read_bytes()).hexdigest(),
                  stars=[dict(center=a,leaves=leaves) for a,leaves in selected],
                  removed_root_lines=sorted(removed),points=c['points']+gain)
 args.output.parent.mkdir(parents=True,exist_ok=True)
 args.output.write_text(json.dumps(certificate,indent=2)+'\n')
