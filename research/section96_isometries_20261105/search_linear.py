import itertools, numpy as np, subprocess, re, sys
ng=int(sys.argv[1]) if len(sys.argv)>1 else 4
m=int(sys.argv[2]) if len(sys.argv)>2 else 4
s=6
u=int(sys.argv[3]) if len(sys.argv)>3 else 15
v=int(sys.argv[4]) if len(sys.argv)>4 else 240
rounds=int(sys.argv[5]) if len(sys.argv)>5 else 500
seed=int(sys.argv[6]) if len(sys.argv)>6 else 961361+ng
assert rounds>0 and seed>=0
assert m==4 and 0<u<256 and 0<v<256 and u!=v
assert 0<=s<=12-m
assert 1 <= m <= 10
assert 1 <= ng <= 24
assignment=np.arange(24)%ng
from scipy.optimize import milp, Bounds, LinearConstraint
from scipy.sparse import csr_matrix
rng=np.random.default_rng(seed);g=np.array([7,1,4,4,2,7,3,3,4,7,2,1]);G=np.zeros((12,24),dtype=int)
for i in range(12):G[i,i:i+12]=g;G[i,23]=3
offset=np.zeros((3,24),dtype=int)
gray=np.array([sum(((((z>>2)&1)^((z&1)*x)^(((z>>1)&1)*y))<<(2*x+y)) for x in (0,1) for y in (0,1)) for z in range(8)])
wt=np.array([int(i).bit_count() for i in range(16)]);menu=[]
pairs=[(0,15),(3,12),(5,10),(6,9)]
for perm in itertools.permutations(range(4)):
 for flips in range(16):
  mapping={}
  for i,(a,b) in enumerate(pairs):
   pair=pairs[perm[i]];f=(flips>>i)&1;mapping[a]=pair[f];mapping[b]=pair[1-f]
  menu.append([mapping[int(x)] for x in gray])
menu=np.unique(np.array(menu,dtype=np.uint8),axis=0);menu=np.r_[menu,menu^1];P=len(menu);assert P==768
# Exact preservation of every within-cell distance is essential.
for item in menu:
 assert all(wt[item[a]^item[b]]==wt[gray[a]^gray[b]] for a in range(8) for b in range(8))
allowed=np.ones((ng,P));allowed[np.unique(assignment[:m]),P//2:]=0
rows=[];lb=[];ub=[]
for j in range(ng):
 row=np.zeros((ng,P));row[j,:]=1;rows.append(row.ravel());lb.append(1);ub.append(1)
for it in range(rounds):
 r=milp(rng.normal(size=ng*P),integrality=np.ones(ng*P),bounds=Bounds(np.zeros(ng*P),allowed.ravel()),constraints=LinearConstraint(csr_matrix(np.asarray(rows)),lb,ub),options={'time_limit':120})
 if r.x is None:print('no incumbent',it,r.message,flush=True);break
 choices=np.argmax(r.x.reshape(ng,P),axis=1);mapping=menu[choices[assignment]];mapping.tofile(f'search96_sectionlinear_{ng}_{m}_{u}_{v}_current.bin')
 matched=None
 for odd in (0,):
  for ref in range(1):
   p=subprocess.run(['./oracle_linear',f'search96_sectionlinear_{ng}_{m}_{u}_{v}_current.bin','200',str(m),str(u),str(v)],capture_output=True,text=True)
   if p.returncode != 0:raise RuntimeError(f'oracle failed ({p.returncode}): {p.stderr}')
   matched=re.search(r'WITNESS A ([0-7]{12}) B ([0-7]{12}) distance (\d+)',p.stdout)
   if matched: break
  if matched:break
 if not matched:
  print('UNRESOLVED',it,'choices',choices.tolist(),p.stdout.strip(),flush=True);np.save(f'search96_sectionlinear_{ng}_{m}_{u}_{v}_survivor_{it}.npy',mapping);continue
 a=np.array(list(matched[1]),dtype=int);b=np.array(list(matched[2]),dtype=int);assert sum((int(a[k])%2)<<(k-m) for k in range(m,12)) in (0,u,v,u^v);x=(a@G)%8;y=(b@G+sum((ref>>i&1)*offset[i] for i in range(3)))%8
 d=int(wt[mapping[np.arange(24),x]^gray[y]^odd].sum());assert d==int(matched[3]) and d<24
 assert np.all(wt[mapping[np.arange(m),x[:m]]]==2) and np.all(wt[gray[y[:m]]]==2)
 per_block=wt[menu[:,x].T^gray[y,None]^odd].astype(int)
 for j in range(m):per_block[j]+=24*(wt[menu[:,x[j]]]!=2)
 row=np.zeros((ng,P),dtype=int)
 for j in range(24):row[assignment[j]]+=per_block[j]
 rows.append(row.ravel());lb.append(24);ub.append(np.inf)
 print('witness',it,'distance',d,'cuts',len(rows)-ng,flush=True)
