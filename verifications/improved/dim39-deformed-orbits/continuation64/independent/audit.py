import json,hashlib,math,itertools,subprocess,time,pathlib,sys
HERE=pathlib.Path(__file__).resolve().parent
SRC=pathlib.Path(sys.argv[1]).resolve() if len(sys.argv)>1 else HERE.parent
if not __debug__:raise RuntimeError('Run without -O.')
t=time.time()
c=json.loads((SRC/'orbit_certificate.json').read_text()); b=json.loads((SRC/'base_certificate.json').read_text())
sha=hashlib.sha256((SRC/'orbit_certificate.json').read_bytes()).hexdigest()
assert sha=='bfea4089709e11e7b810da37f751868552ed0c39d34256c410e386531f5bc401'
V=c['representatives'];G=c['group_words'];L=c['norm_lower'];U=c['norm_upper'];d=c['dimension']
assert d==39 and len(V)==5127 and len(G)==64 and len(set(G))==64 and 0 in G
assert all(type(g)==int and 0<=g<2**d for g in G)
assert all(g^h in G for g in G for h in G)
M=max(abs(x) for v in V for x in v);assert 2*d*M*M<2**63 and L*L<2**63
rays=set();norms=[];slacks=[10**100]*3
for v in V:
 assert len(v)==d and all(type(x)==int and x!=0 for x in v)
 norm=sum(x*x for x in v);norms.append(norm);assert L*L<=norm<=U*U
 av=sorted(map(abs,v),reverse=True)
 ss=[2*L*L-sum(av[:8])**2,L*L-2*sum(av[:2])**2,2*L*L-(2*av[0]+sum(av[1:5]))**2]
 assert min(ss)>=0
 slacks=[min(x,y) for x,y in zip(slacks,ss)]
 gcd=math.gcd(*v);sg=sum((x<0)<<i for i,x in enumerate(v));key=(tuple(abs(x)//gcd for x in v),min(sg^g for g in G))
 assert key not in rays;rays.add(key)
# Read supports as left-to-right coordinates, independently check each edited row.
lines=(SRC/'source/supports_3424.txt').read_text().splitlines();assert len(lines)==3424
supports=[]
for line in lines:
 assert len(line)==39 and set(line)<=set('01')
 supports.append(tuple(i for i,x in enumerate(line) if x=='1'))
for edit in b['support_replacements']:
 assert list(supports[edit['remove_index']])==edit['remove']
 supports[edit['remove_index']]=tuple(edit['add'])
assert all(len(s)==8 and len(set(s))==8 and all(0<=i<39 for i in s) for s in supports)
assert len(set(supports))==3424
masks=[sum(1<<i for i in s) for s in supports];maxintersection=0
for i,a in enumerate(masks):
 for bb in masks[:i]:
  overlap=(a&bb).bit_count();assert overlap<=4;maxintersection=max(maxintersection,overlap)
# Explicitly construct heavy and retained root vectors.
removed={tuple(sorted((s['center'],l))) for s in b['heavy_stars'] for l in s['leaves']}
roots=[]
for i,j in itertools.combinations(range(39),2):
 if (i,j) in removed:continue
 for a,bb in itertools.product([-4,4],repeat=2):
  r=[0]*39;r[i]=a;r[j]=bb;roots.append(tuple(r))
heavy=[]
for star in b['heavy_stars']:
 a=star['center'];leaves=star['leaves'];assert len(set([a]+leaves))==5
 for signs in itertools.product([-2,2],repeat=4):
  if sum(x<0 for x in signs)%2:continue
  for center in [-4,4]:
   v=[0]*39;v[a]=center
   for i,x in zip(leaves,signs):v[i]=x
   heavy.append(tuple(v))
assert len(roots)==len(set(roots))==2908 and len(heavy)==len(set(heavy))==64
assert all(sum(x*x for x in v)==32 for v in roots+heavy)
assert set(roots).isdisjoint(heavy)
dot=lambda a,bb:sum(x*y for x,y in zip(a,bb))
maxhh=max(dot(a,bb) for i,a in enumerate(heavy) for bb in heavy[:i]);assert maxhh<=16
maxhr=max(dot(a,r) for a in heavy for r in roots);assert maxhr<=16
# Maximize middle signs exactly using parity, not just a geometric assumption.
maxhm=-1000
for h in heavy:
 for s in supports:
  values=[2*h[k] for k in s];mx=sum(map(abs,values))
  if sum(x<0 for x in values)%2:mx-=2*min(map(abs,values))
  assert mx<=16;maxhm=max(maxhm,mx)
# Middle-middle: distinct supports intersection<=4, equal supports even parity differs>=2.
# Root-root: distinct same-support roots <=0, distinct supports overlap<=1.
# Root-middle: absolute dot<=4*2*2=16. All sparse norms²=32.
# Dense versus sparse: top8 bounds middle; top2 bounds roots;
# rearrangement 2*largest+next4 bounds every heavy vector, all group actions preserve abs values.
sparse=len(supports)*128+len(roots)+len(heavy);dense=len(V)*len(G)
assert sparse==441244 and dense==328128 and sparse+dense==c['claimed_size']==769372
print('Structural and all sparse checks passed; full relative product check starting.',flush=True)
subprocess.run(['g++','-O3','-march=native','-std=c++17',str(HERE/'relative.cpp'),'-o',str(HERE/'relative')],check=True)
data=f'{len(V)} {d} {len(G)} {L}\n'+' '.join(map(str,G))+'\n'+'\n'.join(' '.join(map(str,v)) for v in V)+'\n'
p=subprocess.run([str(HERE/'relative')],input=data,text=True,stdout=subprocess.PIPE,check=True)
ct,mx,i,j,g=map(int,p.stdout.split());assert ct==64*len(V)*(len(V)+1)//2-len(V)==841315065
result=dict(passed=True,certificate_sha256=sha,dimension=d,count=sparse+dense,dense=dense,sparse=sparse,relative_products=ct,maximum_dot=mx,maximum_witness=[i,j,g],allowed_dot=L*L//2,norm_squared_range=[min(norms),max(norms)],cross_slacks=slacks,max_support_intersection=maxintersection,heavy_heavy_max=maxhh,heavy_root_max=maxhr,heavy_middle_max=maxhm,coordinate_max=M,elapsed_seconds=time.time()-t)
(HERE/'audit_result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
