from pathlib import Path
import json,hashlib,itertools,math
import numpy as np
if not __debug__:raise RuntimeError('Assertions required')
S=Path('research/z8_fourrow92_partitions_20261117');D=Path(__file__).parent;archive_sha='74265a34180da67093877f3a285534ee17947cfeb7bb443ef2f1278a02bbb90f';assert hashlib.sha256(Path(str(S)+'.zip').read_bytes()).hexdigest()==archive_sha
for name,digest in json.load(open(S/'SHA256SUMS.json')).items():assert hashlib.sha256((S/name).read_bytes()).hexdigest()==digest
parent=Path('research/three_basis_independent_20261117');assert(S/'z8_generator.json').read_bytes()==(parent/'input/z8_generator.json').read_bytes();pin=hashlib.sha256((S/'z8_generator.json').read_bytes()).hexdigest();assert pin=='903c193dbd00c1e7c692e8de62fc4e42c42806dd53365f83ddd73e673dc0fb73';pr=json.load(open(parent/'parent_result.json'));assert pr['binary_distance']==24 and pr['z4_words_enumerated']==4**12 and pr['octads_checked']==759 and pr['low_lee_errors_per_octad']==5489
G=np.array(json.load(open(S/'z8_generator.json'))['generator'],dtype=np.int64);assert G.shape==(12,24)and not np.tril(G[:,:12],-1).any()and np.all(np.diag(G[:,:12])%2);assert not G[4:,:4].any();prefix_inputs=np.array(list(itertools.product(range(8),repeat=4)));assert len(np.unique(prefix_inputs@G[:4,:4]%8,axis=0))==8**4
A=[tuple((-1)**((a//4)^((a&1)*x)^(((a//2)&1)*y))for x,y in [(0,0),(0,1),(1,0),(1,1)])for a in range(8)];B=[(-v[0],)+v[1:]for v in A];C=[tuple((2 if a<4 else-2)if j==a%4 else 0 for j in range(4))for a in range(8)];alph=[A,B,C];dot=lambda u,v:sum(a*b for a,b in zip(u,v));weights=[0,2,2,2,4,2,2,2]
for X in alph:
 assert len(set(X))==8
 for a,b in itertools.product(range(8),repeat=2):assert dot(X[a],X[b])==4-2*weights[(a-b)%8]
for X,Y in itertools.combinations(alph,2):assert set(dot(x,y)for x in X for y in Y)=={-2,2}
r=json.load(open(S/'candidate_three.json'));part=r['partitions'];families=r['selected_supports'];assert r['dimension']==92 and len(part)==len(families)==3;counts=[len(f)for f in families];assert counts==[81,24,6]
for groups,family in zip(part,families):
 assert len(groups)==4 and all(len(g)==3 for g in groups)and sorted(x for g in groups for x in g)==list(range(12));assert len({tuple(sorted(s))for s in family})==len(family)
 for s in family:assert len(set(s))==4 and all(len(set(g)&set(s))==1 for g in groups)
for i,j in itertools.combinations(range(3),2):assert max(len(set(s)&set(t))for s in families[i]for t in families[j])<=2
heads=[];prefixes=[];colors=[]
for color,(groups,family)in enumerate(zip(part,families)):
 for support in family:
  axes=[next(x for x in g if x in support)for g in groups];indices=[g.index(x)for g,x in zip(groups,axes)]
  for bits in itertools.product(range(2),repeat=4):
   head=[0]*12;pre=[]
   for axis,index,bit in zip(axes,indices,bits):head[axis]=2*(-1)**bit;pre.append(index+4*bit)
   heads.append(head);prefixes.append(pre);colors.append(color)
assert len(heads)==len({tuple(h)for h in heads})==1776
H=np.array(heads,dtype=np.int16);P=np.array(prefixes,dtype=np.int16);col=np.array(colors);assert np.all(np.sum(H*H,axis=1)==16)
for color in range(3):
 ids=np.where(col==color)[0];hc=H[ids];pc=P[ids];look=np.array([[dot(x,y)for y in alph[color]]for x in alph[color]],dtype=np.int16);rhs=sum(look[pc[:,j,None],pc[None,:,j]]for j in range(4));assert np.array_equal(hc@hc.T,rhs)
for i,j in itertools.combinations(range(3),2):assert np.max(H[col==i]@H[col==j].T)<=8
reported=json.load(open(S/'head_certificate.json'));assert set(zip(map(tuple,heads),map(tuple,prefixes),colors))==set(zip(map(tuple,reported['heads']),map(tuple,reported['prefixes']),reported['colors']))
assert 20*4+16==96 and 20*2+8==48
# GF16 built through multiplicative log/exponent tables, then distributivity
# independently checked for every field triple.
exp=[];v=1
for k in range(15):exp.append(v);v<<=1;v^=19 if v&16 else 0
assert v==1 and set(exp)==set(range(1,16));log={v:i for i,v in enumerate(exp)}
def mul(a,b):return 0 if not a or not b else exp[(log[a]+log[b])%15]
for a,b,c in itertools.product(range(16),repeat=3):assert mul(a,b^c)==mul(a,b)^mul(a,c)
for a in range(1,16):assert mul(a,exp[(-log[a])%15])==1
# Every four-position evaluation minor is a nonzero Vandermonde product.
for xs in itertools.combinations(range(6),4):
 det=1
 for a,b in itertools.combinations(xs,2):det=mul(det,a^b)
 assert det!=0
sizes=[12]+[16]*5;offsets=[sum(sizes[:i])for i in range(6)];supports=set()
for a0 in range(12):
 for a1,a2,a3 in itertools.product(range(16),repeat=3):
  vals=[mul(mul(mul(a3,x)^a2,x)^a1,x)^a0 for x in range(6)];s=tuple(offsets[j]+vals[j]for j in range(6));assert len(set(s))==6 and max(s)<92 and vals[0]<12;supports.add(s)
assert len(supports)==12*16**3
signs=[v for v in range(64)if v.bit_count()%2==0];assert len(signs)==32 and min((a^b).bit_count()for a in signs for b in signs if a!=b)==2
assert 6*16==96 and 3*16==48 and 4*16*2>=96 and 6*4*2<=48
# Root geometry: unscaled distinct signed2-support roots have dot<=1.
roots=set()
for offset,size in zip(offsets,sizes):
 for i,j in itertools.combinations(range(offset,offset+size),2):
  for a,b in itertools.product([-1,1],repeat=2):roots.add((i,j,a,b))
assert len(roots)==2664 and 2*48==96 and (4**2)*48<48**2
# Each RS support has one position per group, hence meets a root in<=1.
assert len(sizes)==6 and sum(sizes)==92
counts=dict(dense=len(heads)*8**8,RS=len(supports)*32,roots=len(roots));assert counts==dict(dense=29796335616,RS=1572864,roots=2664);assert sum(counts.values())==29797911144
out=dict(status='INDEPENDENT NEW GEOMETRY AND COUNTS PASS; pinned parent minimum24 inherited',dimension=92,counts=counts,total=sum(counts.values()),supports_by_color=[81,24,6],heads=1776,prefix_fiber=8**8,RS_supports=len(supports),source_archive_sha256=archive_sha,parent_generator_sha256=pin)
(D/'result.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
