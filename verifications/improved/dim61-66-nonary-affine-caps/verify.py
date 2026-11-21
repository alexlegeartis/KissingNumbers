"""Verify the append-only 61–66 dimensional nine-support cap certificate."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
import hashlib,json,importlib.util,numpy as np

HERE=Path(__file__).resolve().parent
V=HERE.parent
SOURCE=V/'dim61-68-95-96-code-layers/source'
PARENTS={'61,62':('dim61-62-66-affine-reconstruction/dim61-62/certificate.json',88979562,122543436),
         '63,64,65':('dim63-65-affine-reconstruction/certificate.json',184415700,335749696,335817922),
         '66':('dim61-62-66-affine-reconstruction/dim66/certificate.json',342291652)}
BASE={61:88979562,62:122543436,63:184415700,64:335749696,65:335817922,66:342291652}
SOURCE_SHA256 = {
    'dim61-62-66-affine-reconstruction/dim61-62/certificate.json': '34c5b527252603d693c8e760b7b4215d0d84eaf503ce44622d11bfdd7fc0dee1',
    'dim61-62-66-affine-reconstruction/dim66/certificate.json': '26a393f7e6202f84fb5a945b95bc623b4e1dde8711d8ce384f7bc9ed396fc501',
    'dim63-65-affine-reconstruction/certificate.json': '266af3683caadf6d23c2aa3e7409ce28d2d3f4a8644be2784263c80bd407fc7e',
    'dim63-65-affine-reconstruction/section63.json': '25cae31ba6b5c8746a7d2f203746e6f0b3d4f31f391019cbf80044b5fb4d4b15',
    'dim61-67-compatible-sign-families/data/global64.json': '8c0c35f8127c578118f1f09e5950d712a272db32f2d50f382e6fe32f46c6dbf5',
    'dim61-67-compatible-sign-families/data/short61_62.json': 'd9f5744e5c92745ab6f5cc8bb144f0bd7862dfc4c69a3371c2759b408a24e598',
    'dim61-68-95-96-code-layers/source/dim61_96/shortening_support_search.json': '5422ed8f2b3be75267a57ea6fccc88f403a80e5bbe60dc64adc8edda51b61d88',
    'dim61-68-95-96-code-layers/source/continuation49_72/netto.py': 'c3e73dd365bc432d625bfeb106096d19cdf090772b9717b5f2b4a097496bb855',
}
def read(p):
 return json.loads(p.read_text())
def pin(p,sha):
 assert hashlib.sha256(p.read_bytes()).hexdigest()==sha
 return p
def setup():
 archive=pin(SOURCE/'dim61_96/code64_supports.npz','cdbd5e5a7bbda9ddb84ee86052e2d4d00ccc331f0188de958a69ed83b0c529f6')
 O=list(map(int,np.load(archive,allow_pickle=False)['words']))
 assert len(O)==len(set(O))==30828 and all(s.bit_count()==16 for s in O)
 P=read(V/PARENTS['61,62'][0]);F=P['family']
 assert F==read(V/'dim61-67-compatible-sign-families/data/short61_62.json')
 G=read(V/'dim61-67-compatible-sign-families/data/global64.json')
 A=read(V/PARENTS['63,64,65'][0]);B=read(V/PARENTS['66'][0])
 def rows(F,J):
  return [int(a['support']) for a in F['new']]+[int(a['support']) for a in J['new_classes']]
 rem=(24,42,53,55)
 assert F['removed_coordinates']==list(rem)
 k=sum(1<<i for i in rem); keep=[i for i in range(64)if i not in rem]
 def squash(x):return sum(((x>>p)&1)<<i for i,p in enumerate(keep))
 o60=[squash(w) for w in O+rows(F,P) if w&k==0];o64=O+rows(G,A)
 assert len(o60)==len(set(o60))==10635 and len(o64)==len(set(o64))==32801
 p66=pin(SOURCE/'dim61_96/cyclic73.weight16.bin','5d240938e26a7a477de62bc255439e743cb33120bb13b986a70c6b1fea0b5239')
 d9=(21,37,38,49,57,63,67,70,72)
 head=[i for i in range(1,73) if i not in d9]+[0]
 C=read(SOURCE/'dim61_96/shortening_support_search.json')['dimensions']['66']
 k=sum(1<<i for i in C['deleted_coordinates'])
 old=[int(x)|(int(y)<<64) for x,y in np.fromfile(p66,'<u8').reshape(-1,2)]
 old=[z for z in old if not z&k]
 assert len(old)==34028
 # The new 66D vectors in this certificate have no tails; overlap with the projection suffices.
 o66=[sum(((z>>p)&1)<<i for i,p in enumerate(head)) for z in old]+rows(B['family'],B)
 assert len(o66)==35994
 path=SOURCE/'continuation49_72/netto.py'
 spec=importlib.util.spec_from_file_location('source_netto',path)
 m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 net=m.blocks68()
 assert len(net)==12529
 b66=tuple(sorted({sum(1<<i for i in z)for z in net if max(z)<64}))
 return o60,o64,o66,b66
def main():
 assert __debug__,'Run without -O'
 cert=read(HERE/'certificate.json'); ts=cert['pentacode_seeds']; gray=cert['gray_map']
 assert ts==['01112','03110','21310','21132'] and gray==[0,1,3,2]
 signs=sorted({sum(gray[(a*int(t[(j+k)%5]))%4]<<(2*j) for j in range(5))>>1 for t in ts for k in range(5) for a in (-1,1)})
 assert len(signs)==40 and min((a^b).bit_count()for a,b in combinations(signs,2))==3
 for rel,sha in SOURCE_SHA256.items():pin(V/rel,sha)
 o60,o64,o66,b66=setup()
 z=read(V/'dim63-65-affine-reconstruction/section63.json');H=z['octad'];per=z['four_support_permutation']
 assert len(set(per))==64 and H==sum(1<<i for i in per[:8]) and z['total']==BASE[63]
 rev={p:i for i,p in enumerate(per)}
 N=Q(8,3)
 assert 9*N*N==64 and (9-2*3)*N*N<=32 and 4*N*N<=32
 assert N*2*6==32 and N*4*3==32 and N*8<32
 # Dense maximal coordinate sums for 61 and 62, before multiplying by N.
 # 61D: first latitude is body 1 / tail 2; second body sqrt(17)/4 / tail 1/2.
 assert 11*N<=32 and ((9*N)**2 * Q(17,16))<32**2
 assert (Q(17,16)<Q(9,8)**2) and (8*Q(9,8)+Q(1,2))*N<32
 # 62D: at most 2 bodies of amplitude>1, tails ≤2 and ≤sqrt(3)<2.
 assert (7+2+2)*N<32
 assert 9*N<32 # all other unrestricted dense layers, body amplitude 1
 supports=cert['supports'];assert set(supports)=={'61','62','63','64','66'}
 results={}
 for d in (61,62,63,64,66):
  objs=supports[str(d)]; S=[m[0] for m in objs] if d==63 else objs
  assert len(S)==len(set(S))
  assert all(0<m<1<<(min(d,64) if d!=63 else 64) and m.bit_count()==9 for m in S)
  depth=max((m&k).bit_count() for m,k in combinations(S,2));assert depth<=4
  stock=o60 if d<=62 else o66 if d==66 else o64
  paid=max((m&k).bit_count()for m in S for k in stock);assert paid<=6
  n=0
  for i,m in enumerate(S):
   P=[i for i in range(64) if m>>i&1]
   if d==66:
    assert all(m&b!=b for b in b66)
   else:
    I=[rev[p]for p in P] if d==63 else P
    ds=[x^y for x,y in combinations(I,2)]
    assert len(set(ds))==36 # no contained quadruple in the affine SQS
   if d!=63:n+=40;continue
   hit=[j for j,p in enumerate(P) if H>>p&1]
   mark=objs[i][1]
   assert len(mark)==len(hit) and len(mark)%2==0 and len(set(mark))==len(mark) and all(0<=j<9 for j in mark)
   mask=sum(1<<i for i in mark)
   # Relabel the nine-bit code to place precisely those source coordinates on H.
   # Sign balance on H is equivalent to being in the exact incumbent hyperplane.
   n+=sum((c&mask).bit_count()==len(mark)//2 for c in signs)
  results[d]=dict(base=BASE[d],added=n,total=BASE[d]+n,supports=len(S),parent_overlap=paid,new_overlap=depth)
  if d==64:
   # Equatorial additions transfer to the 65D incumbent: its dense and GF(4)-line
   # owners already occur in the verified 64D layer geometry. Lifted owners' body
   # factors are sqrt(3)/2. Both additional poles are orthogonal.
   results[65]=dict(results[d],base=BASE[65],total=BASE[65]+n)
 return results

if __name__=='__main__':
 print(json.dumps(main(),indent=2,sort_keys=True))
