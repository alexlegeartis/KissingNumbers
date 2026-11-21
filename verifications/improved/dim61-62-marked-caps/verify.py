"""Exact append/replace verification for distinguished-amplitude supports."""
from pathlib import Path
from fractions import Fraction as Q
import importlib.util,itertools,json,sys
HERE=Path(__file__).resolve().parent
def upstream():
 v=HERE.parent/'dim61-66-nonary-affine-caps'/'verify.py'
 if not v.is_file():v=HERE.parent/'mid_49_69_20261121/verifications/improved/dim61-66-nonary-affine-caps/verify.py'
 spec=importlib.util.spec_from_file_location('predecessor_nonary',v)
 p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)
 for rel,sha in p.SOURCE_SHA256.items():p.pin(p.V/rel,sha)
 p.pin(v.with_name('certificate.json'),'d2f7743669ca5ef39a6f7a58fdffc859bbc84550c815fcd89a6d958cb8b3f6de')
 return p
def phase():
 # Conway–Sloane (1994), Theorem 4: cyclic shifts and negation over Z/4;
 # exact phase check is independent of that theorem.
 seeds="001122 002211 010212 020121 012021 011332 012313 013123 013231 021133 311111 200000 022222".split()
 gray=(0,1,3,2)
 A=sorted({sum(gray[sg*int(t[(i+k)%6])%4]<<(2*i) for i in range(6))
           for t in seeds for k in range(6) for sg in (1,-1)})
 assert len(A)==144 and min((a^b).bit_count()for a,b in itertools.combinations(A,2))==4
 return [(z,a) for z in (1,-1) for a in A]

def newcount(A,H,mark,lo):
 I=[i for i,p in enumerate(lo)if p&H]
 return sum(2*z*bool(H&mark)+sum(1-2*(a>>i&1)for i in I)==0 for z,a in A)
def main():
 assert __debug__
 p=upstream()
 o60,*_=p.setup()
 C=p.read(p.HERE/'certificate.json'); A=phase()
 gray=C['gray_map']; seeds=C['pentacode_seeds']
 J={sum(gray[(z*int(t[(i+j)%5]))%4]<<(2*i)for i in range(5))>>1 for t in seeds for j in range(5)for z in (1,-1)}
 assert len(J)==40 and min((a^b).bit_count()for a,b in itertools.combinations(J,2))==3
 best_file=Path(sys.argv[1]) if len(sys.argv)>1 else HERE/'certificate.json'
 data=json.loads(best_file.read_text());out={}
 assert set(data)=={'61','62'}
 refill_path=HERE/'refill.json'
 extra_stock=json.loads(refill_path.read_text())
 assert set(extra_stock)==set(data)
 assert 4**2+12*2**2==64 and 8*4==32
 assert 64-2*4*4==32 and 64-2*16==32
 assert 28**2*Q(17,16)<32**2 and 28<=32 # body-only: excludes all exceptional tail coordinates.
 N=Q(8,3)
 assert 9*N*N==64 and 4*N*N<32 and (9-2*3)*N*N<32
 assert 6*2*N==32 and 3*4*N==32 and 8*N<32
 assert 11*N<32 and (9*N)**2*Q(17,16)<32**2
 assert (8*Q(9,8)+Q(1,2))*N<32  # 61D radical latitude with occupied tail.
 for d0,pairs in data.items():
  d=int(d0); D=60
  assert d in (61,62)
  H=0
  S=[(int(m,16),int(l,16))for m,l in pairs]
  assert S and len(S)==len(set(S))
  assert all(m.bit_count()==1 and l.bit_count()==12 and m&l==0 and (m|l)<1<<D for m,l in S)
  bits=lambda l:[1<<i for i in range(D)if l>>i&1]
  count=sum(newcount(A,H,m,bits(l))for m,l in S)
  assert count==len(S)*288
  removed=set()
  for m,l in S:
   assert max(2*(m&k).bit_count()+(l&k).bit_count()for k in o60)<=8
   idx=m.bit_length()-1
   L=[u.bit_length()-1 for u in bits(l)]
   # A four-support block violates iff it contains the mark and three low entries.
   assert all(i^j^k!=idx for i,j,k in itertools.combinations(L,3))
   for a,k in enumerate(C['supports'][d0]):
    if 2*(m&k).bit_count()+(l&k).bit_count()>6:removed.add(a)
  for (m,l),(a,b) in itertools.combinations(S,2):
   assert 4*(m&a).bit_count()+2*((m&b).bit_count()+(l&a).bit_count())+(l&b).bit_count()<=8
  extra=list(map(int,extra_stock.get(d0,[])))
  if extra:
   assert d in (61,62)
   assert len(extra)==len(set(extra))
   kept=[k for i,k in enumerate(C['supports'][d0])if i not in removed]
   assert not set(kept)&set(extra)
   for i,k in enumerate(extra):
    assert k.bit_count()==9 and 0<k<1<<d
    assert all((k&r).bit_count()<=4 for r in kept+extra[:i])
    assert all((k&r).bit_count()<=6 for r in o60)
    J=[j for j in range(d)if k>>j&1];diff=[a^b for a,b in itertools.combinations(J,2)]
    assert len(set(diff))==36
    assert all(2*(m&k).bit_count()+(l&k).bit_count()<=6 for m,l in S)
   # The inherited checked 40-word phase works again for these supports.
   count+=40*len(extra)
  lost=40*len(removed)
  gain=count-lost
  assert gain>0
  out[d]=dict(total=(p.BASE[d]+40*len(C['supports'][d0]))+gain,gain=gain,added=count,lost=lost,removed_supports=len(removed),templates=len(S),refilled_supports=len(extra))
 return out
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
