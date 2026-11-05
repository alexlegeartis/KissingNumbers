"""Replay an octad exchange after independently verifying the complete seed."""
import argparse,itertools,json,runpy,sys
from pathlib import Path
p=Path(__file__).resolve().parent;root=p.parents[1];sys.path.insert(0,str(root));from discovery.core import digest
runpy.run_path(str(p/'verify_seed.py'))
a=argparse.ArgumentParser();a.add_argument('proposal');args=a.parse_args();c=json.load(open(p/'seed32.json'));F=c['full_even_supports'];E=[tuple(x) for x in c['exceptional_antipodal_lines']];z=json.load(open(args.proposal));data=z['data'];assert data['backend']=='finite_field_signed_octad_exchange';assert data['source_supports_sha256']==digest(F) and data['source_exceptional_sha256']==digest(E)
new=[tuple(x) for x in data['candidate_lines']];removed=[tuple(x) for x in data['removed_lines']];assert len(set(new))==len(new) and len(set(removed))==len(removed)
def valid(s,w):return type(s) is int and type(w) is int and 0<=s<1<<32 and s.bit_count()==8 and not(w&~s) and not(w&(s&-s))
def bad(s,v,t,w):
 ov=s&t;return abs(ov.bit_count()-2*(ov&(v^w)).bit_count())>4
assert all(valid(*x) for x in new+removed)
Fset=set(F);Eset=set(E)
for s,w in removed:assert (s in Fset and w.bit_count()%2==0) or (s,w) in Eset
for s,w in new:assert s not in Fset and all(s!=t for t,v in E)
assert all(not bad(*a,*b) for a,b in itertools.combinations(new,2))
rm=set(removed)
for s,v in new:
 for t in F:
  if (s&t).bit_count()<=4:continue
  pos=[i for i in range(32) if t>>i&1]
  for bits in range(256):
   if bits&1 or bits.bit_count()%2:continue
   w=sum(1<<pos[j] for j in range(8) if bits>>j&1)
   assert not bad(s,v,t,w) or (t,w) in rm
 for t,w in E:assert not bad(s,v,t,w) or (t,w) in rm
gain=len(new)-len(removed);assert z['facts']['gain_lines']==gain;print('PASS exact exchange:',352800+2*gain,'points, gain',2*gain)
