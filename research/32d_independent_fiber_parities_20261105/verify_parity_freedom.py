"""Exact parity freedom on the independently verified 32–35D certificates."""
import itertools,json,importlib.util
from pathlib import Path
root=Path(__file__).resolve().parents[2]/'verifications/improved/dim32-35-signed-layers/source/classic'
spec=importlib.util.spec_from_file_location('classic_signed_verifier',root/'verify.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
expected={32:(352800,1704,1704),33:(370368,1848,1848),34:(390982,2013,1878),35:(415310,2202,2063)}
for n,(points,total,free_expected) in expected.items():
 p=root/f'certificate{n}.json';c=json.loads(p.read_text());assert module.verify(p)['verified_points']==points
 F=c['full_even_supports'];E=[tuple(x) for x in c['exceptional_antipodal_lines']]
 assert len(F)==total and len(set(F))==total and all(s.bit_count()==8 for s in F)
 assert all((s&t).bit_count()<=4 for s,t in itertools.combinations(F,2))
 assert all((s&t).bit_count()<=4 for s in F for t,v in E)
 for parity in (0,1):
  signs=[x for x in range(256) if x.bit_count()%2==parity];assert len(signs)==128
  assert all((x^y).bit_count()>=2 for x,y in itertools.combinations(signs,2))
 heavy=c.get('heavy_pair_extension');heads=[]
 if heavy:
  a,b=heavy['coordinates']
  for bs,ts in itertools.product((-1,1),repeat=2):
   q=[bs]*n;q[a]=3*ts;q[b]=-3*ts;heads.append(q)
 free=0
 for s in F:
  allowed=[];pos=[j for j in range(n) if s>>j&1]
  for parity in (0,1):
   okay=True
   for q in heads:
    v=[q[j] for j in pos];mx=sum(map(abs,v))
    if sum(t<0 for t in v)%2!=parity:mx-=2*min(map(abs,v))
    dot=2*mx
    if 4*dot*dot>32*(n+16):okay=False
   allowed.append(okay)
  assert allowed[0]
  free+=allowed[1]
 assert free==free_expected
 print(f'PASS {n}D: {points} points; {free} independently flippable full fibers; 2^{free} configurations')
