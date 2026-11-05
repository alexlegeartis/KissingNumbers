"""Search whole signed-octad families using nonlinear GF(32) permutations."""
import argparse,json,pathlib,sys,tempfile,importlib.util
ROOT=pathlib.Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from discovery.core import RunContext
from discovery.engines.support_recombination import search,certify
p=argparse.ArgumentParser();p.add_argument('--dimension',type=int,choices=(32,33,34,35),default=32);p.add_argument('--seconds',type=float,default=120);p.add_argument('--output',type=pathlib.Path);args=p.parse_args();d=args.dimension
location=ROOT/'verifications/improved/dim32-35-signed-layers/source/classic';source=json.loads((location/f'certificate{d}.json').read_text());spec=importlib.util.spec_from_file_location('independent_signed_octad_verifier',location/'verify.py');verifier=importlib.util.module_from_spec(spec);spec.loader.exec_module(verifier);verifier.verify(location/f'certificate{d}.json')
z=json.loads(pathlib.Path(__file__).with_name('permutations.json').read_text())
def mul(a,b):
 out=0
 while b:
  if b&1:out^=a
  a<<=1
  if a&32:a^=37
  b>>=1
 return out
def power(a,k):
 out=1
 while k:
  if k&1:out=mul(out,a)
  a=mul(a,a);k//=2
 return out
for (a,b,c,e),perm in zip(z['metadata'],z['permutations']):
 assert perm==[power(x,a)^mul(c,power(x,b))^mul(e,x) for x in range(32)]
perms=[v+list(range(32,d)) for v in z['permutations']];protected=sorted(set(v[0] for v in source['exceptional_antipodal_lines']))
ctx=RunContext(ROOT,{'max_rank_work':500000000,'max_pairs':10000000,'max_codewords':50000,'solver_seconds':max(args.seconds,1)})
a=search([] ,{'dimension':d,'source':source['full_even_supports'],'protected_supports':protected,'permutations':perms,'minimum_gain':1,'seconds':args.seconds},ctx);code=certify([a],{},ctx)
print(json.dumps({'dimension':d,'old_support_count':len(source['full_even_supports']),'candidate_count':a.data['candidate_count'],'solver_status':a.data['solver_status'],'verified_support_count':code.facts['cardinality']},indent=2))
if a.data['selected']:
 source['full_even_supports']=a.data['selected'];source['comparison']={'kissing_points':source['points']};source['points']+=128*(len(a.data['selected'])-len(a.data['source']))
 with tempfile.TemporaryDirectory() as t:
  path=pathlib.Path(t)/'certificate.json';path.write_text(json.dumps(source));print(json.dumps({'independent_kissing_verification':verifier.verify(path)},indent=2))
 if args.output:args.output.write_text(json.dumps(source,indent=2))
