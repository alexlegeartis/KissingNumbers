#!/usr/bin/env python3
"""Replay both exact checkers against mathematical corruptions in temp copies.

Standard library only. No production input or checker is modified. The altered
embedding tests refresh the checksum in an isolated checker copy so rejection
must reach the algebraic Gram calculation rather than only a hash guard.
"""
from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ThreadPoolExecutor
import copy,hashlib,json,shutil,subprocess,sys,tempfile

HERE=Path(__file__).resolve().parent
PIN='5f180d5d6c737010083e72210e16b92926ec2e79ed9a0891d0960a018b2d7842'

def main():
 if not __debug__:raise RuntimeError('Run without -O.')
 root=json.loads((HERE/'certificate.json').read_text())
 orth=json.loads((HERE/'orthonormal_certificate.json').read_text())
 embedding=json.loads((HERE/'embedding.json').read_text())
 controls=[]
 def case(checker,cert,name,mutate=None,alter_embedding=False,optimized=False):
  c=copy.deepcopy(cert)
  if mutate is not None:mutate(c)
  controls.append(dict(checker=checker,certificate=c,name=name,alter_embedding=alter_embedding,optimized=optimized))
 for checker,c in (('verify.py',root),('verify_orthonormal.py',orth)):
  case(checker,c,'unchanged_valid_certificate')
  case(checker,c,'false_claimed_bound',lambda d:d.__setitem__('bound',d['bound']+1))
  case(checker,c,'wrong_shell_population',lambda d:d.__setitem__('population',d['population']+1))
  case(checker,c,'false_moment_rhs',lambda d:d['rhs'].__setitem__(0,d['rhs'][0]+1))
  case(checker,c,'false_population_objective',lambda d:d['objective'].__setitem__(0,d['objective'][0]+1))
  case(checker,c,'changed_embedded_normal_hash_refreshed',alter_embedding=True)
  case(checker,c,'assertions_disabled',optimized=True)
 case('verify.py',root,'duplicate_selected_slice',lambda c:c['offsets'].__setitem__(1,c['offsets'][0]))
 case('verify.py',root,'slice_outside_proved_domain',lambda c:c['offsets'].__setitem__(0,[5,0,0,0]))
 case('verify.py',root,'omitted_profile_cell',lambda c:c['cells'].pop())
 case('verify.py',root,'false_moment_coefficient',lambda c:c['rows'][0][0].__setitem__(1,c['rows'][0][0][1]+1))
 case('verify.py',root,'false_dual_multiplier',lambda c:c['multipliers'].__setitem__(0,str(F(c['multipliers'][0])+1)))
 case('verify.py',root,'negative_squared_scale',lambda c:c['alpha_squared'].__setitem__(0,'-1'))
 def cross_root(c):
  # Profile zero has projected norm squared 1/4; this keeps its image norm
  # and same-slice feasibility, but violates a genuine cross-slice pair.
  c['heights'][0]='-1/2';c['alpha_squared'][0]='1'
 case('verify.py',root,'cross_slice_failure_with_norms_preserved',cross_root)
 case('verify_orthonormal.py',orth,'duplicate_selected_slice',lambda c:c['profiles'].__setitem__(1,c['profiles'][0]))
 case('verify_orthonormal.py',orth,'wrong_dual_lattice_parity',lambda c:c['profiles'].__setitem__(0,[0,0,0,1]))
 case('verify_orthonormal.py',orth,'omitted_profile_cell',lambda c:c['groups'][0].pop())
 case('verify_orthonormal.py',orth,'false_moment_coefficient',lambda c:c['rows'][0].__setitem__(0,c['rows'][0][0]+1))
 case('verify_orthonormal.py',orth,'false_dual_multiplier',lambda c:c['dual_multipliers'].__setitem__(1,c['dual_multipliers'][1]+1))
 case('verify_orthonormal.py',orth,'negative_dual_denominator',lambda c:c.__setitem__('dual_denominator',-c['dual_denominator']))
 case('verify_orthonormal.py',orth,'cross_slice_failure_with_norms_preserved',lambda c:c['heights'].__setitem__(0,'-1/2'))
 with tempfile.TemporaryDirectory(prefix='gamma72-69-controls-') as td:
  td=Path(td)
  def execute(item):
   folder=td/(item['checker'].replace('.py','')+'-'+item['name']);folder.mkdir()
   code=(HERE/item['checker']).read_text();e=copy.deepcopy(embedding)
   if item['alter_embedding']:
    e['normals'][0][0]+=1
    raw=(json.dumps(e,indent=2)+'\n').encode();newhash=hashlib.sha256(raw).hexdigest()
    assert PIN in code;code=code.replace(PIN,newhash)
   else:raw=(HERE/'embedding.json').read_bytes()
   (folder/'embedding.json').write_bytes(raw);(folder/item['checker']).write_text(code)
   (folder/'test.json').write_text(json.dumps(item['certificate']))
   command=[sys.executable]+(['-O'] if item['optimized'] else [])+[str(folder/item['checker']),str(folder/'test.json')]
   result=subprocess.run(command,capture_output=True,text=True,timeout=60)
   expected_pass=item['name']=='unchanged_valid_certificate'
   assert (result.returncode==0)==expected_pass,(item['checker'],item['name'],result.stdout,result.stderr)
   return dict(checker=item['checker'],case=item['name'],expected='PASS' if expected_pass else 'REJECT',
               actual='PASS' if result.returncode==0 else 'REJECT',returncode=result.returncode,
               checksum_refreshed=item['alter_embedding'])
  # At most two processes, inheriting the caller's CPU affinity.
  with ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(execute,controls))
 report={'status':'PASS','valid_baselines':2,'rejected_corruptions':len(results)-2,
         'source_sha256':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in ('verify.py','verify_orthonormal.py','certificate.json','orthonormal_certificate.json','embedding.json')},
         'cross_slice_witness':{'changed_slice_index':0,'changed_height':'-1/2','squared_scale':'1','other_slice_index':5,'squared_geometric_gap':'-31/16','within_slice_norms_preserved':True},
         'tests':results,'scope':'Corruption coverage for the exact finite proof. Published Gamma72 minimum, shell population, design and catalogue identity remain explicit premises.'}
 (HERE/'rejection_controls.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps(report,indent=2))

if __name__=='__main__':main()
