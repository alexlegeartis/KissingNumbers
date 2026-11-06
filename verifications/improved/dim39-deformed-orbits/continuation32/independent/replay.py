"""Replay the independent implementation from source, including two corruption tests."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,tempfile
import numpy as np
H=Path(__file__).resolve().parent
if not __debug__:raise RuntimeError('Run without -O.')
env=dict(os.environ,OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='1')
def run(*args):return subprocess.run(args,cwd=H,env=env,check=True,capture_output=True,text=True)
run(sys.executable,str(H/'prepare_dense.py'))
run('g++','-O3','-march=native','-fopenmp','-std=c++17',str(H/'dense_direct.cpp'),'-o',str(H/'dense_direct'))
p=run(str(H/'dense_direct'),str(H/'dense_input.bin'));d=json.loads(p.stdout)
c=json.loads((H.parent/'orbit_certificate.json').read_text());V=c['representatives'];G=c['group_words'];L=c['norm_lower'];n=len(V)
assert d['passed'] and d['products_checked']==n*(n+1)//2*len(G)-n==1683122416
assert d['maximum_integer_dot']<=L*L//2 and d['allowed_integer_dot']==L*L//2
assert d['minimum_squared_norm']>=L*L and d['maximum_squared_norm']<=c['norm_upper']**2
a,b,g=d['witness'];dot=sum(x*y*(1-2*((g>>j)&1)) for j,(x,y) in enumerate(zip(V[a],V[b])));na=sum(x*x for x in V[a]);nb=sum(x*x for x in V[b]);assert dot==d['maximum_integer_dot'] and 4*dot*dot<=na*nb
d.update(certificate_sha256=hashlib.sha256((H.parent/'orbit_certificate.json').read_bytes()).hexdigest(),dense_points=n*len(G),simd_scalar_cross_checks=256*len(G),unbounded_integer_witness_check={'dot':dot,'norm_i':na,'norm_j':nb,'left_4_dot_squared':str(4*dot*dot),'right_norm_product':str(na*nb)})
(H/'dense_result.json').write_text(json.dumps(d,indent=2)+'\n')
raw=(H/'dense_input.bin').read_bytes();controls=[]
for name,code in [('duplicate_representative',6),('zero_norm',5)]:
 with tempfile.TemporaryDirectory(prefix='kissing39_binary_control_') as directory:
  altered=bytearray(raw);rows=np.frombuffer(altered,dtype=np.int32,count=n*40,offset=40).reshape(n,40)
  if name=='duplicate_representative':rows[1]=rows[0]
  else:rows[0]=0
  path=Path(directory)/'bad.bin';path.write_bytes(altered)
  p=subprocess.run([str(H/'dense_direct'),str(path)],env=env,capture_output=True,text=True)
  assert p.returncode==code,(name,p.returncode,p.stderr)
  controls.append({'control':name,'rejected':True,'exit_code':code})
(H/'negative_controls.json').write_text(json.dumps(controls,indent=2)+'\n')
print(run(sys.executable,str(H/'audit_sparse39.py'),'--source',str(H.parent),'--output',str(H)).stdout)
print(run(sys.executable,str(H/'audit_sparse_witnesses39.py'),'--directory',str(H)).stdout)
s=json.loads((H/'sparse_dense_sparse_audit39.json').read_text());w=json.loads((H/'actual_extreme_pair_witnesses39.json').read_text());assert s['status']==w['status']=='PASS';assert s['counts']['sparse_total']==441244 and s['counts']['total_conditional_on_independent_dense_distinctness']==769436
assert s['counts']['dense_total_conditional_on_independent_dense_distinctness']==d['dense_points']
assert s['dense_sparse']['middle']['literal_signed_sums_checked']==8848056320
result={'status':'PASS','dimension':39,'points':769436,'certificate_sha256':d['certificate_sha256'],'dense_relative_products_checked':d['products_checked'],'literal_dense_middle_signed_sums_checked':s['dense_sparse']['middle']['literal_signed_sums_checked'],'dense_maximum_integer_dot':d['maximum_integer_dot'],'middle_maximum_signed_sum':s['dense_sparse']['middle']['maximum_signed_sum'],'middle_exact_squared_slack':s['dense_sparse']['middle']['nonnegative_integer_slack'],'complete_pair_family_coverage':True,'exact_arithmetic_only':True,'producer_code_imported_or_executed':False,'negative_controls_passed':len(controls)}
(H/'complete_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
