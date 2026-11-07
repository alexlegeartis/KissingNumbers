"""Run the complete portable 32D continuation proof and rejection controls."""
from pathlib import Path
import json,os,subprocess,sys
D=Path(__file__).resolve().parent
if not __debug__: raise SystemExit('Run without -O: assertions are required.')
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
def run(script,*args):
 subprocess.run([sys.executable,'-B',str(D/script),*args],cwd=D,env=env,check=True)
for script in ('verify.py','verify_coordinates.py','independent_audit.py','verify_pool_bound.py','verify_delta.py','negative_controls.py'):
 run(script)
a=json.loads((D/'verification.json').read_text())[0]
b=json.loads((D/'coordinate_verification.json').read_text())[0]
c=json.loads((D/'report.json').read_text())
p=json.loads((D/'pool_bound.json').read_text())
n=json.loads((D/'negative_controls.json').read_text())
assert a['all_checks_pass'] and b['all_checks_pass'] and c['status']=='PASS'
assert a['verified_points']==b['verified_points']==c['certificate_audit']['total']==352832
assert a['certificate_sha256']==c['certificate_sha256']=='69c3f368a4aae5956ca75a8569f761224148831e8a85edd5972a0b3b057af646'
assert c['baseline_audit']['total']==352800 and c['delta']['net_gain']==32
assert len(c['negative_controls'])==9 and all(x['rejected'] for x in c['negative_controls'])
assert len(n)==8 and all(x['rejected'] for x in n)
assert p['all_checks_pass'] and p['upper_bound']==p['achieved']==1664
bad=subprocess.run([sys.executable,'-B',str(D/'independent_audit.py'),'--expected-sha256','0'*64],cwd=D,env=env,capture_output=True,text=True)
assert bad.returncode!=0 and 'certificate SHA256' in bad.stderr
print(json.dumps({'status':'PASS','dimension':32,'total':352832,'gain':32,'independent_geometric_controls':9,'hash_control_rejected':True,'producer_control_rejections':8,'fixed_25_support_capacity':1664}))
