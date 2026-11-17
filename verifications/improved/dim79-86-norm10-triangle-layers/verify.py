"""Full exact proof replay, with explicit published lattice/design premises."""
from pathlib import Path
import hashlib,json,runpy,subprocess,sys
if not __debug__:raise RuntimeError('Assertions required')
p=Path(__file__).resolve().parent
for name,digest in json.loads((p/'MANIFEST.json').read_text()).items():
 assert hashlib.sha256((p/name).read_bytes()).hexdigest()==digest,('hash',name)
a=runpy.run_path(str(p/'inherited_source.py'));theta,c=a['theta_coefficients']();n8,n10=theta[4:6]
assert (n8,n10)==(6218175600,15281788354560)
small=a['solve_moments'](n8,80,{});shell=a['solve_moments'](n10,100,{6:small[4],10:2})
assert shell=={0:5118989848896,1:7238164864000,2:2503639594500,3:395496128000,4:24994840000,5:500966912,6:2112250,10:2}
claims=json.loads((p/'claims.json').read_text())
for n in range(79,87):
 subprocess.run([sys.executable,str(p/'verify_one.py'),str(n),'_marginal'],check=True)
 subprocess.run([sys.executable,str(p/'check_recurrence.py'),str(n)],check=True)
 a=json.loads((p/f'audit_result{n}_marginal.json').read_text());assert a['kissing_bound']==claims[str(n)]
print('COMPLETE 79-86 TRIANGLE LAYER PROOF PASS; published Gamma72 and shell-design premises retained')
