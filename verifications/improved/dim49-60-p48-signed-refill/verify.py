"""Read-only complete signed-refill verification with manifest guards."""
from pathlib import Path
import hashlib,importlib.util,json
from verify_single_class import check as check_single
from verify_signed_family import check_base,verify_dimension,PARENT
from final import results
HERE=Path(__file__).resolve().parent

def manifest():
 expected=json.loads((HERE/'SHA256.json').read_text())
 actual={p.relative_to(HERE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in HERE.rglob('*') if p.is_file() and p.name!='SHA256.json' and '__pycache__' not in p.parts}
 assert actual==expected,'Signed-refill manifest mismatch'
 return len(actual)

def main():
 if not __debug__:raise SystemExit('Assertions must be enabled.')
 n=manifest();single=check_single(HERE/'p48_signed_class_lns.npz',cpp=True)
 B,G,pairs=check_base();reports=[]
 for d in range(49,61):
  report=verify_dimension(d,B,G);reports.append(report);print(json.dumps(report),flush=True)
 spec=importlib.util.spec_from_file_location('signed_refill_tail_check',PARENT/'tail_and_lift.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);tails=module.verify()
 counts=results();assert [r['bound'] for r in reports]==[r['total'] for r in counts]
 assert manifest()==n
 print(json.dumps(dict(status='PASS',single_class=single,base_pairs=pairs,dimensions=reports,tails=tails,signed_count='52416000 + sum of signed class sizes + pole count',manifest_files=n),indent=2))

if __name__=='__main__':main()
