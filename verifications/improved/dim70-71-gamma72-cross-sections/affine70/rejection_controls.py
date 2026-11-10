import json,copy,subprocess,sys,tempfile
from pathlib import Path
p=Path(__file__).parent;src=json.load(open(p/'certificate.json'));mutations={'bound':lambda c:c.update(bound=c['bound']+1),'numerator':lambda c:c.update(numerator=str(int(c['numerator'])+1)),'gram':lambda c:c['gram'][0].__setitem__(1,3),'rhs':lambda c:c['rhs'].__setitem__(0,c['rhs'][0]+1),'row':lambda c:c['rows'][0][0].__setitem__(1,c['rows'][0][0][1]+1),'offset':lambda c:c['offsets'][0].__setitem__(0,8),'cell':lambda c:c['cells'][0].__setitem__(0,8),'dual':lambda c:c['multipliers'].__setitem__(0,str(int(c['multipliers'][0])+1))};out=[]
for name,mutate in mutations.items():
 c=copy.deepcopy(src);mutate(c)
 with tempfile.TemporaryDirectory() as td:
  f=Path(td)/'bad.json';f.write_text(json.dumps(c));r=subprocess.run([sys.executable,str(p/'verify.py'),str(f)],capture_output=True,text=True);assert r.returncode!=0,name;out.append(dict(mutation=name,rejected=True))
(p/'rejection_controls.json').write_text(json.dumps(out,indent=2));print(out)
