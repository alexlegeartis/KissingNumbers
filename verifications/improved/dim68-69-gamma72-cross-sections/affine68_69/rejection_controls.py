import json,copy,subprocess,sys,tempfile
from pathlib import Path
p=Path(__file__).parent;out=[]
for K in [3,4]:
 src=json.load(open(p/f'certificate{K}.json'))
 mutations={'bound':lambda c:c.update(bound=c['bound']+1),'numerator':lambda c:c.update(numerator=str(int(c['numerator'])+1)),'gram':lambda c:c['gram'][0].__setitem__(1,-3),'group':lambda c:c['groups'][0].pop(),'moment':lambda c:c['tags'][0][1].__setitem__(0,11),'offset':lambda c:c['offsets'][0].__setitem__(0,8),'objective':lambda c:c['objective'].__setitem__(0,c['objective'][0]+1),'dual':lambda c:c['multipliers'].__setitem__(0,str(int(c['multipliers'][0])+1))}
 for name,mutate in mutations.items():
  c=copy.deepcopy(src);mutate(c)
  with tempfile.TemporaryDirectory() as td:
   f=Path(td)/'bad.json';f.write_text(json.dumps(c));r=subprocess.run([sys.executable,str(p/'verify.py'),str(K),str(f)],capture_output=True,text=True);assert r.returncode!=0,(K,name);out.append(dict(dimension=72-K,mutation=name,rejected=True));print(out[-1],flush=True)
(p/'rejection_controls.json').write_text(json.dumps(out,indent=2))
