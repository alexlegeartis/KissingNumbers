"""Semantic corruption tests for the exact averaged proof checker."""
import json,sys,subprocess,tempfile,copy
from pathlib import Path
p=Path(__file__).parent
results=[]
for k in (3,4):
 original=json.loads((p/f'certificate{k}.json').read_text())
 for name in ('bound','numerator','gram','group','moment','offset','objective','dual','population','dimension'):
  bad=copy.deepcopy(original)
  if name=='bound':bad['bound']+=1
  elif name=='numerator':bad['numerator']=str(int(bad['numerator'])+1)
  elif name=='gram':bad['gram'][0][0]+=1
  elif name=='group':bad['groups'][0].append(bad['groups'][1][0])
  elif name=='moment':bad['tags'][0][1][0]+=2
  elif name=='offset':bad['offsets'][0]=[8]*k
  elif name=='objective':bad['objective'][0]+=1
  elif name=='dual':bad['multipliers'][0]=str(int(bad['multipliers'][0])+1)
  elif name=='population':bad['population']+=1
  elif name=='dimension':bad['dimension']+=1
  with tempfile.TemporaryDirectory() as tmp:
   f=Path(tmp)/'bad.json';f.write_text(json.dumps(bad));run=subprocess.run([sys.executable,'-B',str(p/'verify.py'),str(k),str(f)],capture_output=True,text=True)
   assert run.returncode!=0,(k,name,'Corrupt certificate accepted')
  results.append({'dimension':72-k,'corruption':name,'result':'rejected'})
report={'status':'PASS','rejections':len(results),'cases':results};(p/'negative_controls.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
