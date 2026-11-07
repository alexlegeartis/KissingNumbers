import copy,importlib.util,json,tempfile
from pathlib import Path
D=Path(__file__).resolve().parent
mods=[]
for name in ['verify','verify_coordinates']:
 spec=importlib.util.spec_from_file_location(name,D/(name+'.py'));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);mods.append((name,m))
base=json.loads((D/'certificate32_352832.json').read_text());cases=[]
x=copy.deepcopy(base);x['points']+=1;cases.append(('false_count',x))
x=copy.deepcopy(base);x['exceptional_antipodal_lines'][0]=x['exceptional_antipodal_lines'][1];cases.append(('duplicate_exceptional_line',x))
x=copy.deepcopy(base);m,s=x['exceptional_antipodal_lines'][0];b=(m&-m);b=(m^b)&-(m^b);x['exceptional_antipodal_lines'][0]=[m,s^b];cases.append(('incompatible_sign_flip',x))
x=copy.deepcopy(base);x['dense_generator'][0]=x['dense_generator'][1];cases.append(('rank_deficient_dense_generator',x))
results=[]
with tempfile.TemporaryDirectory(dir=D) as td:
 for label,c in cases:
  p=Path(td)/(label+'.json');p.write_text(json.dumps(c))
  for name,m in mods:
   try:m.verify(p)
   except Exception as e:results.append({'case':label,'checker':name,'rejected':True,'exception':type(e).__name__})
   else:raise RuntimeError(f'{label} unexpectedly accepted by {name}')
print(json.dumps(results));(D/'negative_controls.json').write_text(json.dumps(results,indent=2)+'\n')
