from pathlib import Path
import json,numpy as np
ROOT=Path(__file__).resolve().parents[3];P=ROOT/'verifications/improved/dim61-67-compatible-sign-families';C=ROOT/'research/dim32plus_continue_20261113/codes';OUT=Path(__file__).resolve().parent
def span(rows):
 a=[0]
 for h in rows:a += [x^h for x in a]
 return a
def expand(s,w):
 return sum(((w>>j)&1)<<i for j,i in enumerate(i for i in range(64) if s>>i&1))
f=json.loads((P/'data/global64.json').read_text());print(f.keys())
old=list(map(int,np.fromfile(C/'old_supports.bin',dtype='<u8')));data={r['support']:r for r in f['old']};new=f['new'];lines=[f'{len(old)} {len(new)} 0']
for s in old+[r['support'] for r in new]:
 r=data.get(s) or next((r for r in new if r['support']==s),None)
 if r:phase=expand(s,r['shift']);checks=[expand(s,h) for h in span(r['checks']) if h.bit_count()==8]
 else:phase=0;checks=[]
 lines.append(' '.join(map(str,[s,phase,len(checks),*checks])))
(OUT/'input.txt').write_text('\n'.join(lines)+'\n');(OUT/'context.json').write_text(json.dumps(dict(old_order=old,extra_full=[],published_new=[r['support'] for r in new])))
