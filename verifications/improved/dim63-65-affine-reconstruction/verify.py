from pathlib import Path
import json,hashlib,sys
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2];P=ROOT/'verifications/improved/dim61-67-compatible-sign-families';SOURCE=ROOT/'verifications/improved/dim61-68-95-96-code-layers/source'
sys.path.insert(0,str(P));from subcode67 import code,expand,span,need

def verify():
 c=json.loads((HERE/'certificate.json').read_text());f=json.loads((P/'data/global64.json').read_text())
 with np.load(SOURCE/'dim61_96/code64_supports.npz',allow_pickle=False) as z:parent=list(map(int,z['words']))
 fixed_old={r['support']:code(r,True) for r in f['old']};fixed_new=[code(r,True) for r in f['new']]
 modified={r['support']:code(r,True) for r in c['modified_parent_classes']};added=[code(r) for r in c['new_classes']]
 need(len(parent)==len(set(parent))==30828,'parent count');need(len(fixed_new)==1940,'base added count')
 need(set(modified)<=set(parent) and not set(modified)&set(fixed_old),'modified identity')
 need(len({a['s'] for a in added})==len(added),'new unique');need(not {a['s'] for a in added}&(set(parent)|{a['s'] for a in fixed_new}),'new absent')
 old=fixed_old|modified;checks=0;critical=[]
 for s in modified:
  need(all(s==t or (s&t).bit_count()<=8 for t in parent),'modified old contacts')
  need(all((s&b['s']).bit_count()<=8 for b in fixed_new),'modified fixed new contacts')
 def pair(a,b):
  nonlocal checks
  checks+=1;m=(a['s']&b['s']).bit_count()
  if m<=8:return
  need(m<=10 and 'dual' in b,'overlap repair')
  hs=[h for h in a['dual']&b['dual'] if h.bit_count()==8 and h&(a['s']&b['s'])==h and ((a['phase']^b['phase'])&h).bit_count()%2]
  need(bool(hs),'opposite check');critical.append([a['s'],b['s'],min(hs)])
 for i,a in enumerate(added):
  for s in parent:pair(a,old.get(s,dict(s=s)))
  for b in fixed_new:pair(a,b)
  for b in added[:i]:pair(a,b)
 base=(1<<28)+(30828+1940)*2048+10416*16+128
 gain=sum(a['count'] for a in added)
 need(base==335711104 and base+gain==c['total'],'count');need(gain==38592,'recovered gain')
 out=dict(status='PASS',base=base,gain=gain,total64=base+gain,total65=base+gain+68226,new_classes=len(added),modified_classes=len(modified),support_pairs_checked=checks,critical_repairs=critical,certificate_sha256=hashlib.sha256((HERE/'certificate.json').read_bytes()).hexdigest(),scope='Fresh exact finite extension. Unchanged parent code and owner lift inherited from separately replayed base package.')
 (HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items() if k!='critical_repairs'})
if __name__=='__main__':verify()
