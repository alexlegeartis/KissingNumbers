"""Exact sufficient middle-layer bound for specified weight-eight supports."""
import json,numpy as np
from pathlib import Path

def middle_bound(source,V,L,group):
 source=Path(source);base=json.loads((source/'base_certificate.json').read_text())
 supports=[]
 for line in (source/'source/supports_3424.txt').read_text().splitlines():
  if len(line)!=39 or set(line)-set('01'):raise ValueError('Invalid support row')
  supports.append([i for i,v in enumerate(line) if v=='1'])
 for edit in base['support_replacements']:
  if supports[edit['remove_index']]!=edit['remove']:raise ValueError('Support edit mismatch')
  supports[edit['remove_index']]=edit['add']
 if len(supports)!=3424 or any(len(s)!=8 or len(set(s))!=8 or min(s)<0 or max(s)>=39 for s in supports):raise ValueError('Support shape')
 S=np.asarray(supports,dtype=np.int64)
 even=np.zeros(len(S),dtype=bool) if group is None else np.asarray([all((sum(1<<int(j) for j in support)&int(g)).bit_count()%2==0 for g in group) for support in S])
 minimum=2*L*L;maximum=0;checked=0
 if 2*(8*int(np.abs(V).max()))**2>=2**63:raise ValueError('Arithmetic safety')
 for start in range(0,len(V),64):
  rows=V[start:start+64][:,S];magnitude=np.abs(rows);values=magnitude.sum(axis=2)
  odd=((rows<0).sum(axis=2)%2==1)&even[None,:]
  values-=2*magnitude.min(axis=2)*odd
  slack=2*L*L-values*values
  if int(slack.min())<0:raise ValueError('Middle support angle')
  minimum=min(minimum,int(slack.min()));maximum=max(maximum,int(values.max()));checked+=values.size
 return dict(minimum_middle_support_slack=minimum,maximum_middle_absolute_sum=maximum,middle_support_representative_checks=checked,middle_points_removed=0)
