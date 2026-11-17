import json,copy
from pathlib import Path
p=Path(__file__).resolve().parent
def check(a):
 assert len({tuple(z['h'])for z in a})==len(a)
 for i,z in enumerate(a):
  h=z['h'];r=128-sum(x*x for x in h);T=3*r+z['u']**2+3*z['v']**2
  assert r>0 and 0<T<=384
  for w in a[:i]:
   k=w['h'];S=3*(128-sum(x*x for x in k))+w['u']**2+3*w['v']**2
   P=3*(64-sum(x*y for x,y in zip(h,k)))+z['u']*w['u']+3*z['v']*w['v']
   assert P<=0 or 4*P*P<=T*S
 assert sum(z['count']for z in a)==1361134720
 receipt=json.load(open(p/'independent_result.json'))['rows'];d={tuple(z['h']):z['replayed']for z in receipt}
 assert all(d[tuple(z['h'])]>=z['count']for z in a)
a=json.load(open(p/'candidate.json'));check(a)
cases=[]
b=copy.deepcopy(a);b[0]['u']=100;cases.append(b)
b=copy.deepcopy(a);b.append(b[0]);cases.append(b)
b=copy.deepcopy(a);b[0]['count']+=1;cases.append(b)
b=copy.deepcopy(a);b[0]['h']=[100,0,0,0];cases.append(b)
for b in cases:
 try:check(b)
 except (AssertionError,KeyError):pass
 else:raise AssertionError('Corruption accepted')
print('PASS exact geometry, population receipt binding, four corrupted inputs rejected')
