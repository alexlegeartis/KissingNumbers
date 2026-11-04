exec(open('research/sol_discovery_20261104/orbit_probe.py').read().split('for n in [')[0])
def partitions(n,m=2):
 if n==0:yield []
 for i in range(m,n+1):
  for p in partitions(n-i,i):yield [i]+p
n=12
for cycles in partitions(n):
 perm=[];at=0
 for s in cycles:perm.extend(list(range(at+1,at+s))+[at]);at+=s
 remaining=set(itertools.combinations(range(n),4));orbits=[];invoices=[]
 while remaining:
  b=min(remaining);orb=set();cur=b
  while cur not in orb:orb.add(cur);cur=tuple(sorted(perm[i] for i in cur))
  remaining.difference_update(orb);ts=[t for s in orb for t in itertools.combinations(s,3)]
  if len(ts)==len(set(ts)):orbits.append(sorted(orb));invoices.append(ts)
 if not orbits:continue
 ti={t:i for i,t in enumerate(itertools.combinations(range(n),3))};rr=[];cc=[]
 for j,ts in enumerate(invoices):
  for t in ts:rr.append(ti[t]);cc.append(j)
 M=coo_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(ti),len(orbits))).tocsc()
 r=milp(-np.array(list(map(len,orbits))),integrality=np.ones(len(orbits)),bounds=Bounds(0,1),constraints=LinearConstraint(M,0,1),options={'time_limit':.3})
 selected=[] if r.x is None else np.flatnonzero(r.x>.5);count=sum(len(orbits[i]) for i in selected)
 print(json.dumps({'cycles':cycles,'supports':count,'status':int(r.status)}),flush=True)
 if count>=51:break
