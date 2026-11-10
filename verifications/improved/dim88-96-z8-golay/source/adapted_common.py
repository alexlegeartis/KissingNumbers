"""Structural guards used by production replays and semantic rejection tests."""
def check_pattern(key,pattern,equal):
 sizes=list(map(int,key.split(',')))
 if equal:
  assert len(sizes)==5 and sum(sizes) in (3,4)
  weights=[w for w,n in enumerate(sizes) for _ in range(n)]
  if len({w%2 for w in weights})>1:
   assert pattern['count']==0
   return None
 else:
  assert len(sizes)==2
  weights=[2]*sizes[0]+[4]*sizes[1]
 groups=pattern['positions'];flat=sum(groups,[])
 assert [len(g) for g in groups]==weights
 assert all(type(x)==int and 0<=x<24 for x in flat) and len(set(flat))==len(flat)
 assert 0<=pattern['count']<=16384
 return groups

def place_signs(cells,selected,groups,equal):
 perrow={t:[i for i in cells if i//4==t] for t in selected}
 rows=sorted((t for t in selected if equal or perrow[t]),key=lambda t:(len(perrow[t]),t))
 assert len(rows)==len(groups)
 mapping={}
 for t,g in zip(rows,groups):
  assert len(perrow[t])==len(g)
  mapping.update(zip(perrow[t],g))
 mapping.update(zip(sorted(set(cells)-set(mapping)),sorted(set(range(24))-set(mapping.values()))))
 assert set(mapping)==set(cells) and sorted(mapping.values())==list(range(24))
 assert [[mapping[c] for c in perrow[t]] for t in rows]==groups
 return mapping
