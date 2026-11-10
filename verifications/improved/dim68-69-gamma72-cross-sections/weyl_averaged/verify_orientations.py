import json,itertools,sys
from pathlib import Path
p=Path(__file__).parent;K=int(sys.argv[1]);c=json.loads(Path(sys.argv[2] if len(sys.argv)>2 else p/f'certificate{K}.json').read_text());G=c['gram'];assert __debug__
I=tuple(int(i==j) for i in range(K) for j in range(K))
def mul(A,B):return tuple(sum(A[i*K+k]*B[k*K+j] for k in range(K)) for i in range(K) for j in range(K))
def tr(A):return tuple(A[j*K+i] for i in range(K) for j in range(K))
gen=[]
# A rows describe new normal vectors in the old normal basis.
for k in range(K):
 A=list(I)
 for j in range(K):A[j*K+k]-=G[j][k]//4
 gen.append(tuple(A))
for a,b in ([(1,2),(2,3)] if K==4 else []):
 A=[list(I[i*K:(i+1)*K]) for i in range(K)];A[a],A[b]=A[b],A[a];gen.append(tuple(sum(A,[])))
gen.append(tuple(-v for v in I))
flat=tuple(sum(G,[]))
assert all(mul(mul(A,flat),tr(A))==flat for A in gen)
seen={I};todo=[I]
while todo:
 A=todo.pop()
 for S in gen:
  B=mul(S,A)
  if B not in seen:seen.add(B);todo.append(B)
assert len(seen)=={3:48,4:1152}[K]
assert all(mul(mul(A,flat),tr(A))==flat for A in seen)
# Every orientation is an integral basis change, so its norm-8 rows are actual lattice vectors.
# Direct action-only orbit graph (no population translation identities).
cells=list(map(tuple,c['cells']));idx={v:i for i,v in enumerate(cells)};remaining=set(range(len(cells)));orbits=[]
while remaining:
 i=min(remaining);v=cells[i]
 orb={idx[tuple(sum(A[j*K+k]*v[k] for k in range(K)) for j in range(K))] for A in seen}
 orbits.append(sorted(orb));remaining-=orb
assert orbits==c['groups']
R={'orientation_group_order':len(seen),'orbit_count':len(orbits),'orbit_sizes':list(map(len,orbits)),'all_generators_integral_gram_isometries':True,'all_groups_are_actual_orientation_orbits':True,'shell_invariance_required':False}
(p/f'orientation{K}.json').write_text(json.dumps(R,indent=2));print(json.dumps(R,indent=2))
