import itertools,json
from pathlib import Path
D=Path(__file__).resolve().parent
P0=[tuple(range(3*i,3*i+3))for i in range(4)];S0=set(tuple(sorted(x))for x in itertools.product(*P0));outside=[s for s in itertools.combinations(range(12),4)if len({x//3 for x in s})<=2];assert len(outside)==90;idx={s:i for i,s in enumerate(outside)}
def parts(rem):
 if not rem:yield[];return
 a=min(rem)
 for b,c in itertools.combinations(sorted(rem-{a}),2):
  g=(a,b,c)
  for tail in parts(rem-set(g)):yield[g]+tail
masks={};hist={}
for P in parts(set(range(12))):
 A=tuple(sorted(idx[s]for x in itertools.product(*P)if(s:=tuple(sorted(x)))in idx));masks.setdefault(A,P);hist[len(A)]=hist.get(len(A),0)+1
items=sorted(masks.items(),key=lambda x:-len(x[0]));perms=list(itertools.permutations(range(4)))
def canon(P):
 mat=[[len(set(a)&set(b))for b in P]for a in P0]
 return min(tuple(sorted(tuple(row[j]for j in pc)for row in mat))for pc in perms)
reps={}
for A,P in items:reps.setdefault(canon(P),(A,P))
reps=sorted(reps.values(),key=lambda x:-len(x[0]));print('PARTITIONS',sum(hist.values()),'UNIQUE',len(items),'ORBITS',len(reps),'HIST',hist,flush=True)
conf=[[len(set(a)&set(b))>=3 for b in outside]for a in outside]
def indep(A,B):
 match={}
 def dfs(a,seen):
  for b in B:
   if conf[a][b]and b not in seen:
    seen.add(b)
    if b not in match or dfs(match[b],seen):match[b]=a;return True
  return False
 for a in A:dfs(a,set())
 LA=set(A)-set(match.values());RB=set();queue=list(LA)
 for a in queue:
  for b in B:
   if conf[a][b]and b not in RB:
    RB.add(b)
    if b in match and match[b]not in LA:LA.add(match[b]);queue.append(match[b])
 return sorted(LA),sorted(set(B)-RB),len(A)+len(B)-len(match)
best=30;tested=0
for A,Pa in reps:
 if len(A)+len(items[0][0])<=best:continue
 for B,Pb in items:
  if len(A)+len(B)<=best:break
  AA,BB,n=indep(A,B);tested+=1
  if n>best:
   best=n;out={'dimension':92,'supports_count':81+n,'dense_count':(81+n)*16*8**8,'partitions':[P0,Pa,Pb],'selected_supports':[list(sorted(S0)),[outside[k]for k in AA],[outside[k]for k in BB]],'tested_pairs':tested};(D/'candidate_stronger.json').write_text(json.dumps(out,indent=2)+'\n');print('BEST',81+n,[81,len(AA),len(BB)],tested,flush=True)
print('FINAL',81+best,'tested',tested,flush=True)
(D/'stronger_search_report.json').write_text(json.dumps({'fixed_first_color_supports':81,'maximum_total_supports':81+best,'pair_tests':tested,'partition_count':sum(hist.values()),'eligible_histogram':hist,'symmetry_orbits':len(reps),'scope':'All balanced3+3+3+3 extra partitions, fixed complete first-color palette; computational search report not an independent optimality audit.'},indent=2)+'\n')
