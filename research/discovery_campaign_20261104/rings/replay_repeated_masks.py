"""Independent replay of a bound for PAIRWISE-DISJOINT translate unions only."""
from pathlib import Path
import json,time,hashlib
OUT=Path(__file__).resolve().parent
start=time.monotonic()
def require(ok,msg):
 if not ok:raise ValueError(msg)
def gray(a):
 return sum((((a>>3)&1)^((a&1)*(x&1))^(((a>>1)&1)*((x>>1)&1))^(((a>>2)&1)*((x>>2)&1)))<<x for x in range(8))
T=[gray(a)for a in range(16)]
for a in range(16):
 for b in range(16):require((T[a]^T[b]).bit_count()==(0 if a==b else 8 if (a-b)%16==8 else 4),'Gray metric')
# The group is Z16 x Z16: two ring messages, one from each codeword.
subtraction=[[(((s//16-u//16)%16)*16+(s%16-u%16)%16)for u in range(256)]for s in range(256)]
last=[[((scalar*(s//16)%16)*16+scalar*(s%16)%16)for s in range(256)]for scalar in range(16)]
actual=[[0]*256 for _ in range(16)]
for mask in range(256):
 one=[(T[a]^T[b]^mask).bit_count()for a in range(16)for b in range(16)]
 two=[min(one[u]+one[subtraction[s][u]]for u in range(256))for s in range(256)]
 three=[min(two[u]+one[subtraction[s][u]]for u in range(256))for s in range(256)]
 for scalar in range(16):actual[scalar][mask]=min(three[s]+one[last[scalar][s]]for s in range(256))
data=json.loads((OUT/'repeated_mask_dp.json').read_text());counts=[]
require(data.get('bound_scope')=='pairwise-disjoint translates only','unsupported overlapping-union scope')
require(data.get('unrestricted_union_status')=='unresolved; prior unrestricted upper bound retracted','unrestricted branch must remain unresolved')
require([r['scalar']for r in data['results']]==list(range(16)),'scalar coverage')
for r in data['results']:
 require(r.get('bound_scope')=='pairwise-disjoint translates only','row scope')
 scalar=r['scalar'];require(r['cross_distances']==actual[scalar],'distance table')
 H=set(r['incompatible_subspace']);require(0 in H and all(a^b in H for a in H for b in H),'not a subspace')
 require(all(actual[scalar][m]<8 for m in H),'compatible difference inside clique part')
 bound=256//len(H);require(r['max_translates_upper_bound']==bound and r['disjoint_translate_point_count_upper_bound']==4096*bound,'wrong coset bound')
 require(4096*bound<352800,'no competitive exclusion');counts.append(4096*bound)
result=dict(status='PASS_FOR_DISJOINT_TRANSLATES_ONLY',unrestricted_union_status='unresolved',scalar_extensions=16,repeated_masks=256,largest_disjoint_translate_point_upper_bound=max(counts),comparison32=352800,scope='Pairwise-disjoint repeated8bit-mask translate unions only. General unions after point deduplication remain unresolved; distance0 does not establish a positive-distance failure.',elapsed_seconds=time.monotonic()-start,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(OUT/'repeated_mask_independent_replay.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
