"""Distance tables and bounds for PAIRWISE-DISJOINT repeated-mask translate unions.
Zero distance proves overlap, not invalidity of the deduplicated union.
"""
from pathlib import Path
import json,time,hashlib
from collections import Counter
import numpy as np
OUT=Path(__file__).resolve().parent
start=time.monotonic();I=np.arange(16);T=np.array([sum((((a>>3)&1)^sum(((a>>j)&1)*((x>>j)&1)for j in range(3))%2)<<x for x in range(8))for a in range(16)])
local=np.array([i.bit_count()for i in range(256)]);pairs=T[:,None]^T[None,:]
D=np.zeros((16,256),dtype=int)
for m in range(256):
 c=local[pairs^m]
 two=np.array([[np.min(c+c[np.ix_((s-I)%16,(t-I)%16)])for t in range(16)]for s in range(16)])
 three=np.array([[np.min(two+c[np.ix_((s-I)%16,(t-I)%16)])for t in range(16)]for s in range(16)])
 for scalar in range(16):D[scalar,m]=np.min(three+c[np.ix_(scalar*I%16,scalar*I%16)])
results=[]
for scalar in range(16):
 good=[m for m in range(256)if D[scalar,m]>=8]
 bad_even=[m for m in range(256)if m.bit_count()%2==0 and D[scalar,m]<8]
 even_good=all(m.bit_count()%2==0 for m in good)
 closed=all(a^b in bad_even for a in bad_even for b in bad_even)
 incompatible=set(m for m in range(256)if D[scalar,m]<8);subspace=[0]
 for v in range(1,256):
  if v not in subspace and all(v^h in incompatible for h in subspace):subspace += [v^h for h in subspace]
 bound=256//len(subspace)
 results.append(dict(scalar=scalar,cross_distances=D[scalar].tolist(),cross_distance_histogram=dict(Counter(map(int,D[scalar]))),all_compatible_differences_even=even_good,bad_even_subspace=bad_even if closed else None,incompatible_subspace=subspace,max_translates_upper_bound=bound,disjoint_translate_point_count_upper_bound=4096*bound if bound else None,bound_scope="pairwise-disjoint translates only"))
 print(scalar,'bad_even',len(bad_even),'even',even_good,'closed',closed,'max translates',bound,flush=True)
receipt=dict(description='All256 repeated8bit masks; distance includes zero for overlap. Coset bounds apply ONLY to pairwise-disjoint translate unions. Prior unrestricted upper bound retracted.',bound_scope='pairwise-disjoint translates only',unrestricted_union_status='unresolved; prior unrestricted upper bound retracted',results=results,elapsed_seconds=time.monotonic()-start,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(OUT/'repeated_mask_dp.json').write_text(json.dumps(receipt,indent=2)+'\n');print('elapsed',receipt['elapsed_seconds'])
