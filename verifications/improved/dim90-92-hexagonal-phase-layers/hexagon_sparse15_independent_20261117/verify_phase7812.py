from pathlib import Path
from collections import Counter
import hashlib,json,time
P=Path(__file__).resolve().parent
f=P.parent/'hexagon_sparse15_phase7812_20261117/phase7812.txt'
b=f.read_bytes();assert hashlib.sha256(b).hexdigest()=='d5ba3940239da9750271082790a2ac3f62aa52a14c4ee66727f9a66121f593e2'
words=b.decode().splitlines();assert len(words)==7812 and len(set(words))==7812
assert all(len(w)==15 and set(w)<=set('012') for w in words)
planes=[(sum((c=='1')<<j for j,c in enumerate(w)),sum((c=='2')<<j for j,c in enumerate(w))) for w in words]
h=[0]*16
for i,(a,b) in enumerate(planes):
 for c,d in planes[:i]:h[((a^c)|(b^d)).bit_count()]+=1
assert sum(h)==7812*7811//2 and sum(h[:5])==0 and h[5]>0
r={'status':'PASS','independent_algorithm':'two binary indicator planes with integer XOR/OR population counts for every unordered pair','phase_size':7812,'pairs':sum(h),'distance_histogram':h,'minimum':5,'sparse45':130944*7812,'dense90_union':2*3**21+130944*7812,'gain_over_previous_internal90':2*3**21+130944*7812-20961198018,'gain_over_public92':2*3**21+130944*7812-21896439920}
(P/'phase7812_independent.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
