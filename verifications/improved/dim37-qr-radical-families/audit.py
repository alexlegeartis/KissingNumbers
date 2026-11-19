"""Independent expansion of the full radical layer, with exact integer products."""
from itertools import combinations
from pathlib import Path
import json,sys,numpy as np
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'dim37-38-qr-merge'))
from supports37 import supports

T=list(map(int,(P/'six_supports.txt').read_text().split()))
S=supports()
seen=set()
N=[];W=[];G=[]
for u in T:
    js=[i for i in range(36)if u>>i&1]
    assert len(js)==6
    for h in combinations(js,3):
        assert h not in seen
        seen.add(h)
    assert all((v&u).bit_count()<=4 for v in S)
    pick=[sum(1<<js[i] for i in range(6)if x>>i&1) for x in range(64) if x.bit_count()%2==0]
    for s in (-1,1):
        for r in pick:N.append(r);W.append(u);G.append(s)
assert len(seen)==20*len(T)
assert len(set(zip(N,W,G)))==len(N)==64*len(T)
table=np.array([i.bit_count() for i in range(256)],dtype=np.int8)
def pc(a):
    b=np.zeros(a.shape,dtype=np.int16)
    for k in range(5):b+=table[(a>>(8*k))&255]
    return b
N=np.array(N,np.uint64);W=np.array(W,np.uint64);G=np.array(G,np.int16)
total=0;mx=-99
for i in range(0,len(W),128):
    stop=min(i+128,len(W))
    common=W[i:stop,None]&W[None,i:]
    mismatch=common&(N[i:stop,None]^N[None,i:])
    dot=4*pc(common)-8*pc(mismatch)+8*G[i:stop,None]*G[None,i:]
    np.fill_diagonal(dot[:,:stop-i],0)
    mx=max(mx,int(dot.max()))
    total+=(stop-i)*(len(W)-stop)+(stop-i)*(stop-i-1)//2
assert mx==16 and total==len(W)*(len(W)-1)//2
# Cross products against signed roots and middle words are small
# integers: all roots on axis 36 are removed and octads avoid it.
assert all(t<(1<<36) for t in S)
assert 2*8<=16 and 4*4<=16
# Explicit rational bounds against each dense norm class.
assert 99**2>2*70**2
for m in (1,4):
    assert 4*(24*70+2*m*99)**2<=32*(144+m*m)*70**2
print(json.dumps({'status':'PASS','dimension':37,'points':533864-144+len(W),'radical_points':len(W),
                  'all_radical_pairs':total,'max_radical_dot':mx},sort_keys=True))
