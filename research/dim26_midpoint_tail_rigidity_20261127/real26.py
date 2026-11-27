from tail26 import *
g,pinned,n=bank_graph()
G=max(range(n),key=lambda i: int(np.sum(g==i)))
idxG=[np.flatnonzero((g==G)&(pcol==i)) for i in range(2)]
print(json.dumps({"giant":G,"counts":[len(i) for i in idxG],"pinned":G in pinned}),flush=True)
C=candidates(); Z=np.array(list(C),np.int32)
hist=Counter(); okay=[]
mat=P.astype(np.float32).T.copy()
for j in range(0,len(Z),512):
    V=Z[j:j+512]
    d=(V.astype(np.float32)@mat).astype(np.int32)
    keep=d.max(axis=1)<=128
    if not keep.any():continue
    D0,D1=(d[keep][:,k].max(axis=1) for k in idxG)
    for a,b in zip(D0,D1):hist[(int(a),int(b))]+=1
    # Coordinates |A·t| and |B·t| have squared sum 1024. Doubling to d2 scale gives 4096.
    good=(np.minimum(128-D0,64)**2+np.minimum(128-D1,64)**2)>=4096
    okay.extend(V[keep][good].tolist())
(here/"real26_geometric_candidates.json").write_text(json.dumps(okay)+"\n")
print(json.dumps({"hist":[[list(a),b] for a,b in sorted(hist.items())],"geometric_candidates":len(okay)}),flush=True)
