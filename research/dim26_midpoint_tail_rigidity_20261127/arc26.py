from tail26 import *
from math import isqrt
g,pinned,n=bank_graph()
vq=[tuple(q) for q in Q.tolist()]
fl=[tuple(q) for q in FT.tolist()]
rows=[np.flatnonzero((g==j)&(pcol==i)) for j in range(n) for i in range(2)]
variants=json.load(open(here/"real26_geometric_candidates.json"))

def sign_le(a,b,r):
    if b==0 or r==0:return a<=0
    if a<=0 and b<=0:return True
    if a>0 and b>0:return False
    if b>0:return b*b*r<=a*a
    return a*a<=b*b*r

def checker(b2,econ,HH):
    # e=(ax,ay,bx,by,R,d), t=(a+b*sqrt(R))/d
    xa,ya,xb,yb,r,d=map(int,econ)
    H=[int(v) for v in HH]
    def le(v,h):
        x,y=map(int,v)
        return sign_le(2*(x*xa+y*ya)-h*d,2*(x*xb+y*yb),r)
    def ab(v,h):
        return h>=0 and le(v,h) and le((-int(v[0]),-int(v[1])),h)
    # Check all paired components; most donor rays fail first here.
    X=H[-2*n:]
    rotations=[]
    for j in range(n):
        h,k=X[2*j:2*j+2]
        if ab(vq[0],h) and ab(vq[1],k):continue
        if j not in pinned and ab(vq[1],h) and ab(vq[0],k):rotations.append(j);continue
        return None
    for v,h in zip(fl,H): 
        if not le(v,h):return None
    return rotations

visited=0
for b2 in np.asarray(variants,np.int64):
    dots=P@b2
    HH0=128*Fden-FB@b2
    HH1=[int(128-dots[z].max()) if len(z) else 10000 for z in rows]
    HH=np.r_[HH0,HH1]
    cuts=[]
    for q,h in zip(fl,HH0):
        L=q[0]*q[0]+q[1]*q[1]
        if not L:
            if h<0: break
            continue
        if h<0 and h*h>128*L:break
        if h*h<=128*L:cuts.append((q,int(h)))
    else:
        cuts.extend((q,int(h)) for h in HH1 if h<=64 for u in vq for q in (u,(-u[0],-u[1])))
        raw=[(4,4,0,0,0,1)]
        for (u,v),h in set(cuts):
            L=u*u+v*v
            R=128*L-h*h
            if R<0:continue
            raw.extend((h*u,h*v,s*(-v),s*u,R,2*L) for s in (-1,1))
        for e in raw:
            visited+=1
            rr=checker(b2,e,HH)
            if rr is not None:
                ans=dict(b2=b2.tolist(),tail_representation=list(e),flipped_components=rr)
                print(json.dumps(ans),flush=True)
                (here/"arc26_winner.json").write_text(json.dumps(ans)+"\n")
                raise SystemExit
print(json.dumps(dict(bodies=len(variants),exact_angles_tested=visited,gain=0)),flush=True)
