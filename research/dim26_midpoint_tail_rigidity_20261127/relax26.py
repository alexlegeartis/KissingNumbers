from pin_angles import *
from math import isqrt
from fractions import Fraction
assert all(checks)
bank,pinned,G=bank_graph()
giant=max(range(G),key=lambda a:int((bank==a).sum()))
assert giant in pinned
force=set(pinned)|set(bank[~hard[g]])
print(json.dumps(dict(forced=sorted(map(int,force)),giant=int(giant),soft_at=Counter(map(int,bank[~hard[g]])))),flush=True)
allowed=json.load(open(here/'real_donor_survivors.json'))
fl=[tuple(v) for v in FT.tolist()]
vq=[tuple(v) for v in Q.tolist()]
rows=[np.where(hard[g] & (bank==j) & (pcol==i))[0] for j in range(G) for i in range(2)]
soft=[np.where(~hard[g] & (pcol==i))[0] for i in range(2)]
def le(a,b,r):
    if b==0 or r==0:return a<=0
    if a<=0 and b<=0:return True
    if a>0 and b>0:return False
    if b>0:return b*b*r<=a*a
    return a*a<=b*b*r
def test(b2,econ,HH0,HH1,S):
    x,y,u,v,rad,den=econ
    def cmp(q,h):
        a,b=q
        return le(2*(a*x+b*y)-h*den,2*(a*u+b*v),rad)
    def ab(q,h):
        return h>=0 and cmp(q,h) and cmp((-q[0],-q[1]),h)
    def assigned():
        for j in range(G):
            h,k=HH1[2*j:2*j+2]
            if ab(vq[0],h) and ab(vq[1],k): continue
            if j not in force and ab(vq[1],h) and ab(vq[0],k):continue
            return False
        return True
    if not ab(vq[0],S[0]) or not ab(vq[1],S[1]): return False
    return assigned() and all(cmp(q,h) for q,h in zip(fl,HH0))
surv=[];nc=0
for b2 in np.array(allowed,np.int64):
    d=P@b2
    H0=list(map(int,128*Fden-FB@b2))
    H1=[int(128-d[z].max()) if len(z) else 10000 for z in rows]
    SH=[int(136-d[z].max()) if len(z) else 10000 for z in soft]
    cp=[(q,h) for q,h in zip(fl,H0) if q!=(0,0)]
    cp.extend((q,h) for h in H1+SH if h<=64 for u in vq for q in (u,(-u[0],-u[1])))
    E=[(4,4,0,0,0,1)]
    for (u,v),h in set(cp):
        L=u*u+v*v;r=128*L-h*h
        if r>=0:E.extend((h*u,h*v,s*(-v),s*u,r,2*L) for s in (-1,1))
    good=[]
    for e in E:
        nc+=1
        if test(b2,e,H0,H1,SH):good.append(e)
    if good:surv.append((b2.tolist(),good))
print(json.dumps(dict(tested=nc,relaxed_survivors=len(surv))),flush=True)
(here/'relax26_survivors.json').write_text(json.dumps(surv)+'\n')
