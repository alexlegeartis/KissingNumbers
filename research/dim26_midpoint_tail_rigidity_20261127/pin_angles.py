from component_lib import *
from fractions import Fraction
from collections import Counter
from math import gcd,isqrt
caps=64*Fden[None,:]-P@FB.T
zint=np.arange(len(FT))<len(R)
H=np.array([caps[:,zint&(FT[:,i]!=0)&(FT[:,1-i]==0)].min(axis=1) for i in range(2)]).T
H[pcol==1]=H[pcol==1][:,::-1]
Hmin=np.full((n,2),127,np.int64);np.minimum.at(Hmin,g,H)
hard=(Hmin==32).all(axis=1)
assert hard.sum()==954
qs=np.vstack((FT,np.c_[-FT[:,1],FT[:,0]]))
checks=[]
for j in np.flatnonzero(~hard):
 hs=np.r_[64*Fden-stat[2*j],64*Fden-stat[2*j+1]]
 have=set()
 for (x,y),hc in zip(qs.tolist(),hs):
  alpha,beta=x+y,y-x
  if alpha<0:alpha,beta=-alpha,-beta
  if alpha<0 or beta==0:continue
  h=int(hc)
  assert h>=0
  a,b=alpha,abs(beta)
  if 4*b>h and (4*(8*a+b))**2>65*h*h:
   if alpha==0:have.update((-1,1))
   else:have.add((1 if beta>0 else -1))
 found=len(have)==2
 checks.append(found)
assert all(checks)
print('SOFT_ANGLE_PIN',Counter(checks),flush=True)
if __name__=="__main__":
    C=candidates(); Z=np.array(list(C),np.int32)
    mat=P.astype(np.float32).T.copy()
    bank,_,_=bank_graph()
    qA=np.flatnonzero(hard[g]&(pcol==0)&(bank==bank[np.argmax(np.bincount(bank)[bank])]))
    qB=np.flatnonzero(hard[g]&(pcol==1)&(bank==bank[np.argmax(np.bincount(bank)[bank])]))
    sc=Counter();dotsset=set();surv=[]
    for k in range(0,len(Z),512):
     V=Z[k:k+512]
     dz=(V.astype(np.float32)@mat).astype(np.int32)
     ok=dz.max(axis=1)<=128
     u=V[ok];d=dz[ok]
     if not len(u):continue
     nostr=d[:,hard[g]].max(axis=1)<128
     sc['no_hard_exact']+=int(nostr.sum())
     dd=d[nostr]
     if not len(dd):continue
     dotsset.update(np.unique(dd[dd>64]).tolist())
     ss=(128-dd[:,qA].max(axis=1))**2+(128-dd[:,qB].max(axis=1))**2>=4096
     sc['passes_giant_hard_caps']+=int(ss.sum())
     surv.extend(u[nostr][ss].tolist())
    print('SURV',sc,'critical dots2',sorted(dotsset),flush=True)
    (here/'real_donor_survivors.json').write_text(json.dumps(surv)+'\n')
