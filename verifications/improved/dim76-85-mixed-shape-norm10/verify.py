from pathlib import Path
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, factorial, prod, ceil, isqrt
import json,sys
import numpy as np
if '--generate' in sys.argv:
 from scipy.optimize import linprog
HERE=Path(__file__).resolve().parent
N=15281788354560; N8=6218175600
P={0:5118989848896,1:7238164864000,2:2503639594500,3:395496128000,4:24994840000,5:500966912,6:2112250,10:2}
assert sum(P.values())==N
assert all(t%2==0 for j,t in P.items() if j)

# Published lineage: Nebe Gamma72 even unimodular/minimum 8; Venkov's 11-design theorem.
def independently_replay_parent_counts():
    # Weight-36 modular basis E4**9, E4**6*Delta, E4**3*Delta**2, Delta**3.
    def mul(a,b):
        return [sum(a[j]*b[i-j] for j in range(i+1)) for i in range(6)]
    def power(a,k):
        w=[1,0,0,0,0,0]
        for _ in range(k):
            w=mul(w,a)
        return w
    e4=[1]+[240*sum(d**3 for d in range(1,n+1) if n%d==0) for n in range(1,6)]
    delta=[0,1,0,0,0,0]
    for n in range(1,6):
        delta=mul(delta,[(-1)**(i//n)*comb(24,i//n) if i%n==0 else 0 for i in range(6)])
    theta=power(e4,9)
    coefficients=[]
    for i in range(1,4):
        a=mul(power(e4,9-3*i),power(delta,i))
        assert a[:i]==[0]*i and a[i]==1
        coefficients.append(-theta[i])
        theta=[s+coefficients[-1]*t for s,t in zip(theta,a)]
    assert coefficients==[-2160,965520,-27302400] and theta==[1,0,0,0,N8,N]
    def populations(pop, squared_norm_product, known):
        # Evaluate Lagrange basis at a squared product; its expectation isolates one bin.
        def moment(r):
            return F(pop*prod(range(1,2*r,2))*squared_norm_product**r,prod(72+2*j for j in range(r)))
        counts={}
        for b in range(6):
            z=[F(1)]
            for a in range(6):
                if b==a:continue
                nz=[F(0)]*(len(z)+1)
                for i,c in enumerate(z):
                    nz[i] -= c*a*a/(b*b-a*a)
                    nz[i+1] += c/(b*b-a*a)
                z=nz
            count=sum(c*(moment(i)-sum(p*j**(2*i)for j,p in known.items()))for i,c in enumerate(z))
            assert count.denominator==1 and count>=0
            counts[b]=int(count)
        return counts|known
    small=populations(N8,80,{})
    # Norm-10 vectors x with |<v,x>|=6 biject with norm-8 u=v-x with <v,u>=4,
    # using one orientation in each case. The remaining exceptional bin is ±v at product10.
    assert populations(N,100,{6:small[4],10:2})==P

# (integer tail shape, squared scale on the Gamma head).
PALETTES={
    74:((4,1),F(26,41)),
    75:((2,2,1),F(7,11)),
    76:((2,1,1),F(8,13)),
    77:((2,2,1,1),F(16,25)),
    78:((2,2,1,1),F(16,25)),
    79:((2,1,1,1,1),F(12,19)),
    80:((2,1,1,1,1),F(12,19)),
    81:((1,1,1,1),F(4,7)),
    82:((1,1,1,1),F(4,7)),
    83:((1,1,1,1),F(4,7)),
    84:((1,1,1,1),F(4,7)),
    85:((1,1,1,1),F(4,7)),
}

def oddfac(n):return prod(range(n,0,-2))
def wick(r,s,j):
 return sum(comb(r,l)*comb(s,l)*factorial(l)*oddfac(r-l-1)*oddfac(s-l-1)*10**((r+s-2*l)//2)*j**l for l in range(min(r,s)+1) if (r-l)%2==0 and(s-l)%2==0)
def tails_of(k,shape):
    # Enumerate exactly all signed permutations of the shape, without duplicates.
    groups=sorted(Counter(shape).items(),reverse=True)
    def absolute(inds,g):
        if g==len(groups):
            yield ()
            return
        val,c=groups[g]
        for used in combinations(inds,c):
            for tail in absolute(tuple(i for i in inds if i not in used),g+1):
                yield tuple((i,val) for i in used)+tail
    out=[]
    for absol in absolute(tuple(range(k)),0):
        for signs in product((-1,1),repeat=len(shape)):
            t=[0]*k
            for (i,v),j in zip(absol,signs):
                t[i]=v*j
            out.append(tuple(t))
    expect=2**len(shape)*factorial(k)//(factorial(k-len(shape))*prod(factorial(c)for _,c in groups))
    assert len(out)==len(set(out))==expect
    return np.array(out,dtype='i8')
def profile_data(j):
 points=[]
 for a,b in product(range(-6,7),repeat=2):
  if 1000+2*j*a*b-10*(j*j+a*a+b*b)<0:continue
  if any(0<10*(x*x+y*y+z*z)+2*(j*x*y+a*x*z+b*y*z)<8
         for x,y,z in product(range(-2,3),repeat=3) if z and (x or y)):continue
  points.append((a,b))
 exc=[(10,j),(-10,-j),(j,10),(-j,-10)]
 exps=[(a,b)for a in range(11)for b in range(11-a)if(a+b)%2==0]
 rhs=[]
 for a,b in exps:
  m=(a+b)//2
  total=F(N*10**m*wick(a,b,j),prod(72+2*i for i in range(m)))
  rhs.append((total-sum(x**a*y**b for x,y in exc))/6**(a+b))
 rhs.extend(F(P[abs(val)]//(1 if val==0 else 2)-sum(pt[axis]==val for pt in exc))
         for axis in range(2)for val in range(-6,7))
 def feat(a,b):
  p=[F(a**r*b**s,6**(r+s))for r,s in exps]
  return p+[F(int((a,b)[axis]==val))for axis in range(2)for val in range(-6,7)]
 return points,exc,rhs,[feat(a,b)for a,b in points]
def low_factor(D,T,blocks=512):
 # Same Boßecker–Rautenbach recurrence as pinned 79–86D proof, but tighter rational blocks.
 # Fixed-point downward rounding is done at every stage rather than storing enormous Fractions.
 lo=isqrt(2*T)+2
 while lo*(lo-1)<=2*T:lo+=1
 assert lo<=D
 # Geometric step from threshold plus fine uniform interior. First few blocks may offer little gain.
 scale=10**32
 g=0
 width=max(1,(D-lo+blocks)//blocks)
 for a in range(lo,D+1,width):
  b=min(D,a+width-1);L=b-a+1;den=a*a+1-2*T
  B=(1-F(a+1,den))/b;C=F(2*(a+T),a*den)
  assert B>=0 and 0<C<1
  nxt=max(F(0),1-L*C)*g + scale*L*B/(1+L*C)
  g=nxt.numerator//nxt.denominator
 return 1+F(g,scale)
def calc(k,generate=False):
 shape,X=PALETTES[72+k]
 ts=tails_of(k,shape);M=len(ts);r=np.array(list(shape)+[0]*(k-len(shape)),dtype='i8')
 norm=sum(t*t for t in shape)
 Q=8-10*X
 assert Q>0 and 25*X<=16 and 4*(shape[0]+shape[1])**2*Q/norm<=16
 assert np.all(np.sum(ts*ts,axis=1)==norm)
 H=X.numerator*norm
 B=8*X.denominator-10*X.numerator
 LIMIT=4*X.denominator*norm
 def edge(j,d):return abs(H*j+B*d)>LIMIT
 dotr=ts@r
 npos=len(shape)
 if shape==(1,1,1,1):
  keys=[(int(np.count_nonzero(t[:4])),int(t[:4].sum()))for t in ts]
 else:
  marks={a:np.flatnonzero(r==a)for a in sorted(set(shape))}
  keys=[tuple(tuple(sorted(map(int,tinds[inds])))for inds in marks.values())+
        (tuple(sorted(map(int,abs(tinds[npos:])))),)for tinds in ts]
 hist=Counter(keys)
 reps={key:t for key,t in zip(keys,ts)}
 dotkeys={key:int(v@r)for key,v in reps.items()}
 assert all(dotkeys[key]==d for key,d in zip(keys,dotr))
 D=int(np.count_nonzero(edge(10,dotr)))-1
 same=D
 for j in range(7):
  assert j>=4 or not any(edge(j,d)for d in dotkeys.values())
  D+=sum(count for key,count in hist.items()if edge(j,dotkeys[key]))*P[j]//2
 triangles2=F(same*D)
 file=HERE/f'cert_{k}.json'
 frozen={} if generate else json.loads(file.read_text())['duals']
 new={}
 for j in [4,5,6]:
  points,exc,rhs,A=profile_data(j)
  sc=[F(N)/abs(x) if x else F(1) for x in rhs]
  if generate:
   af=np.asarray([[float(t*s)for t,s in zip(row,sc)]for row in A]).T
   rf=np.asarray([float(t*s/N)for t,s in zip(rhs,sc)])
  for key,cnt in sorted(hist.items()):
   if not edge(j,dotkeys[key]):continue
   label=f'{j}:{key}'
   dott=ts@reps[key]
   def cost(a,b):return F(int(np.sum(edge(a,dotr)&edge(b,dott))),2)
   c=[cost(a,b) for a,b in points]
   if generate:
    sol=linprog(-np.asarray([float(v)for v in c]),A_eq=af,b_eq=rf,method='highs')
    assert sol.success,sol.message
    y=[F(float(-u)).limit_denominator(10**12)*s for u,s in zip(sol.eqlin.marginals,sc)]
    y[0]+=max(F(0),max(t-sum(a*u for a,u in zip(row,y)) for row,t in zip(A,c)))
    new[label]=[str(t) for t in y]
   else:y=list(map(F,frozen[label]))
   assert len(y)==len(rhs)
   assert all(sum(a*u for a,u in zip(row,y))>=t for row,t in zip(A,c))
   common=sum(y_i*x for y_i,x in zip(y,rhs))+sum(cost(a,b)for a,b in exc)
   triangles2+=common*(P[j]//2)*cnt
 T=ceil(triangles2/2); V=N//2*M
 f=low_factor(D,T)
 pairs=ceil(F(V,D+1)*f)
 out=dict(dimension=72+k,vertices=V,degree=D,triangles=T,ratio=float(F(2*T,D*D)),factor=float(f),pairs=pairs,bound=N8+2*pairs+2*k*(k-1))
 if generate:file.write_text(json.dumps(dict(out,duals=new),indent=2)+'\n')
 else:
  saved=json.loads(file.read_text())
  assert all(saved[key]==v for key,v in out.items())
 return out
if __name__=='__main__':
 assert __debug__
 independently_replay_parent_counts()
 dims=[int(d)for d in sys.argv[1:] if not d.startswith('--')] or list(range(76,86))
 for dim in dims:print(calc(dim-72,'--generate' in sys.argv),flush=True)
