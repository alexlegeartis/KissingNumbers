"""Exact local exclusion ball around the 842-point near-miss.

For nonnegative pair weights summing to one, let g be their average dot
product, a_i=sum_j lambda_ij x_j, h_i=a_i.x_i, f_i=a_i-h_i*x_i.
If every unit y_i is at chord distance <=r from unit x_i, then

 max y_i.y_j >= g - r sum_i ||f_i|| - r^2 (1+g),

provided all h_i>=0. All bounds below use integers and outward rounding.
This excludes only this explicit neighborhood, not general 842 arrangements.
"""
from common import *
import math
from fractions import Fraction

def interval_product(a,b):
    v=[a[0]*b[0],a[0]*b[1],a[1]*b[0],a[1]*b[1]]
    return min(v),max(v)

def ceildiv(a,b):return -((-a)//b)

def unit_intervals(row,D):
    ratios=[float(v).as_integer_ratio() for v in row]
    den=max(q for p,q in ratios);p=[a*(den//b) for a,b in ratios]
    squared=sum(v*v for v in p);lo=math.isqrt(squared*D*D);hi=lo+(lo*lo!=squared*D*D)
    result=[]
    for v in p:
        result.append((v*D*D//(hi if v>=0 else lo),ceildiv(v*D*D,lo if v>=0 else hi)))
    return result

def certify(path,out):
    data=np.load(path);x=data['coordinates'];n,d=x.shape;assert n==842
    weights=[int(round(float(v)*2**48)) for v in data['weight']]
    assert min(weights)>=0;L=sum(weights);assert L>0
    D=2**80;xx=[unit_intervals(row,D) for row in x]
    aa=[[[0,0] for _ in range(d)] for _ in range(n)];hh=[[0,0] for _ in range(n)];G=[0,0]
    for i,j,w in zip(data['i'],data['j'],weights):
        if not w:continue
        i=int(i);j=int(j);dot=[0,0]
        for k in range(d):
            p=interval_product(xx[i][k],xx[j][k]);dot[0]+=p[0];dot[1]+=p[1]
            for z in [0,1]:aa[i][k][z]+=w*xx[j][k][z];aa[j][k][z]+=w*xx[i][k][z]
        assert dot[0]>0
        for z in [0,1]:G[z]+=w*dot[z];hh[i][z]+=w*dot[z];hh[j][z]+=w*dot[z]
    assert all(h[0]>=0 for h in hh)
    F=0
    for i in range(n):
        squared=0
        for k in range(d):
            hx=interval_product(hh[i],xx[i][k])
            low=aa[i][k][0]*D*D-hx[1];high=aa[i][k][1]*D*D-hx[0]
            squared+=max(abs(low),abs(high))**2
        root=math.isqrt(squared);F+=root+(root*root!=squared)
    g_lo=Fraction(G[0],L*D*D);g_hi=Fraction(G[1],L*D*D);f_hi=Fraction(F,L*D**3)
    def lower(r):return g_lo-r*f_hi-r*r*(1+g_hi)
    low,high=0,100000
    while low+1<high:
        mid=(low+high)//2
        if lower(Fraction(mid,10**6))>Fraction(1,2):low=mid
        else:high=mid
    radius=Fraction(low,10**6);bound=lower(radius)
    assert radius>0 and bound>Fraction(1,2)
    result=dict(n=n,dimension=d,input_sha256=sha(path),weight_scale=2**48,
                interval_scale_bits=80,weight_sum=L,
                g_lower=str(g_lo),g_upper=str(g_hi),tangent_force_sum_upper=str(f_hi),
                radius=str(radius),max_cosine_lower_bound=str(bound),
                diagnostic_radius=float(radius),diagnostic_lower_bound=float(bound),
                diagnostic_force_sum=float(f_hi),diagnostic_g=float(g_lo),
                exact_integer_check=True,
                scope='Every labelled unit configuration with each point at chord distance <= radius from the normalized input rays has a pair with dot product > 1/2.')
    save_json(out,result);return result

if __name__=='__main__':
    reserve()
    for label,filename in [('improved','local_certificate.json'),('published','published_local_certificate.json')]:
        r=certify(ROOT/'evidence'/f'stress_{label}842.npz',ROOT/filename)
        print(json.dumps(dict(configuration=label,radius=r['radius'],lower_bound=r['diagnostic_lower_bound'],
                              exact_integer_check=r['exact_integer_check'])),flush=True)
