"""Read-only, source-free exact checks for the support witness and triangle bound."""
from fractions import Fraction
import itertools
import json
import math
from source_utils import ROOT

def golay_words():
    generator=sum(1<<i for i in (0,2,4,5,6,10,11))
    result=[]
    for m in range(1<<12):
        word=0
        for j in range(12):
            if (m>>j)&1:word^=generator<<j
        result.append(word|((word.bit_count()&1)<<23))
    assert len(set(result))==4096
    assert {w.bit_count()for w in result}=={0,8,12,16,24}
    return result

def verify_supports():
    supports=[int(s,16)for s in (ROOT/'support31_best.hex').read_text().split()]
    assert len(supports)==len(set(supports))==928
    octads={w for w in golay_words()if w.bit_count()==8}
    assert len(octads)==759
    tail=((1<<7)-1)<<24
    for s in supports:
        assert 0<=s<(1<<31) and s.bit_count()==8
        assert (s&tail).bit_count()>=2 or s in octads
    maximum=max((a&b).bit_count()for a,b in itertools.combinations(supports,2))
    assert maximum<=4
    return {'support_count':928,'maximum_support_intersection':maximum,'dimension':31,'resulting_construction_count':98304+2*31*30+128*928,'improves_current_record':False}

def verify_triangle_bound():
    certificate=json.loads((ROOT/'even24_d10_exact_dual.json').read_text())
    coefficients={int(k):Fraction(v)for k,v in certificate['dual_coefficients'].items()}
    assert coefficients=={1:Fraction(27,32),2:Fraction(1,4),3:Fraction(3,32)}
    def K(k,x):return sum((-1)**j*math.comb(x,j)*math.comb(24-x,k-j)for j in range(max(0,k-24+x),min(k,x)+1))
    def f(x):return Fraction(1)+sum(c*K(k,x)for k,c in coefficients.items())
    for x in range(25):assert f(x)==Fraction((10-x)*(14-x)*(16-x),8)
    assert all(c>=0 for c in coefficients.values())
    assert all(f(d)<=0 for d in range(10,25,2))
    assert f(0)==280
    return {'binary_length':24,'allowed_distances':list(range(10,25,2)),'exact_bound_per_fixed_phase':280,'six_fixed_phase_bound':1680,'scope':'Applies to the six fixed orthogonal-plane triangle phases, not arbitrary triangle orientations.'}

if __name__=='__main__':print(json.dumps({'signed_octad_witness':verify_supports(),'triangle_phase_bound':verify_triangle_bound()},indent=2))
