"""Exact second sparse layer: seven-site q=2√2 and six ordinary signs."""
from itertools import combinations,product
from pathlib import Path
import json
import numpy as np
from supports37 import supports
P=Path(__file__).resolve().parent

def replay_six(X):
    data=json.loads((P/'six_recipe.json').read_text())
    R=data['ordinary_supports']
    assert len(R)==len(set(R)) and all(0<x<1<<36 and x.bit_count()==6 for x in R)
    assert max((u&v).bit_count() for u,v in combinations(R,2))<=2
    S=supports()
    assert max((x&y).bit_count() for x in R for y in S)<=4
    # Coefficients of the special coordinate are ±2, in a sqrt(2) basis.
    # Product with a prior integer star is a+b*sqrt(2), bounded by both
    # absolute terms; use (16-a)^2 >= 2b^2 with 16-a nonnegative.

    A=np.abs(X[:,36])
    for w in R:
        cols=[i for i in range(36) if w>>i&1]
        lhs=16-2*np.abs(X[:,cols]).sum(1)
        assert (lhs>=0).all() and (lhs*lhs>=8*A*A).all()
    # Worst dense products: 24+2m sqrt(2) against norm-squared 144+m^2.
    # Use the strict rational upper sqrt(2)<99/70.
    assert 99**2>2*70**2
    for m in (1,4):
        assert 4*(24*70+2*m*99)**2<=32*(144+m*m)*70**2
    # Every root on the special coordinate was removed in the first stage.
    centers=X[np.abs(X[:,36])==4]
    assert np.any(centers[:,:36]!=0,axis=0).all()
    # All new pair types are norm-32. For equal special signs, same support
    # even parity has d>=2, and distinct support overlap <=2.
    # Different special signs have raw dot <= -8+24.
    assert 8+4*(6-2*2)==8+4*2==-8+4*6==16
    n=64*len(R)
    return dict(supports_six=len(R),seven_site_points=n,total=534008+n)
