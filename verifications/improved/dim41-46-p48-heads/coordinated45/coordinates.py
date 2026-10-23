"""Coordinates in R^45 from exact coefficients and rational head data.

The exact definition is (a L/sqrt(3), h) / sqrt(r+|h|^2), where
L L^T equals the integer body Gram matrix. Numerical export uses NumPy.
"""
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
def family_coordinates(index,coefficients):
    m=np.load(ROOT/'configuration.npz');g=np.load(ROOT/'geometry.npz');a=np.asarray(coefficients,dtype=np.int64)
    c=int(m['classes'][index]);r=int(m['r3'][index]);q=g['body_gram']
    assert np.all(a%3==g['body_residues'][c]%3)
    assert np.all(np.sum((a@q)*a,axis=1)==3*r)
    P=m['span_projector_numerator']/m['span_projector_denominator'];ev,U=np.linalg.eigh(P);U=U[:,ev>.5]
    h=m['head_numerators'][index]/m['head_denominator'];body=a@np.linalg.cholesky(q)/np.sqrt(3)
    return np.c_[body,np.repeat((h@U)[None],len(a),axis=0)]/np.sqrt(r+h@h)
