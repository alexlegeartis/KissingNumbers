"""Analytical controls for the audit, and a counterexample to tolerance acceptance."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='1'
import ast,itertools,json,sys
from pathlib import Path
import numpy as np
from archive_exact import check as check48
from audit import riesz
from cube49 import value_gradient
from alternative_bridges import bridge_supports
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'dim12_quaternionic_20261104'))
from verify import verify


def main():
    tests=[];rows=[]
    for support in [(0,1,2,3),(0,1,4,5),(2,3,4,5)]:
        for signs in itertools.product([-1,1],repeat=4):
            x=np.zeros(6);x[list(support)]=signs;rows.append(x)
    x=np.array(rows);assert check48(x)['passed'];tests.append('Leech48_exact_positive_control')
    bad=x.copy();bad[-1]=bad[0];assert not check48(bad)['passed'];tests.append('duplicate_rejected')
    bad=x.copy();k=np.flatnonzero(bad[0])[0];bad[0,k]*=1+2**-40
    assert not check48(bad)['passed'];tests.append('tiny_coordinate_violation_rejected')
    assert check48(np.load(ROOT/'strict48_counterexample.npy'))['passed'];tests.append('strict48_verified')
    rng=np.random.default_rng(941);h=1e-6
    x=rng.normal(size=(17,6));d=rng.normal(size=x.shape);d/=np.linalg.norm(d)
    _,g=value_gradient(x,.5)
    fd=(value_gradient(x+h*d,.5)[0]-value_gradient(x-h*d,.5)[0])/(2*h)
    assert abs(fd-np.sum(g*d))<1e-7;tests.append('cube49_gradient_finite_difference')
    x=rng.normal(size=(17,12));d=rng.normal(size=x.shape);d/=np.linalg.norm(d)
    _,g=riesz(x,8.)
    fd=(riesz(x+h*d,8.)[0]-riesz(x-h*d,8.)[0])/(2*h)
    assert abs(fd-np.sum(g*d))<1e-7;tests.append('Riesz_gradient_finite_difference')
    supports,*_=bridge_supports()
    assert all(len(set(a)&set(b))<=2 for a,b in itertools.combinations(supports,2))
    tests.append('bridge_support_intersections_exact')
    # Reproduce the public notebook's spectral reconstruction and predicate on
    # an invalid code whose violation lies inside its explicit tolerance.
    raw=np.zeros((20,12));raw[0,:2]=1;raw[1,[0,2]]=[1,1-1e-7]
    raw[2:11,3:]=np.eye(9);raw[11:,3:]=-np.eye(9)
    raw/=np.linalg.norm(raw,axis=1,keepdims=True);gram=raw@raw.T
    w,v=np.linalg.eigh(gram);idx=np.argsort(w)[::-1][:12]
    recovered=v[:,idx]*np.sqrt(np.maximum(w[idx],0));recovered/=np.linalg.norm(recovered,axis=1,keepdims=True)
    gg=recovered@recovered.T;np.fill_diagonal(gg,-np.inf);maximum=float(gg.max())
    tolerance_accept=maximum<=.5+1e-7
    np.save(ROOT/'checker_false_positive.npy',recovered)
    certificate=verify(ROOT/'checker_false_positive.npy')
    assert tolerance_accept and not certificate['passed'];tests.append('tolerance_false_positive_demonstrated')
    polish=(ROOT/'inputs'/'polish_841.py').read_text();tree=ast.parse(polish);scaled=False
    for node in ast.walk(tree):
        if isinstance(node,ast.Assign) and isinstance(node.value,ast.BinOp) and isinstance(node.value.op,ast.Div):
            if isinstance(node.value.left,ast.Name) and node.value.left.id=='lr' and isinstance(node.value.right,ast.Constant) and node.value.right.value==10:scaled=True
    assert scaled;assert 'assert args.N == 841' in (ROOT/'inputs'/'search_841_riesz.py').read_text()
    code=dict(checker_tolerance=1e-7,false_positive_max_cos=maximum,
              exact_false_positive_certificate=certificate,
              polish_actual_lr_is_printed_lr_divided_by_ten=True,
              released_search_hardcodes_N841=True,
              note='These are code limitations, not a refutation of the valid 841-point coordinate file.')
    (ROOT/'code_findings.json').write_text(json.dumps(code,indent=2)+'\n')
    result=dict(passed=True,checks=tests)
    (ROOT/'checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
