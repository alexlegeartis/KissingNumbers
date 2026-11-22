#!/usr/bin/env python3
"""Exact weighted binary distance in the 62D omitted fibre."""
import hashlib, importlib.util, json
from fractions import Fraction as Q
from pathlib import Path
HERE=Path(__file__).resolve().parent
def load(path, sha):
    assert hashlib.sha256(path.read_bytes()).hexdigest()==sha
    spec=importlib.util.spec_from_file_location('pinned_parent',path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    return mod
def verify():
    assert __debug__
    old=load(HERE.parent/'dim61-62-omitted-fibre-spikes/verify.py',
             '29749fd2a4019dacfd67c0ce89e744e70e627c4a270aee6be08c487505a4405f')
    rec=next(x for x in old.verify()['dimensions'] if x['dimension']==62)
    assert rec['lower_bound']==124625036 and rec['points_per_axis']==32768
    src=old.pin('dim61-67-compatible-sign-families/verify_family.py','222c034a630bfab0db9df1249797d0e633436bfc')
    m=importlib.util.spec_from_file_location('old_dense',src)
    a=importlib.util.module_from_spec(m);m.loader.exec_module(a)
    P,K=old.split(a.generator())
    assert len(K)==24
    rows={int(i):tuple(J) for i,J in json.loads((HERE/'rows.json').read_text()).items()}
    assert set(rows)<=set(old.body) and len(rows)==60
    a=Q(9,2);b2=Q(64)/(59+a*a)
    assert (59+a*a-2*20)*b2<=32
    assert (59+a*a-2*(16+a*a-1))*b2<=32
    translates={int(i):t for i,t in json.loads((HERE/'translations.json').read_text()).items()}
    assert set(translates)==set(rows)==set(old.body)
    used=bytearray(1<<24);counts={}
    for p,J in sorted(rows.items()):
        assert len(J)==16 and all(0<x<1<<24 for x in J)
        t=translates[p];assert 0<=t<1<<24
        C=old.span(J)
        assert len(C)==len(set(C))==65536
        Z=old.span([old.choose_rows(K,x) for x in J])
        assert len(Z)==len(set(Z))==65536
        assert all(old.prefix(w)==0 for w in Z)
        assert all(w.bit_count()>=16 and (w.bit_count()>=20 or bool(w>>p&1))
                   for w in Z[1:])
        k=0
        for x in C:
            i=x^t
            if not used[i]:
                used[i]=1
                k+=1
        counts[p]=k
    assert len(counts)==60
    assert (60-2*16+2*(a-1))*b2<=32
    n=sum(counts.values())
    assert n==sum(used)
    return dict(dimension=62, previous=rec['lower_bound'], gain=n-rec['inserted'],
                lower_bound=rec['construction_base']+n-rec['removed_poles'],
                inserted=n,removed_poles=rec['removed_poles'],
                retained_per_physical_axis=counts,duplicates=60*65536-n)
if __name__=='__main__':print(json.dumps(verify(),indent=2))
