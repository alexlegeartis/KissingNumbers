"""Exact comparison of the 30 placements of the paper's second 60-point block."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from collections import Counter
import itertools,json
from pathlib import Path
import numpy as np
from audit import clebsch_coefficients,leq,FACTORS

ROOT=Path(__file__).resolve().parent


def place_right(a,b,first,last):
    aa=a.copy();bb=b.copy();order=[first,*[i for i in range(6) if i not in [first,last]],last]
    aa[60:120,6:12]=0;bb[60:120,6:12]=0
    aa[60:120,np.array(order)+6]=a[60:120,6:12]
    bb[60:120,np.array(order)+6]=b[60:120,6:12]
    return aa,bb


def main():
    a,b=clebsch_coefficients();rows=[]
    color={tuple(sorted(edge)):r for r,matching in enumerate(FACTORS) for edge in matching}
    for first,last in itertools.permutations(range(6),2):
        aa,bb=place_right(a,b,first,last)
        pp=aa@aa.T+2*bb@bb.T;qq=aa@bb.T+bb@aa.T
        valid=leq(pp,qq,72);np.fill_diagonal(valid,True);assert valid.all()
        contact=(pp==72)&(qq==0);np.fill_diagonal(contact,False)
        degree=np.sum(contact,axis=1)
        rows.append(dict(first=first,last=last,special_edge_color=color[tuple(sorted([first,last]))],
                         contact_pairs=int(contact.sum()//2),
                         contact_degree_histogram={str(k):v for k,v in sorted(Counter(degree.tolist()).items())},
                         bridge_degree_histogram={str(k):v for k,v in sorted(Counter(degree[120:].tolist()).items())}))
    original=next(r for r in rows if (r['first'],r['last'])==(0,5))
    alternative=next(r for r in rows if (r['first'],r['last'])==(0,1))
    assert original['contact_degree_histogram']!=alternative['contact_degree_histogram']
    report=dict(valid_exact_840_placements=30,original=original,alternative=alternative,
                different_degree_histograms_prove_nonisometry=True,all_placements=rows)
    (ROOT/'placements.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='all_placements'},indent=2))


if __name__=='__main__':main()
