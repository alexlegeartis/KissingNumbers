"""Certify an archived 48-system outside the paired-floor construction."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
import hashlib
from pathlib import Path
import numpy as np
from archive_exact import check
ROOT=Path(__file__).resolve().parent


def main():
    raw=np.load(ROOT/'runs'/'archive_valid'/'valid_004.npy')
    assert check(raw)['passed'];p=2**19-1;q=2**20;vectors=[];norms=[]
    for row in raw:
        ratios=[float(v).as_integer_ratio() for v in row];e=max(d.bit_length()-1 for _,d in ratios)
        z=[n<<(e-(d.bit_length()-1)) for n,d in ratios];vectors.append(z);norms.append(sum(v*v for v in z))
    for i,z in enumerate(vectors):
        assert all(q*q*v*v<=p*p*norms[i] for v in z)
        for j in range(i):
            dot=sum(a*b for a,b in zip(z,vectors[j]))
            assert dot<=0 or q*q*dot*dot<=p*p*norms[i]*norms[j]
    normalized=raw/np.linalg.norm(raw,axis=1,keepdims=True)
    g=normalized@normalized.T;np.fill_diagonal(g,-np.inf)
    np.save(ROOT/'strict48_counterexample.npy',raw)
    report=dict(exact=True,threshold_numerator=p,threshold_denominator=q,
                threshold=float(p/q),margin_below_half=2**-20,
                max_cos_diagnostic=float(g.max()),max_coordinate_diagnostic=float(abs(normalized).max()),
                interpretation='Exact dyadic rays, normalized mathematically',
                no_contact_at_cosine_half=True,
                consequence='Not isometric to any Theorem 1 paired-floor configuration, which has 16 pairs at cosine 1/2.',
                independent_open_perturbation_radius=2**-22,local_dimension=48*5)
    report['source_archive_member']=json.loads((ROOT/'archive_exact.json').read_text())['valid'][4]['name']
    report['coordinate_sha256']=hashlib.sha256((ROOT/'strict48_counterexample.npy').read_bytes()).hexdigest()
    (ROOT/'strict48_certificate.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
