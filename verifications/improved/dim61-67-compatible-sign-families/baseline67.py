"""Fresh replay of the unchanged full-class family used by the 67D extension."""
from pathlib import Path
import hashlib,json,tempfile
import numpy as np
from verify_family import SOURCE,check
from expanded_geometry import check_expanded

def verify_baseline67():
    if not __debug__:raise RuntimeError('Run without -O: certificate assertions are required.')
    source=SOURCE/'shared_parity/expanded67_68.npz'
    certificate=dict(format='explicit-hamming-sign-classes-v1',removed_coordinates=[],odd_total_parity=False,
        original_supports_retained=30828,source_archive_sha256=hashlib.sha256((SOURCE/'dim61_96/code64_supports.npz').read_bytes()).hexdigest(),
        discovery_seed_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
    with np.load(source,allow_pickle=False) as archive:
        for kind in ('old','new'):
            certificate[kind]=[dict(support=int(s),checks=list(map(int,h)),shift=int(t)) for s,h,t in zip(archive[kind+'_supports'],archive[kind+'_duals'],archive[kind+'_shifts'])]
    assert len(certificate['new'])==len(certificate['old'])==1792
    with tempfile.TemporaryDirectory(prefix='full-baseline67-') as directory:
        path=Path(directory)/'certificate.json';path.write_text(json.dumps(certificate))
        family=check(path,cpp=True);expanded=check_expanded(path,67,cpp=True)
    assert expanded['total']==355090822
    return dict(status='PASS',family=family,expanded=expanded)

if __name__=='__main__':print(json.dumps(verify_baseline67(),indent=2))
