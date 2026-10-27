"""Derive achieved counts from explicit coordinate certificates (not a proof replay)."""
from pathlib import Path
import numpy as np
import hashlib
ROOT=Path(__file__).resolve().parent
DATA=ROOT/'source/27d_local_replacements_20261026'
CONFIG={26:(None,201108),27:(99,209532),28:(154,216680),29:(230,226520),30:(325,238808),31:(385,246638)}
def results():
    rows=[]
    for n,(g,expected) in sorted(CONFIG.items()):
        base='certificate26_200918.npz' if n==26 else 'base%d_groups%d.npz'%(n,g)
        ext='extension26_201108.npz' if n==26 else 'extension%d_groups%d.npz'%(n,g)
        with np.load(DATA/base,allow_pickle=False) as b, np.load(DATA/ext,allow_pickle=False) as e:
            assert str(e['base_certificate_sha256'])==hashlib.sha256((DATA/base).read_bytes()).hexdigest()
            removed=b['removed_leech']; inserted=b['new_points']; free=e['free_heads']
            assert removed.shape[1]==24 and inserted.shape[1]==free.shape[1]==n
            assert all(x.dtype.kind in 'iu' for x in (removed,inserted,free))
            pure=len(e['pure_tail_points']) if n==26 else 2*(n-24)*(n-25)
            if n>=30:
                with np.load(DATA/('pure_tail%d_exceptional.npz'%n),allow_pickle=False) as t:
                    assert t['rational'].shape==t['radical'].shape==((72,6) if n==30 else (126,7))
                    assert t['rational'].dtype.kind in 'iu' and t['radical'].dtype.kind in 'iu'
                    assert int(t['radicand'])==(3 if n==30 else 2)
                    pure=len(t['rational'])
            total=196560-len(removed)+len(inserted)+len(free)+pure
        assert total==expected
        rows.append(dict(dimension=n,total=total))
    return rows
if __name__=='__main__':
    import json
    print(json.dumps(results()))
