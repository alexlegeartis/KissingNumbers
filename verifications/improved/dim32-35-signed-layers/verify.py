"""Replay two independent exact proofs for every current configuration."""
from pathlib import Path
import json,sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[2]/'common'))
from checkpoint import check_manifest,temporary_source,replay
from final import results

if __name__=='__main__':
    check_manifest(ROOT)
    expected={r['dimension']:r['total'] for r in results()}
    with temporary_source(ROOT) as source:
        for family in ('classic','multi35'):
            replay(source/family,'verify.py')
            replay(source/family,'verify_coordinates.py')
            symbolic=json.loads((source/family/'verification.json').read_text())
            coordinate=json.loads((source/family/'coordinate_verification.json').read_text())
            for rows in (symbolic,coordinate):
                observed={}
                for row in rows:
                    assert row.get('all_checks_pass',row.get('verified',False))
                    d=row['dimension'];n=row.get('verified_points',row.get('points'))
                    assert d not in observed;observed[d]=n
                target=dict(expected) if family=='classic' else {35:expected[35]}
                if family=='classic':target[35]=415310  # retained earlier witness
                assert observed==target
    print(json.dumps(dict(status='PASS',results=[dict(dimension=d,total=n) for d,n in sorted(expected.items())])))
