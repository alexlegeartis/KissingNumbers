"""Replay two independent exact proofs for every current configuration."""
from pathlib import Path
import json,shutil,sys
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
                target={32:352800,33:370368,34:390982,35:415310} if family=='classic' else {35:415314}
                assert observed==target
        target=source.parent/'root_trades'
        shutil.copytree(ROOT/'root_trades',target,ignore=shutil.ignore_patterns('__pycache__'))
        for script in ('verify.py','verify_coordinates.py'):replay(target,script)
        symbolic=json.loads((target/'verification.json').read_text())
        coordinate=json.loads((target/'coordinate_verification.json').read_text())
        for rows in (symbolic,coordinate):
            assert all(row['all_checks_pass'] for row in rows)
            assert {row['dimension']:row['points'] for row in rows}=={d:expected[d] for d in (33,34,35)}
        assert [(row['dimension'],row['certificate_sha256']) for row in symbolic]==[(row['dimension'],row['certificate_sha256']) for row in coordinate]
    print(json.dumps(dict(status='PASS',results=[dict(dimension=d,total=n) for d,n in sorted(expected.items())])))
