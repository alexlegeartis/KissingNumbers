"""Replay two independent exact proofs for every current configuration."""
from pathlib import Path
import json,os,shutil,subprocess,sys,tempfile
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[2]/'common'))
from checkpoint import check_manifest,temporary_source
from final import results

def replay(source, script):
    # Respect the one-thread proof budget without changing other packages.
    path=Path(source)/script
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
    subprocess.run([sys.executable,'-B',str(path)],cwd=path.parent,env=env,check=True)


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
            assert {row['dimension']:row['points'] for row in rows}=={33:370376,34:390990,35:415334}
        assert [(row['dimension'],row['certificate_sha256']) for row in symbolic]==[(row['dimension'],row['certificate_sha256']) for row in coordinate]
    with tempfile.TemporaryDirectory(prefix='rotated34_replay_') as directory:
        current=Path(directory)/'rotated34'
        shutil.copytree(ROOT/'rotated34',current,ignore=shutil.ignore_patterns('__pycache__'))
        replay(current,'verify.py')
        replay(current,'verify_coordinates.py')
        replay(current,'negative_controls.py')
        symbolic=json.loads((current/'verification.json').read_text())
        coordinate=json.loads((current/'coordinate_verification.json').read_text())
        controls=json.loads((current/'negative_control_results.json').read_text())
        assert symbolic['all_checks_pass'] and coordinate['all_checks_pass']
        assert symbolic['points']==coordinate['points']==expected[34]==390992
        assert symbolic['certificate_sha256']==coordinate['input_sha256']['certificate.json']
        assert coordinate['unique_rows']==390992 and coordinate['exact_rank']==34
        assert len(controls)==10 and all(row['rejected'] for row in controls)
    with tempfile.TemporaryDirectory(prefix='exceptional32_replay_') as directory:
        current=Path(directory)/'exceptional32'
        shutil.copytree(ROOT/'exceptional32',current,ignore=shutil.ignore_patterns('__pycache__'))
        replay(current,'replay.py')
        report=json.loads((current/'report.json').read_text())
        assert report['status']=='PASS' and report['certificate_audit']['total']==352832  # retained history
    assert expected[32]==352992
    replay(ROOT/'nonantipodal32','replay.py')
    assert expected[33]==370408
    replay(ROOT/'nonantipodal33','replay.py')
    print(json.dumps(dict(status='PASS',results=[dict(dimension=d,total=n) for d,n in sorted(expected.items())])))
