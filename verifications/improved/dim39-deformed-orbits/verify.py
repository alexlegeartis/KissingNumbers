"""Check complete package integrity and replay retained and current exact proofs."""
from pathlib import Path
import json,shutil,sys,tempfile
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[2]/'common'))
from checkpoint import check_manifest,temporary_source,replay
from final import results

if __name__=='__main__':
    check_manifest(ROOT)
    expected=results()[0]
    with temporary_source(ROOT) as source:
        replay(source,'verify.py')
        receipt=json.loads((source/'orbit_verification.json').read_text())
        assert receipt['all_checks_passed'] and receipt['lower_bound']==769180
    with tempfile.TemporaryDirectory(prefix='kissing39_continuation_') as directory:
        source=Path(directory)/'continuation'
        shutil.copytree(ROOT/'continuation',source,ignore=shutil.ignore_patterns('__pycache__'))
        replay(source,'verify.py')
        replay(source,'verify_direct.py')
        first=json.loads((source/'orbit_verification.json').read_text())
        second=json.loads((source/'audit_result.json').read_text())
        assert first['all_checks_passed'] and second['passed']
        assert first['lower_bound']==second['total']==769308
        assert first['certificate_sha256']==second['certificate_sha256']
        replay(source,'negative_controls.py')
        controls=json.loads((source/'negative_control_results.json').read_text())
        assert len(controls)==8 and all(row['rejected'] for row in controls)
    with tempfile.TemporaryDirectory(prefix='kissing39_continuation64_') as directory:
        source=Path(directory)/'continuation64'
        shutil.copytree(ROOT/'continuation64',source,ignore=shutil.ignore_patterns('__pycache__'))
        replay(source,'verify.py')
        replay(source/'independent','audit.py',str(source))
        replay(source/'independent','negative_controls.py')
        first=json.loads((source/'orbit_verification.json').read_text())
        second=json.loads((source/'independent/audit_result.json').read_text())
        assert first['all_checks_passed'] and second['passed']
        assert first['lower_bound']==second['count']==769372
        assert first['certificate_sha256']==second['certificate_sha256']
    with tempfile.TemporaryDirectory(prefix='kissing39_actual_support_') as directory:
        source=Path(directory)/'continuation32'
        shutil.copytree(ROOT/'continuation32',source,ignore=shutil.ignore_patterns('__pycache__'))
        replay(source,'verify.py')
        replay(source/'independent','replay.py')
        first=json.loads((source/'orbit_verification.json').read_text())
        second=json.loads((source/'independent/complete_audit.json').read_text())
        assert first['all_checks_passed'] and second['status']=='PASS'
        assert first['lower_bound']==second['points']==769436
        assert first['certificate_sha256']==second['certificate_sha256']
        assert second['negative_controls_passed']==2
        replay(source,'negative_controls.py')
        controls=json.loads((source/'negative_control_results.json').read_text())
        assert len(controls)==4 and all(row['rejected'] for row in controls)
    with tempfile.TemporaryDirectory(prefix='kissing39_replacement28_') as directory:
        source=Path(directory)/'replacement28'
        shutil.copytree(ROOT/'replacement28',source,ignore=shutil.ignore_patterns('__pycache__','_build'))
        replay(source,'verify.py')
        replay(source/'independent','verify_independent.py')
        replay(source/'independent','check_witnesses.py')
        replay(source/'independent','negative_controls.py')
        first=json.loads((source/'orbit_verification.json').read_text())
        second=json.loads((source/'independent/independent_report.json').read_text())
        assert first['all_checks_passed'] and second['all_checks_passed']
        assert first['lower_bound']==second['lower_bound']==expected['total']
        assert first['certificate_sha256']==second['certificate_sha256']
        assert json.loads((source/'independent/witness_checks.json').read_text())['passed']
        controls=json.loads((source/'independent/negative_controls.json').read_text())
        assert controls['all_rejected'] and len(controls['controls'])==4
    print(json.dumps(dict(status='PASS',results=[expected],proof_routes=2,negative_controls=20)))
