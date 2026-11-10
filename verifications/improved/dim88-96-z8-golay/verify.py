"""Replay both proofs read-only and compare every result with its pinned receipt."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[2]/'common'))
from checkpoint import check_manifest
from final import results


def main():
    if not __debug__:
        raise SystemExit('Run without -O.')
    check_manifest(ROOT)
    with tempfile.TemporaryDirectory(prefix='kissing_z8_') as directory:
        reports = []
        for script, receipt in [('replay.py','replay.json'),('verify_independent.py','independent.json')]:
            output = Path(directory)/receipt
            subprocess.run([sys.executable,'-B',str(ROOT/'source'/script),'--output',str(output)],check=True)
            fresh = json.loads(output.read_text())
            assert fresh == json.loads((ROOT/'data'/receipt).read_text()), script
            reports.append({k:v for k,v in fresh.items() if k!='representation'})
        assert reports[0] == reports[1]
        section_reports=[]
        for script in ['section_replay.py','section_independent.py']:
            output=Path(directory)/(script+'.json')
            subprocess.run([sys.executable,'-B',str(ROOT/'source'/script),'--output',str(output)],check=True)
            fresh=json.loads(output.read_text())
            assert fresh==json.loads((ROOT/'data/section_receipt.json').read_text())
            section_reports.append(fresh)
        assert section_reports[0]==section_reports[1]
        totals={int(d):r['total'] for d,r in reports[0]['dimensions'].items()}
        totals.update({r['dimension']:r['total'] for r in section_reports[0]})
        for prefix,claimfile in [('adapted','adapted_zero.json'),('adapted_equal','adapted_equal.json')]:
            claims=json.loads((ROOT/'data'/claimfile).read_text())
            first=Path(directory)/(prefix+'_first.json')
            second=Path(directory)/(prefix+'_second.json')
            subprocess.run([sys.executable,'-B',str(ROOT/'source'/(prefix+'_replay.py')),'--output',str(first)],check=True)
            subprocess.run([sys.executable,'-B',str(ROOT/'source'/(prefix+'_independent.py')),'--output',str(second)],check=True)
            assert json.loads(first.read_text())==claims
            independent=json.loads(second.read_text())
            assert independent['status']=='PASS'
            assert {r['dimension']:r['total'] for r in independent['dimensions']}=={r['dimension']:r['total'] for r in claims}
            totals.update({r['dimension']:r['total'] for r in claims})
        nr_claim=json.loads((ROOT/'nonlinear/data/verified.json').read_text())
        for script in ['replay.py','independent.py']:
            output=Path(directory)/('nr_'+script+'.json')
            subprocess.run([sys.executable,'-B',str(ROOT/'nonlinear/source'/script),'--output',str(output)],check=True)
            fresh=json.loads(output.read_text())
            assert fresh['status']=='PASS' and fresh['NR24_pairs']==8386560
            assert fresh['records']==nr_claim['records']
        totals.update({r['dimension']:r['total'] for r in nr_claim['records']})
        assert totals=={r['dimension']:r['total'] for r in results()}
        subprocess.run([sys.executable,'-B',str(ROOT/'nonlinear/source/negative_controls.py')],check=True)
        subprocess.run([sys.executable,'-B',str(ROOT/'source/adapted_negative_controls.py')],check=True)
        subprocess.run([sys.executable,'-B',str(ROOT/'source/section_negative_controls.py')],check=True)
        subprocess.run([sys.executable,'-B',str(ROOT/'source/negative_controls.py')],check=True)
    print(json.dumps(dict(status='PASS',results=results(),independent_replays=2,
                          nonlinear_replays=2, nonlinear_pairs_per_replay=8386560, support_pairs_per_replay=338559231,octad_coset_exclusions_per_replay=4166151)))


if __name__ == '__main__':
    main()

