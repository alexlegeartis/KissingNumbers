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
        assert {int(d):r['total'] for d,r in reports[0]['dimensions'].items()} == {
            r['dimension']:r['total'] for r in results()}
        subprocess.run([sys.executable,'-B',str(ROOT/'source/negative_controls.py')],check=True)
    print(json.dumps(dict(status='PASS',results=results(),independent_replays=2,
                          support_pairs_per_replay=338559231,octad_coset_exclusions_per_replay=4166151)))


if __name__ == '__main__':
    main()
