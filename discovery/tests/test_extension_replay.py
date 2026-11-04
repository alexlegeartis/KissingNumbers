"""The bounded negative proof must reject empty or disabled evidence."""
from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from discovery.benchmarks.replay_extensions import replay

SOURCE=Path(__file__).parents[1]/'benchmarks/all_linear_extensions_20261104.json'

class ExtensionReplayTests(unittest.TestCase):
    def test_rejects_empty_missing_duplicate_and_rank_deficient_proofs(self):
        original=json.loads(SOURCE.read_text())
        cases=[{'runs':[]}]
        for label in ('missing','duplicate','constraints','rank','zero_message'):
            d=deepcopy(original)
            if label=='missing':d['runs']=d['runs'][:1]
            elif label=='duplicate':d['runs'][1]['factor_index']=0
            elif label=='constraints':d['runs'][0]['obstruction_constraints']=[]
            elif label=='rank':d['runs'][0]['generator_without_extension'][1]=d['runs'][0]['generator_without_extension'][0][:]
            else:d['runs'][0]['obstruction_constraints'][0]['message_coefficients']=[0]*9
            cases.append(d)
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'bad.json'
            for data in cases:
                path.write_text(json.dumps(data))
                with self.subTest(data=data['runs'][:1]),self.assertRaises(ValueError):replay(path)
    def test_optimized_python_cannot_skip_empty_obstruction(self):
        data=json.loads(SOURCE.read_text());data['runs'][0]['obstruction_constraints']=[]
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'bad.json';path.write_text(json.dumps(data))
            result=subprocess.run([sys.executable,'-O','-c',
               'from discovery.benchmarks.replay_extensions import replay; import sys; replay(sys.argv[1])',str(path)],
               cwd=Path(__file__).parents[2],text=True,capture_output=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('ValueError',result.stderr)

if __name__=='__main__':unittest.main()
