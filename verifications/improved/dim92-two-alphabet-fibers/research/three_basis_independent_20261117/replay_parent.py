from pathlib import Path
import importlib.util,json
P=Path(__file__).resolve().parent;S=P/'input'
spec=importlib.util.spec_from_file_location('source_independent',S/'parent_independent.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
_,report=m.parent_proof(json.loads((S/'z8_generator.json').read_text()))
assert report==json.loads((S/'parent_receipt.json').read_text())['z8']
(P/'parent_result.json').write_text(json.dumps(report,indent=2));print('PASS full inherited parent minimum replay; source independent implementation reused, not a newly authored parent solver.')
