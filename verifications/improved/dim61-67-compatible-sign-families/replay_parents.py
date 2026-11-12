"""Replay the unchanged 65D owner geometry and the 66/67D four-support design."""
from pathlib import Path
from itertools import combinations
import importlib.util,json,subprocess,tempfile

HERE=Path(__file__).resolve().parent
REPOSITORY=HERE.parents[2]

def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result);return result

def verify_parents():
    if not __debug__:raise RuntimeError('Run without -O: certificate assertions are required.')
    owners=HERE.parent/'dim64-65-signed-parity'
    construction=module(owners/'construction.py','compatible_sign_owner_geometry')
    with tempfile.TemporaryDirectory(prefix='compatible-sign-owners-') as temporary:
        td=Path(temporary);binary=td/'enumerate'
        subprocess.run(['g++','-O3','-std=c++17',str(owners/'enumerate_binary.cpp'),'-o',str(binary)],check=True)
        construction.HERE=td;construction.ENUMERATOR=binary
        report=construction.construction65()
    assert report['owner_class']['total_points']==68224
    assert report['base64_total']==331737984
    source=HERE.parent/'dim61-68-95-96-code-layers/source/continuation49_72'
    netto=module(source/'netto.py','compatible_sign_netto')
    blocks=netto.blocks68();assert len(blocks)==len(set(blocks))==12529
    triples=set()
    for block in blocks:
        assert len(block)==len(set(block))==4 and 0<=min(block)<=max(block)<68
        for triple in combinations(block,3):assert triple not in triples;triples.add(triple)
    assert len(triples)==50116
    counts={d:sum(max(b)<d for b in blocks) for d in (66,67)}
    assert counts=={66:11088,67:11792}
    return dict(status='PASS',owner65_points=68224,owner65_dense_words=1<<28,
                netto68_blocks=12529,distinct_triples=50116,four_support_counts=counts)

if __name__=='__main__':print(json.dumps(verify_parents(),indent=2))
