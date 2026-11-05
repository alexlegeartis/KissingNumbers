"""Reconstruct and independently replay a 33D neutral signed-line exchange."""
import json,sys
from pathlib import Path
root=Path(__file__).resolve().parents[2];sys.path.insert(0,str(root))
from discovery.core import RunContext
from discovery.engines.layers import signed_supports
p=root/'verifications/improved/dim32-35-signed-layers/source/classic/certificate33.json';c=json.loads(p.read_text());d=json.loads(Path(__file__).with_name('delta.json').read_text());old=c['exceptional_antipodal_lines'];remove=d['removed_exceptional_line_indices'];add=d['added_antipodal_lines'];assert len(remove)==len(add)==67 and len(set(remove))==67 and all(type(i) is int and 0<=i<len(old) for i in remove)
remaining=[x for i,x in enumerate(old) if i not in set(remove)];original={tuple(x) for x in old};assert len({tuple(x) for x in add})==67 and not original.intersection(map(tuple,add))
certificate={'format':'signed_support_layers_v1','dimension':33,'dense_length':32,'dense_generator':c['dense_generator'],'full_even_supports':c['full_even_supports'],'exceptional_antipodal_lines':remaining+add,'expected_cardinality':370368};a=signed_supports([] ,{'certificate':certificate},RunContext(root,{'max_pairs':8000000}));assert a.facts['cardinality']==370368 and a.verification.status=='exact';print('PASS: 370368 points; 67 old antipodal lines replaced by 67 new lines')
