"""Recompute the exact signed-point symmetric difference, without search code."""
import hashlib,json
from pathlib import Path
D=Path(__file__).resolve().parent
old=json.loads((D/'control32_352800.json').read_bytes())
new=json.loads((D/'certificate32_352832.json').read_bytes())
delta=json.loads((D/'delta.json').read_text())
assert hashlib.sha256((D/'control32_352800.json').read_bytes()).hexdigest()=='3a96dc79d21d07557a221f2cdd6664d2d5270e93d7c98711d7b962a2ba2b41d2'
assert old['dimension']==new['dimension']==32
assert old['dense_generator']==new['dense_generator']
assert old['full_even_supports']==new['full_even_supports']
def points(c):
 return {(s,n) for s,n in c['exceptional_antipodal_lines']}|{(s,s^n) for s,n in c['exceptional_antipodal_lines']}
a,b=points(old),points(new)
assert sorted(a-b)==list(map(tuple,delta['removed']))
assert sorted(b-a)==list(map(tuple,delta['added']))
assert len(a-b)==delta['removed_points']==380
assert len(b-a)==delta['added_points']==412
assert new['points']-old['points']==delta['gain']==32
assert delta['old_total']==old['points']==352800
assert delta['new_total']==new['points']==352832
assert sorted({s for s,n in a^b})==delta['changed_supports']
assert {s for s,n in a}=={s for s,n in b}
print('PASS: exact 380-point removal and 412-point addition; all other layers unchanged')
