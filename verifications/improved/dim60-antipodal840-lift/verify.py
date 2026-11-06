"""One-command finite replay and exact lift-inequality certificate."""
from pathlib import Path
import subprocess,sys,json,hashlib
from fractions import Fraction as F
if not __debug__:raise SystemExit('Run without -O.')
D=Path(__file__).parent
manifest=json.loads((D/'SOURCE_MANIFEST.json').read_text())
for name,item in manifest['files'].items():
 raw=(D/'source'/name).read_bytes()
 assert len(raw)==item['bytes'] and hashlib.sha256(raw).hexdigest()==item['sha256'],name
for name,args in [('verify_tail.py',[]),('verify_poles.py',[]),('verify_heads.py',['--full'])]:
 subprocess.run([sys.executable,str(D/name),*args],check=True)
t=json.loads((D/'tail_verification.json').read_text());p=json.loads((D/'pole_verification.json').read_text());h=json.loads((D/'head_verification_full.json').read_text())
assert t['points']==840 and t['antipodal_pairs']==h['classes']==420 and p['points']==841
assert t['violations']==p['violations']==0 and h['status']=='PASS'
# Algebraic upper bounds for every pair type, with s=3/4 head energy.
s=F(3,4);assert s*F(1,3)+(1-s)==F(1,2)
assert s*F(1,2)+(1-s)*F(1,2)==F(1,2)
assert s-(1-s)==F(1,2)
assert s*F(1,4)<=F(1,4)  # square of equator-cap maximum <= (1/2)^2
assert 1-s==F(1,4)       # square of pole-cap maximum = (1/2)^2
M=52416000;S=h['selected_lines'];N=M-2*S+4*S+p['points'];assert N==58542821
r={'status':'PASS','dimension':60,'lower_bound':N,'previous_repository_bound':57968343,'gain':N-57968343,'mother_shell_population':M,'selected_mother_lines':S,'direction_count':840,'direction_pairs':420,'independent_poles':841,'exact_tail_pairs':t['pairs_checked'],'exact_pole_pairs':p['pairs_checked'],'exact_base_pairs':h['base_pairs_checked'],'class_isometries':420,'disjointness':'Full canonical 48-coordinate comparison','inherited_inputs':['P48p has minimum squared norm 6 and 52,416,000 minimal vectors; published lattice input, not re-enumerated.'],'scope':'New composition of published constructions; exact finite verification and algebraic pair-type proof, not an expanded all-pairs check of 58 million points. No unrestricted optimality or worldwide-priority claim.'}
(D/'certificate.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
