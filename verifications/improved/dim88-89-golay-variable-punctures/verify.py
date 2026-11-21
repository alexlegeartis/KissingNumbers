import itertools, json, subprocess, sys
from pathlib import Path
P=Path(__file__).resolve().parent
Q=P.parent/'dim89-two-color-latitudes'
dense=json.loads(subprocess.run([sys.executable,str(Q/'verify.py')],check=True,text=True,capture_output=True).stdout)
# The exact dense replay also checks the unchanged sparse input's
# source generator, dodecacode/Julín distances, and all old pattern pairs.
U=Q/'abbreviated_parent'
basis=[sum(((c&1)<<(i+j)) for j,c in enumerate([7,1,4,4,2,7,3,3,4,7,2,1]))+(1<<23)
       for i in range(12)]
C=[0]
for u in basis:C += [v^u for v in C]
assert len(set(C))==4096
params=json.loads((U/'patterns.json').read_text())[0]
assert params['mask']==3 and params['pivot']==3
old11=[(C[m%4096]&~8,C[m//4096]&~8) for m in params['ids']]
old10=[tuple(x)for x in json.loads((U/'weight10.json').read_text())]
new=[tuple(x)for x in json.loads((P/'patterns_extra.json').read_text())]
assert len(new)==353 and len(set(new))==353
mask=(1<<24)-1
W=lambda v: (v[0]|v[1]).bit_count()
D=lambda v,w: ((v[0]^w[0])|(v[1]^w[1])).bit_count()
assert all(0<=a<=mask and 0<=b<=mask and (a|b)&3==0 and W((a,b))==11 for a,b in new)
for a,b in itertools.combinations(new+old11,2):assert D(a,b)>=11
for a,b in itertools.product(new,old10):assert D(a,b)>=11
# Given weights r and s, common occupied plus equal-tag count is r+s-D.
# The row dot is at most 4 or 8 by distinct or equal matchings.  At weight
# 11+11 this proves max raw dot44; at weight 11+10 it proves max40.
assert 4*40**2<=88*80
ext=len(new)*1024*144
out=[(88,12*8**10 + dense['retained_sparse'],12*8**10 + dense['retained_sparse']+ext),
 (89,dense['previous'],dense['lower_bound']+ext)]
actual=[{'dimension':d,'previous':pre,'lower_bound':total,'added_over_previous':total-pre}
for d,pre,total in out]
assert actual==json.loads((P/'claims.json').read_text())
print(json.dumps({'additional_11_patterns':len(new),'additional_11_rays':ext,
 'full_weight_11_patterns':len(new)+len(old11),'old_weight_10_patterns':len(old10),'claims':actual},indent=2))
