"""Turn a quotient-clique candidate into an independently checked span-coset code."""
import importlib.util,json,pathlib,sys,numpy as np
ROOT=pathlib.Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from discovery.core import RunContext
from discovery.engines.binary_span_cosets import build
z=np.load(sys.argv[1]);K=list(map(int,z['dual']));labels=list(map(int,z['chosen']));basis=list(map(int,z['basis']));n=32;r=len(K);k=n-r;assert 0<k<=len(basis) and 0<r<=31
assert all((a&b).bit_count()%2==0 for a in basis[:k] for b in K)
piv={}
for j in range(n):
 syndrome=sum(((K[i]>>j)&1)<<i for i in range(r));word=1<<j
 while syndrome:
  p=syndrome.bit_length()-1
  if p in piv:syndrome^=piv[p][0];word^=piv[p][1]
  else:piv[p]=(syndrome,word);break
assert len(piv)==r
def representative(label):
 s=label;w=0
 while s:
  p=s.bit_length()-1;assert p in piv;s^=piv[p][0];w^=piv[p][1]
 return w
reps=[representative(x) for x in labels];assert all(sum(((word&K[i]).bit_count()%2)<<i for i in range(r))==label for word,label in zip(reps,labels))
ctx=RunContext(ROOT,{'max_codewords':1<<19,'max_pairs':10000000,'max_rank_work':100000000});code=build([] ,{'length':n,'generators':basis[:k],'representatives':reps},ctx);assert code.facts['min_distance']>=8 and code.facts['cardinality']>131072
spec=importlib.util.spec_from_file_location('verify_kissing',pathlib.Path(__file__).with_name('verify_kissing.py'));v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v);print(json.dumps({'exact_binary_code':code.facts,'verified_kissing_counts':v.verify(code)},indent=2))
if len(sys.argv)>2:pathlib.Path(sys.argv[2]).write_text(json.dumps(code.to_dict()))
