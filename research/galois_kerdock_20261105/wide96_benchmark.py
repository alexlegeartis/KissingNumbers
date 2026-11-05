"""Search and exactly certify a small nonlinear 96D translate union."""
from pathlib import Path
import sys,time
root=Path(__file__).resolve().parents[2];sys.path.insert(0,str(root))
from discovery.core import RunContext
from discovery.engines.galois_kerdock import construct,search_even_extensions
from discovery.engines.galois_wide_translates import search
from discovery.certification import verify_artifact
ctx=RunContext(root,{'max_trials':4000,'max_pairs':20000000,'max_gram_entries':40000000,'max_rank_work':1000000000,'max_codewords':1000000,'solver_seconds':40});t=time.time();a=construct([],{'m':5},ctx);b=search_even_extensions([a],{'count':16,'trials':3000,'seed':0,'seconds':30,'family':'all'},ctx);assert b.facts['minimum_distance']==38
c=search([b],{'trials':2000,'seed':150923,'seconds':20,'max_layers':100,'add_d_roots':True},ctx);assert c.facts['dimension']==96 and c.facts['translate_layers']==100 and c.facts['cardinality']==427840 and verify_artifact(c,ctx).status=='exact'
print('PASS 100 nonlinear translates plus D96 roots:',c.facts['cardinality'],'points; seconds',round(time.time()-t,1))
