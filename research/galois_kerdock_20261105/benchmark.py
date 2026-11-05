"""Reproduce nonlinear trace-code search controls, not a new record."""
from pathlib import Path
import sys,time
root=Path(__file__).resolve().parents[2];sys.path.insert(0,str(root))
from discovery.core import RunContext
from discovery.engines.galois_kerdock import construct,find_translates,search_even_extensions,spherical
from discovery.certification import verify_artifact
ctx=RunContext(root,{'solver_seconds':40,'max_trials':4000,'max_rank_work':60000000});t=time.time();nr=construct([],{'m':3},ctx);assert (nr.facts['binary_length'],nr.facts['cardinality'],nr.facts['minimum_distance'])==(16,256,6)
union=find_translates([nr],{'trials':20,'seed':1514,'seconds':5,'add_d_roots':True},ctx);assert union.facts['cardinality']==2528 and verify_artifact(union,ctx).status=='exact'
k=construct([],{'m':5},ctx);assert (k.facts['binary_length'],k.facts['cardinality'],k.facts['minimum_distance'])==(64,4096,28)
z=search_even_extensions([k],{'count':16,'trials':3000,'seed':0,'seconds':30,'family':'all'},ctx);assert (z.facts['binary_length'],z.facts['cardinality'],z.facts['minimum_distance'])==(96,4096,38)
s=spherical([z],{},ctx);assert s.facts['dimension']==96 and s.facts['rank_exact']==96 and s.facts['cardinality']==4096 and verify_artifact(s,ctx).status=='exact'
print('PASS nonlinear (16,256,6) and (64,4096,28); 2528-point 16D union; searched (96,4096,38), 4096-point full-rank 96D sphere; seconds',round(time.time()-t,1))
