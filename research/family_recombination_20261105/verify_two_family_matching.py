"""Verify concrete matching obstructions for each 32D two-family union."""
import json,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[2];here=pathlib.Path(__file__).resolve().parent
source=json.loads((ROOT/'verifications/improved/dim32-35-signed-layers/source/classic/certificate32.json').read_text());F=source['full_even_supports'];E={x[0] for x in source['exceptional_antipodal_lines']};old=set(F);z=json.loads((here/'permutations.json').read_text())
for filename,filter_exceptions in [('two_family_matching.json',True),('two_family_matching_unfiltered.json',False)]:
 cert=json.loads((here/filename).read_text())
 for entry,perm in zip(cert,z['permutations']):
  Q=[]
  for s in F:
   q=sum(1<<perm[j] for j in range(32) if s>>j&1)
   if q not in old and (not filter_exceptions or all((q&t).bit_count()<=4 for t in E)):Q.append(q)
  match=entry['matched_source_indices'];assert len(Q)==entry['candidate_count']==len(match) and len(set(match))==len(match)
  assert all(type(j) is int and 0<=j<len(F) and (Q[i]&F[j]).bit_count()>4 for i,j in enumerate(match))
  print('permutation',entry['index'],'exception-filtered',filter_exceptions,'candidates',len(Q),'matched; support upper bound',len(F))
