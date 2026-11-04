"""Expand the oracle for candidate triples; no sampled distance is a certificate."""
from pathlib import Path
from collections import defaultdict,Counter
import json,time,hashlib
import numpy as np
from search import info_set,H,SEED,RAW
from half_cosets import restricted_screen
OUT=Path(__file__).resolve().parent
start=time.monotonic();rng=np.random.default_rng(SEED+3);G=np.array(json.loads(RAW.read_text())['generator'],dtype=np.int64)
data=json.loads((OUT/'half_cosets.json').read_text());groups=defaultdict(list)
for proposal_index,t in enumerate(data['unverified_triples']):groups[tuple(t['functional'])].append({**t,'proposal_index':proposal_index})
infos=[info_set(G,rng.permutation(24))for _ in range(512)];rows=[];survivors=[];proofs=[]
for ell,triples in groups.items():
 offsets=np.array([v for t in triples for v in [t['offsets'][0],t['offsets'][1],(np.array(t['offsets'][0])-np.array(t['offsets'][1]))%8]],dtype=np.int64)
 d,c=restricted_screen(offsets,G,np.array(ell),infos,rng)
 for index,t in enumerate(triples):
  ds=d[3*index:3*index+3];best=int(ds.min())
  if best>=24:survivors.append({**t,'expanded_proxies':ds.tolist(),'status':'unverified after512information sets'})
  else:
   j=3*index+int(ds.argmin());proofs.append(dict(proposal_index=t['proposal_index'],functional=list(ell),offset=offsets[j].tolist(),parent_message=c[j].tolist(),distance=int(d[j])))
 rows.append(dict(functional=list(ell),triples=len(triples),survivors=int(sum(min(d[3*i:3*i+3])>=24 for i in range(len(triples))))))
 print(rows[-1],flush=True)
receipt=dict(status='all survivors remain uncertified',input_triples=len(data['unverified_triples']),rejected=len(proofs),survivors=len(survivors),information_sets=512,rows=rows,surviving_triples=survivors,rejection_witnesses=proofs,elapsed_seconds=time.monotonic()-start,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(OUT/'triple_results.json').write_text(json.dumps(receipt,indent=2)+'\n');print('SUMMARY',len(proofs),len(survivors),receipt['elapsed_seconds'],flush=True)
