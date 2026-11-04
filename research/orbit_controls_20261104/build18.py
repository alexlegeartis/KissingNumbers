"""Native weighted-Q(sqrt3) adapter for the known8358-point18D control."""
from pathlib import Path
from itertools import combinations
import json,sys,time,hashlib
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from discovery.core import RunContext
from discovery.engines.geometry import quadratic_code

def main():
 start=time.monotonic();out=Path(__file__).resolve().parent;source=ROOT/'verifications/improved/dim18-bent-hexagon/data'
 octads=json.loads((source/'octads.json').read_text())['octads']
 families=json.loads((source/'families.json').read_text())['families_at_angles_0_60_120']
 slots=json.loads((source/'tierB.json').read_text())['slots_deg_to_words']
 COS={0:(2,0),1:(0,1),2:(1,0),3:(0,0),4:(-1,0),5:(0,-1),6:(-2,0),7:(0,-1),8:(-1,0),9:(0,0),10:(1,0),11:(0,1)}
 A=[];B=[]
 def add(body,radius,angle):
  a=body+[0,0];b=[0]*18
  for position,k in ((16,angle),(17,(3-angle)%12)):
   c,d=COS[k];a[position]=radius*c//2;b[position]=radius*d//2
  A.append(a);B.append(b)
 def signs(support,scale):
  for word in range(1<<len(support)):
   if word.bit_count()%2!=1:continue
   row=[0]*16
   for i,j in enumerate(support):row[j]=scale*(-1 if word>>i&1 else 1)
   yield row
 for i,j in combinations(range(16),2):
  for x in (-12,12):
   for y in (-12,12):
    row=[0]*16;row[i]=x;row[j]=y;add(row,0,0)
 for support in octads:
  for body in signs(support,6):add(body,0,0)
 for angle in range(6):
  for support in families[angle%3]:
   for body in signs(support,6):add(body,6,2*angle)
  add([0]*16,12,2*angle)
 for degree,words in slots.items():
  for word in words:add([4*(-1 if word>>(15-i)&1 else 1) for i in range(16)],4,int(degree)//30)
 params={'radicand':3,'rational_coordinates':A,'radical_coordinates':B,'weights':[1]*16+[2,2]}
 context=RunContext(ROOT,{'max_pairs':40000000})
 artifact=quadratic_code([],params,context)
 if artifact.facts['dimension']!=18 or artifact.facts['cardinality']!=8358:raise ValueError('Known18D control count/rank')
 program={'schema_version':1,'name':'Native weighted quadratic18D control','nodes':[{'id':'sphere','op':'geometry.quadratic_code','params':params}],'outputs':['sphere']}
 (out/'control18_program.json').write_text(json.dumps(program,separators=(',',':'))+'\n')
 (out/'control18_artifact.json').write_text(json.dumps(artifact.to_dict(),separators=(',',':'))+'\n')
 report={'facts':artifact.facts,'verification':artifact.verification.to_dict(),'budget':context.budget,'seconds':time.monotonic()-start,
         'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in source.glob('*.json')},
         'new_sota':False,'scope':'Existing8358 construction, now verified through generic weighted quadratic coefficients; no new construction or held-out discovery claim'}
 (out/'control18_summary.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({'facts':artifact.facts,'pairs':artifact.verification.details['pairs_checked'],'seconds':report['seconds']}),flush=True)
if __name__=='__main__':main()
