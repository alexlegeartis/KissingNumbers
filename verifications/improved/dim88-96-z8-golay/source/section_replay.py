"""Exact new-section proof, using packed sign conditions and projected ring images."""
import argparse,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]

def validate_claim(rec,terms):
 T=rec['rows'];equal=rec['section']=='equal';m=len(T)
 assert rec['section'] in ('zero','equal') and len(set(T))==m and all(type(t)==int and 0<=t<24 for t in T)
 assert rec['dimension']==96-m+int(equal)
 assert rec['terms']==terms and rec['total']==sum(terms.values())
 assert rec['gain']==rec['total']-rec['baseline']>0
 U=np.zeros((m,96),dtype=np.int8)
 for i,t in enumerate(T):U[i,4*t:4*t+4]=1
 A=U[1:]-U[0] if equal else U
 cols=[4*t for t in (T[1:] if equal else T)]
 assert np.array_equal(A[:,cols],np.eye(len(A),dtype=np.int8))
 # Every pair-category product bound is exact; radicals are squared.
 assert 96-2*24==4*(24-2*6)==4*12==2*24==16*3==2*4*6==48
 assert 4*6==24 and 16*(6-4)==32
 assert 4*48<48**2 and 16*48<48**2 and 16*96<48**2
 return A.tolist()

def finish(rec,dense,sparse,sign_matrix):
 T=rec['rows'];m=len(T);equal=rec['section']=='equal'
 q=9 if rec['dimension']==88 else 11 if rec['dimension']<=91 else 13
 assert 6*q<=96-4*m
 transversal=int(np.all(sign_matrix[:,T]==sign_matrix[:,[T[0]]],axis=1).sum())*4**9 if equal else 0
 terms=dict(dense=dense,transversal24=transversal,support24=sparse,support6=32*q**4,pure=192+4*m)
 A=validate_claim(rec,terms)
 return dict(dimension=rec['dimension'],total=sum(terms.values()),terms=terms,constraint_rows=A)

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
 assert __debug__
 raw=json.loads((ROOT/'data/z8_generator.json').read_text());G=np.array(raw['generator'],dtype=np.int64)
 with np.load(ROOT/'data/supports.npz',allow_pickle=False) as z:data=dict(z)
 words=[0]
 for row in data['sign_generator']:
  v=sum(int(a)<<j for j,a in enumerate(row));words +=[x^v for x in words]
 matrix=np.array([[(x>>j)&1 for j in range(24)] for x in words],dtype=np.int8)
 cache={}
 for at in range(25):
  for w in range(min(4,24-at)+1):
   vals=[w-2*((x>>at)&((1<<w)-1)).bit_count() for x in words]
   for value in range(-4,5):cache[at,w,value]=sum(1<<i for i,v in enumerate(vals) if v==value)
 output=[]
 for rec in json.loads((ROOT/'data/sections.json').read_text()):
  T=rec['rows'];m=len(T);equal=rec['section']=='equal';sparse=0
  for lo,hi in zip(data['supports_low'],data['supports_high']):
   support=int(lo)|(int(hi)<<64);at=0;conditions=[]
   for r in range(24):
    w=((support>>(4*r))&15).bit_count()
    if r in T:conditions.append((at,w))
    at+=w
   for value in (range(-4,5) if equal else [0]):
    keep=(1<<len(words))-1
    for at,w in conditions:keep&=cache[at,w,value]
    sparse+=keep.bit_count()
  mod=8 if equal else 4;image={(0,)*m}
  for row in G[:,T]:image={tuple((int(v)+a*int(g))%mod for v,g in zip(old,row)) for old in image for a in range(mod)}
  if equal:
   row_sum=[sum(1-2*((a>>2&1)^((a&1)*x)^((a>>1&1)*y)) for x in (0,1) for y in (0,1)) for a in range(8)]
   flip_sum=[row_sum[a]-2*(1-2*(a>>2&1)) for a in range(8)]
   good=sum(len({row_sum[a] for a in v})==1 for v in image)+sum(len({flip_sum[a] for a in v})==1 for v in image)
  else:good=sum(all(v) for v in image)
  assert (8**12*good)%len(image)==0
  output.append(finish(rec,8**12*good//len(image),sparse,matrix))
 args.output.write_text(json.dumps(output,indent=2)+'\n')
 print(json.dumps({'status':'PASS','sections':{v['dimension']:v['total'] for v in output}}))
if __name__=='__main__':main()
