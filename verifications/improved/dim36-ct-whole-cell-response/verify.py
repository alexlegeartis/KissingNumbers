#!/usr/bin/env python3
"""Exact compact verification of the CT pocket refills."""
from pathlib import Path
from collections import defaultdict
import argparse,hashlib,importlib.util,itertools,json,sys
import numpy as np
HERE=Path(__file__).resolve().parent
PARENT_DEFAULT=HERE.parent/"dim36-ct-triangle-refill"
DOUBLE_DEFAULT=HERE.parent/"dim36-ct-double-trine"

def load(name,file):
 spec=importlib.util.spec_from_file_location(name,file)
 mod=importlib.util.module_from_spec(spec);sys.modules[name]=mod;spec.loader.exec_module(mod)
 return mod

def main():
 cli=argparse.ArgumentParser();cli.add_argument("--previous",type=Path,default=PARENT_DEFAULT)
 cli.add_argument("--double",type=Path,default=DOUBLE_DEFAULT)
 args=cli.parse_args();prev=args.previous.resolve();double=args.double.resolve()
 W=json.loads((HERE/"witness.json").read_bytes())
 for k,h in W["SHA256"].items():
  root,file=k.split("/",1)
  path=(prev if root=="triangle" else double)/file
  assert hashlib.sha256(path.read_bytes()).hexdigest()==h,(path,"unexpected ancestor")
 m=load("ct_original",prev/"verify.py")
 g=m.read_geometry();assert g.points==502902
 d=json.loads((double/"witness.json").read_bytes())
 assert d["points"]==502929
 older=set(g.final)
 for r in d["remove_rows"]:
  older.remove(tuple(r))
 for a in d["insert_full_triples"]:
  for p in itertools.product(range(3),repeat=3):
   t=tuple(a)+p
   assert t not in older;older.add(t)
 assert len(older)==99225
 G=g.gram.astype(np.int32)
 lab=np.flatnonzero(g.minima==4)
 lookup={int(a):i for i,a in enumerate(lab)}
 R=g.roots[lab].reshape(-1,12).astype(np.int32)
 T=R@G@R.T
 assert np.all(np.diag(T)==4) and (T-4*np.eye(len(R),dtype=np.int32)).max()<=2
 def convert(rows):
  return np.array([[3*lookup[a]+b for a,b in zip(r[:3],r[3:])] for r in rows],np.int32)
 def decode(mask,n):
  assert type(mask)is int and 0<mask<(1<<n)
  return [i for i in range(n) if mask>>i&1]
 A=set();keys=set()
 for a,b,c,mask in W["add"]:
  assert all(x in lookup for x in [a,b,c])
  assert (a,b,c) not in keys;keys.add((a,b,c))
  for x in decode(mask,27):
   r=(a,b,c,x%3,x//3%3,x//9)
   assert r not in older;A.add(r)
 O=sorted(older);O3=convert(O);N=convert(sorted(A))
 D=sum(T[N[:,i,None],N[:,i]] for i in range(3))
 assert np.all(D.diagonal()==12) and (D-12*np.eye(len(N),dtype=np.int32)).max()<=6
 dead_roots=set();old_root_max=-100
 for lo in range(0,len(N),32):
  X=N[lo:lo+32]
  D=sum(T[X[:,i,None],O3[:,i]] for i in range(3))
  old_root_max=max(old_root_max,int(D.max()))
  dead_roots.update(O[j] for j in np.flatnonzero((D>6).any(0)))
 declared_roots=set(map(tuple,W["remove_root"]))
 assert dead_roots<=declared_roots<=older and len(declared_roots)==len(W["remove_root"])
 dead_roots=declared_roots
 assert (R@G@g.axial.T).max()<=6
 base_drop=set()
 for own,ii,mask in W["remove_base"]:
  assert 0<own<len(g.cells) and sorted(g.minima[g.cells[own]])==[0,4,8]
  s,b,e=m.pairing(g.cells[own],g)
  assert 0<=ii<len(g.cosets[g.cells[own,b]])
  for h in decode(mask,3):
   id=(own,ii,h)
   assert id not in base_drop;base_drop.add(id)
 seen=set();base_pairs=0;base_total=3*len(g.axial)
 def test(blk_a,blk_b,A,B):
  X=R[N[:,blk_a]]@G@A.T;Y=R[N[:,blk_b]]@G@B.T
  return ((X[:,:,None]+Y[:,None,:])>6).any(0)
 for own,cell in enumerate(g.cells[1:],1):
  typ=sorted(g.minima[cell])
  if typ==[4,4,4]:continue
  if typ==[0,6,6]:
   a,b=np.flatnonzero(cell);U,V=g.cosets[cell[a]],g.cosets[cell[b]]
   assert not test(a,b,U,V).any()
   base_pairs+=len(N)*len(U)*len(V);base_total+=len(U)*len(V)
   continue
  assert typ==[0,4,8]
  s,b,e=m.pairing(cell,g)
  V=g.cosets[cell[b]];move=g.switches.get(own,set())
  idx=np.array([tuple(v) in move for v in V],bool)
  for a,l,members in [(s,cell[s],~idx),(e,g.negative[cell[s]],idx)]:
   hi=np.flatnonzero(members)
   if not len(hi):continue
   U=g.roots[l]
   bad=test(a,b,U,V[hi])
   for h,j in zip(*np.nonzero(bad)):seen.add((own,int(hi[j]),int(h)))
   base_pairs+=len(N)*3*len(hi);base_total+=3*len(hi)

 assert seen<=base_drop and base_total==403704
 # Reuse the vacated heavy vector with a new small root in the empty block.
 NEW=[];SK=set()
 for own,ii,small,l,mask in W.get("restore_base",[]):
  cell=g.cells[own]
  assert sorted(g.minima[cell])==[0,4,8] and l in lookup
  s,b,e=m.pairing(cell,g);V=g.cosets[cell[b]]
  assert 0<=ii<len(V) and all((own,ii,ph) in base_drop for ph in range(3))
  assert small!=b and small in range(3)
  for ph in decode(mask,3):
   key=(own,ii,small,l,ph)
   assert key not in SK;SK.add(key)
   Z=np.zeros((3,12),dtype=np.int32)
   Z[b]=V[ii];Z[small]=g.roots[l,ph]
   NEW.append(Z)
 restored_keys={(i,j) for i,j,_,_,_ in SK}
 assert all((i,j) in restored_keys for i,j,_ in base_drop-seen)
 Rr=np.asarray(NEW,dtype=np.int32).reshape(-1,3,12)
 RR=len(Rr);recover_checks=0
 if RR:
  U=sum(Rr[:,b]@G@Rr[:,b].T for b in range(3))
  assert np.all(np.diag(U)==12) and (U-12*np.eye(RR,dtype=np.int32)).max()<=6
  kept=convert(sorted(older-dead_roots))+0
  triple=np.concatenate([kept,N],axis=0)
  mats=[Rr[:,b]@G@R.T for b in range(3)]
  for lo in range(0,RR,16):
   D=sum(mats[b][lo:lo+16,triple[:,b]] for b in range(3))
   assert D.max()<=6
  recover_checks+=RR*len(triple)
  for b in range(3):assert int((Rr[:,b]@G@g.axial.T).max())<=6
  recover_checks+=RR*3*len(g.axial)
  def check_restore(b,c,U,V,keep=None):
   nonlocal recover_checks
   X=Rr[:,b]@G@U.T;Y=Rr[:,c]@G@V.T
   D=X[:,:,None]+Y[:,None,:]
   if keep is not None:
    assert np.all(D[:,keep]<=6)
    recover_checks+=RR*int(keep.sum())
   else:
    assert D.max()<=6
    recover_checks+=RR*len(U)*len(V)
  for own,cell in enumerate(g.cells[1:],1):
   typ=sorted(g.minima[cell])
   if typ==[4,4,4]:continue
   if typ==[0,6,6]:
    b,c=np.flatnonzero(cell);check_restore(b,c,g.cosets[cell[b]],g.cosets[cell[c]])
    continue
   s,b,e=m.pairing(cell,g);V=g.cosets[cell[b]]
   moved=np.array([tuple(v) in g.switches.get(own,set()) for v in V],bool)
   for small,label,members in [(s,cell[s],~moved),(e,g.negative[cell[s]],moved)]:
    j=np.flatnonzero(members)
    if len(j):
     keep=np.array([[(own,int(i),h) not in base_drop for i in j] for h in range(3)],bool)
     check_restore(small,b,g.roots[label],V[j],keep)

 total=base_total+len(older)+len(A)-len(dead_roots)-len(base_drop)+RR
 assert total==W["points"] and total>502929
 out=dict(status="PASS",dimension=36,previous=502929,points=total,
          added=len(A),deleted_base=len(base_drop),deleted_roots=len(dead_roots),restored_base=RR,restored_cross_checks=recover_checks,
          exact_new_old_root_products=len(N)*len(O3),
          exact_new_old_two_block_products=base_pairs,exact_new_new_unordered=len(N)*(len(N)-1)//2)
 print(json.dumps(out,indent=2))

if __name__=="__main__":
 assert __debug__
 main()

