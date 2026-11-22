#!/usr/bin/env python3
"""Replay exact P48 owner migration, post-refill, and the combined tail palettes."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import argparse,json,sys,subprocess
from pathlib import Path
from verify_refill import HERE, run
from base import SHA, np
PIN={
 'tail55/verify.py':'5cdeda9b2a8ac70594219a62130163eea15533972e1655f6358acaa8a035c5c9',
 'tail55/witness.json':'cc4151354d9b6eaf6551275a2d9424eb5b8eba11e0a4c3bf5db0bec31e9aa053',
 'tail60/verify.py':'1a5bf0db29ca00dd48850e80b3c007420df2bb3da42478ea979b2192baf48a4d',
 'tail60/witness.json':'d203b2ec80be235d57a7ff37fcfb63c3c2361dac5df60cb3e528499bb76d70e1',
}
def checked(p,sha):
 b=p.read_bytes()
 assert SHA(b)==sha,p
 return b
def replay(p,s,args=()):
 a=[sys.executable,str(p),'--siblings',str(s),*args]
 return [json.loads(r) for r in subprocess.check_output(a,text=True).splitlines()]
def verify(a):
 if not __debug__:raise SystemExit('Run with assertions enabled.')
 for k,v in PIN.items():checked(HERE/k,v)
 dims=set(a.dimensions or range(50,61));assert dims<=set(range(50,61));s=a.siblings.resolve()
 if 54 in dims:
  path=a.tail/'d54/verify.py'
  checked(path,"fa50996b6c2106511b31ff5786429564e26a3f15a05a4b13afdc854f53429188")
  assert replay(path,s)[0]['lower_bound']==52976779
 ov={}; ds=sorted(dims&{55,57,59})
 if ds:
  wat=json.loads((HERE/'tail55/witness.json').read_text())
  ans=replay(HERE/'tail55/verify.py',s,['--dims',*map(str,ds)])
  dt={55:(53393714,42),57:(54761045,88),59:(56952638,56)}
  assert [r['dimension'] for r in ans]==ds
  for r in ans:
   d=r['dimension'];total,extra=dt[d]
   assert r['lower_bound']==total and r['add']==extra
   ov[d]=np.array(wat[str(d)]['heads'],'i8').reshape(-1,48),extra
 if 60 in dims:
  r,=replay(HERE/'tail60/verify.py',s)
  assert r['lower_bound']==58636706 and r['heavy_subcosets']==20
  w=json.loads((HERE/'tail60/witness.json').read_text())
  H=np.concatenate([np.array(w['frame'],'i8'),np.array(w['sparse'],'i8').reshape(-1,48),
    *[np.array(v,'i8') for u in w['upgrade'] for v in u['tetrahedra'] if v is not None]])
  assert H.shape==(112,48)
  ov[60]=H,r['lower_bound']-58631571-4799
 run(a,ov,json.loads((HERE/'results.json').read_text()))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--siblings',type=Path,default=HERE.parent)
 p.add_argument('--tail',type=Path,default=HERE.parent/'dim54-60-p48-tail-antipodes')
 p.add_argument('dimensions',type=int,nargs='*')
 verify(p.parse_args())
