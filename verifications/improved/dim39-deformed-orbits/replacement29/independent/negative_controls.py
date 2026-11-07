#!/usr/bin/env python3
"""Exercise geometric rejection in the independent C++ engine, after a full replay.
These deliberately bypass the outer frozen-file hash gate to test the geometry.
"""
import hashlib,json,math,pathlib,subprocess,tempfile,time,os
P=pathlib.Path(__file__).resolve().parent
if hasattr(os, 'sched_getaffinity'):
    cpus=sorted(os.sched_getaffinity(0))
    if len(cpus)>1:os.sched_setaffinity(0,cpus[:-1])
lines=(P/'_build/audit_input.txt').read_text().splitlines()
exe=P/'_build/audit_cpp'
if not exe.exists():raise SystemExit('Run verify_independent.py first.')
original=lines[2].split(); cases=[]
zero=lines.copy();zero[2]=' '.join(['0']*39);cases.append(('zero_norm',zero,'coordinate safety bound'))
badnorm=lines.copy();badnorm[2]=' '.join(['1']*39);cases.append(('nonzero_wrong_norm',badnorm,'norm enclosure'))
oversized=lines.copy();row=original.copy();row[0]='100000005';oversized[2]=' '.join(row);cases.append(('coordinate_safety_guard',oversized,'coordinate safety bound'))
rootstart=2+10269+3424
rootdup=lines.copy();rootdup[rootstart+1]=rootdup[rootstart];cases.append(('duplicate_sparse_root',rootdup,'root-root'))
heavystart=rootstart+2908
heavydup=lines.copy();heavydup[heavystart+1]=heavydup[heavystart];cases.append(('duplicate_sparse_heavy',heavydup,'heavy-heavy'))
duplicate=lines.copy();duplicate[3]=duplicate[2];cases.append(('duplicate_dense_direction',duplicate,'dense-dense'))
n=10269; support=[int(x) for x in lines[2+n].split()]
body=math.isqrt((10**16-8*24000000**2)//31)
v=[body]*39
for k in support:v[k]=24000000
norm=sum(x*x for x in v)
assert 99999996**2<=norm<=100000004**2
collision=lines.copy();collision[2]=' '.join(map(str,v));cases.append(('middle_collision_with_valid_norm',collision,'dense-middle'))
count=lines.copy();header=count[0].split();header[0]='10270';count[0]=' '.join(header);cases.append(('wrong_representative_count',count,'unexpected sizes'))
report=[]
with tempfile.TemporaryDirectory(prefix='independent39_controls_') as temp:
    for name,data,reason in cases:
        path=pathlib.Path(temp)/(name+'.txt');path.write_text('\n'.join(data)+'\n')
        start=time.monotonic();run=subprocess.run([str(exe),str(path)],capture_output=True,text=True)
        rejected=run.returncode!=0 and 'FAIL: '+reason in run.stderr
        if not rejected:raise RuntimeError((name,run.returncode,run.stdout,run.stderr))
        report.append({'name':name,'rejected':rejected,'returncode':run.returncode,'reason':reason,'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start})
result={'all_rejected':True,'controls':report,'scope':'Controls bypass frozen-input SHA256 gates to test actual C++ integer geometry. They are not the source of acceptance for the original certificate.'}
(P/'negative_controls.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
