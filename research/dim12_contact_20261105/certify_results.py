"""Package exact, replayable checks for the main numerical and local results."""
from common import *
from local_certificate import certify
import shutil

def cosine_bound(path,p,q):
    x=np.load(path);rows=[]
    for row in x:
        fractions=[float(v).as_integer_ratio() for v in row];den=max(b for a,b in fractions)
        rows.append([a*(den//b) for a,b in fractions])
    norm=[sum(v*v for v in row) for row in rows];assert min(norm)>0
    failed=0;kiss_failed=0
    for i in range(len(rows)):
        for j in range(i):
            dot=sum(a*b for a,b in zip(rows[i],rows[j]))
            if dot>0:
                product=norm[i]*norm[j]
                failed+=q*q*dot*dot>p*p*product
                kiss_failed+=4*dot*dot>product
    return dict(path=str(path),sha256=sha(path),n=len(x),dimension=x.shape[1],
                threshold_numerator=p,threshold_denominator=q,passed=failed==0,
                violating_pairs=failed,kissing_violations=kiss_failed,pairs=len(x)*(len(x)-1)//2,
                method='exact integers on mathematically normalized binary64 rays')

if __name__=='__main__':
    reserve();out=ROOT/'evidence';out.mkdir(exist_ok=True)
    shutil.copyfile(ROOT/'runs'/'minimax842'/'best.npy',out/'improved842.npy')
    upper=cosine_bound(out/'improved842.npy',500900475,1000000000);assert upper['passed'];assert upper['kissing_violations']>0
    save_json(ROOT/'improved842_certificate.json',upper)
    for name,source,certificate in [
        ('stress_improved842.npz',ROOT/'runs'/'stress842'/'stress_1.npz','local_certificate.json'),
        ('stress_published842.npz',ROOT/'runs'/'stress_published842'/'stress.npz','published_local_certificate.json')]:
        data=dict(np.load(source))
        # The exact certificate for the published realization uses the literal
        # input rays, not an intermediate binary64 renormalization of them.
        if 'published' in name:data['coordinates']=np.load(OLD/'data'/'published842.npy')
        np.savez(out/name,**data);r=certify(out/name,ROOT/certificate)
        print(name,r['radius'],r['diagnostic_lower_bound'],flush=True)
    r=json.loads((ROOT/'runs'/'repair841'/'summary.json').read_text())
    passed=[z for z in r if z.get('exact_certificate',{}).get('passed')]
    assert passed;best=min(passed,key=lambda z:z['final']['max_cos'])
    shutil.copyfile(ROOT/'runs'/'repair841'/f"candidate_{best['index']:03d}.npy",out/'recovered841.npy')
    result=verify(out/'recovered841.npy');assert result['passed'];save_json(ROOT/'recovered841_certificate.json',result)
    print(json.dumps(upper),flush=True)
