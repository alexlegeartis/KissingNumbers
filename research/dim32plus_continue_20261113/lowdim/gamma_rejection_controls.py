"""Malformed and norm-preserving geometric controls for the D7 checker."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations_with_replacement
from pathlib import Path
import importlib.util
import json
import subprocess
import sys
import tempfile

HERE=Path(__file__).resolve().parent
GAMMA=HERE.parent/'gamma'
sys.path.insert(0,str(GAMMA))
spec=importlib.util.spec_from_file_location('reviewed_gamma_checker',GAMMA/'verify_extension.py')
checker=importlib.util.module_from_spec(spec);spec.loader.exec_module(checker)


def minimum_gap(c):
    ps=c['profiles'];t=list(map(F,c['heights']))
    dot=lambda a,b:sum(F(x*y,16)for x,y in zip(a,b))
    q=[dot(h,h)for h in ps]
    a=[(8-v*v)/(8-n)for v,n in zip(t,q)]
    assert all(x>0 for x in a) and all(v*v<=2 for v in t)
    return min((4-t[i]*t[j])**2-a[i]*a[j]*(4-dot(ps[i],ps[j]))**2
               for i,j in combinations_with_replacement(range(len(ps)),2))


def main():
    c=json.loads((GAMMA/'d7_compact69.json').read_text())
    baseline69=checker.check(GAMMA/'d7_compact69.json')
    baseline68=checker.check(GAMMA/'d7_compact68.json')
    controls=[]
    with tempfile.TemporaryDirectory(prefix='d7_rejection_controls_')as td:
        directory=Path(td);case=directory/'case.json'
        def reject(name,mutate,original=c,extra=None):
            d=deepcopy(original);mutate(d);case.write_text(json.dumps(d))
            try:checker.check(case)
            except (AssertionError,KeyError,ValueError,RuntimeError,IndexError,TypeError)as exc:
                controls.append(dict(case=name,rejected=True,error_type=type(exc).__name__,detail=str(exc),**(extra or {})))
            else:raise AssertionError('Invalid certificate accepted: '+name)
        reject('duplicate_selected_profile',lambda d:d['profiles'].__setitem__(0,d['profiles'][1]))
        reject('incomplete_profile_domain',lambda d:d['representatives'].pop())
        reject('wrong_orbit_population',lambda d:d['orbit_sizes'].__setitem__(0,d['orbit_sizes'][0]+1))
        reject('odd_moment_power',lambda d:d['tags'][1].__setitem__(6,3))
        reject('degree_beyond_design',lambda d:d['tags'][-1].__setitem__(6,12))
        reject('changed_moment_coefficient',lambda d:d['rows'][0].__setitem__(0,d['rows'][0][0]+1))
        reject('changed_moment_rhs',lambda d:d['rhs'].__setitem__(0,d['rhs'][0]+1))
        reject('changed_objective',lambda d:d['objective'].__setitem__(0,str(F(d['objective'][0])+1)))
        reject('invalid_dual_multiplier',lambda d:d['dual_exact'].__setitem__(0,str(F(d['dual_exact'][0])+1)))
        reject('false_dual_slack',lambda d:d['slack_exact'].__setitem__(0,str(F(d['slack_exact'][0])+1)))
        reject('inflated_exact_mean',lambda d:d.__setitem__('exact_value',str(F(d['exact_value'])+1)))
        reject('inflated_integral_bound',lambda d:d.__setitem__('exact_bound',d['exact_bound']+1))
        reject('inflated_pole_bonus',lambda d:d.__setitem__('with_two_poles',d['with_two_poles']+1))
        mutation=None
        for i in range(len(c['heights'])):
            for height in ['-1','-1/2','0','1/2','1']:
                d=deepcopy(c);d['heights'][i]=height
                gap=minimum_gap(d)
                if gap<0:mutation=(i,height,gap);break
            if mutation:break
        assert mutation is not None
        i,height,gap=mutation
        reject('norm_preserving_cross_layer_failure',lambda d:d['heights'].__setitem__(i,height),
               extra={'changed_height_index':i,'new_height':height,'minimum_squared_gap':str(gap),
                      'all_new_norms_remain_8':True,'pole_height_limit_still_satisfied':True})
        embedding=json.loads((GAMMA/'d7_embedding.json').read_text())
        embedding['normals'][6][0]+=1
        (directory/'d7_embedding.json').write_text(json.dumps(embedding))
        oldhere=checker.HERE;checker.HERE=directory
        reject('changed_actual_integral_normal',lambda d:None)
        checker.HERE=oldhere
    optimized=subprocess.run([sys.executable,'-O',str(GAMMA/'verify_extension.py'),str(GAMMA/'d7_compact69.json')],capture_output=True,text=True)
    assert optimized.returncode!=0 and 'Run without -O' in optimized.stderr
    controls.append({'case':'disabled_assertions','rejected':True,'exit_code':optimized.returncode})
    result={'status':'PASS','valid_baselines':[baseline68['bound'],baseline69['bound']],
            'negative_controls':len(controls),'controls':controls}
    (HERE/'d7_rejection_controls.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
