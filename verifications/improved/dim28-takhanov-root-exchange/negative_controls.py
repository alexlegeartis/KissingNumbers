"""Semantic rejection controls for the standalone geometric proof."""
import json
import numpy as np
from verify import integer_geometry,exchange_check

def reject(name,call):
    try:call()
    except ValueError:return name
    raise RuntimeError('malformed example was accepted: '+name)

def main():
    good=np.array([[8,8,0,0,0],[8,0,8,0,0],[-8,-8,0,0,0]],np.int16)
    integer_geometry(good)
    cases=[]
    cases.append(reject('duplicate',lambda:integer_geometry(np.r_[good,good[:1]])))
    bad=good.copy();bad[0,0]=7
    cases.append(reject('wrong_norm',lambda:integer_geometry(bad)))
    # Both bad points have norm128, while their inner product is96.
    clash=np.array([[8,8,0,0,0],[8,4,4,4,4]],np.int16)
    cases.append(reject('norm_correct_conflict',lambda:integer_geometry(clash)))
    cases.append(reject('noninteger_dtype',lambda:integer_geometry(good.astype(float))))
    empty=np.empty((0,5),np.int16)
    cases.append(reject('false_exchange_identity',lambda:exchange_check(good,good[:2],empty,empty)))
    cases.append(reject('missing_source_removal',lambda:exchange_check(good,good,clash[1:],clash[1:])))
    print(json.dumps(dict(valid_control=True,rejected=cases,count=len(cases))))
if __name__=='__main__':main()
