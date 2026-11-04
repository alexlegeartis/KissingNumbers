from pathlib import Path
import random
import unittest

from discovery.core import Artifact,RunContext
from discovery.engines.geometry import quadratic_code
from discovery.engines.signed_layers import compatible_dot,verify_artifact


class WeightedQuadraticTests(unittest.TestCase):
    def test_metric_changes_normalized_compatibility(self):
        data={'radicand':3,'rational_coordinates':[[1,0],[1,1]],'radical_coordinates':[[0,0],[0,0]],'weights':[1,4]}
        a=quadratic_code([],data,RunContext(Path('.')))
        self.assertEqual(a.verification.status,'exact')
        data['weights']=[4,1]
        with self.assertRaises(ValueError):quadratic_code([],data,RunContext(Path('.')))

    def test_vectorized_weighted_pairs_match_scalar_integer_oracle(self):
        rng=random.Random(20261104)
        for trial in range(150):
            weights=[rng.randrange(1,10) for _ in range(4)];A=[];B=[]
            for i in range(2):
                a=[];b=[]
                for j in range(4):
                    value=rng.randrange(-3,4)
                    if rng.randrange(2):a.append(value);b.append(0)
                    else:a.append(0);b.append(value)
                if not any(a+b):a[0]=1
                A.append(a);B.append(b)
            norms=[sum(w*(x*x+3*y*y) for w,x,y in zip(weights,a,b)) for a,b in zip(A,B)]
            p=sum(w*(x*y+3*u*v) for w,x,y,u,v in zip(weights,A[0],A[1],B[0],B[1]))
            q=sum(w*(x*v+u*y) for w,x,y,u,v in zip(weights,A[0],A[1],B[0],B[1]))
            expected=compatible_dot(p,q,norms[0]*norms[1])
            artifact=Artifact('quadratic_spherical_code',{'radicand':3,'rational_coordinates':A,'radical_coordinates':B,'weights':weights})
            if expected:self.assertEqual(verify_artifact(artifact).status,'exact')
            else:
                with self.assertRaises(ValueError):verify_artifact(artifact)

    def test_second_squaring_and_arbitrary_precision_fallback(self):
        for scale in (2**20,2**80):
            a=Artifact('quadratic_spherical_code',{'radicand':3,'rational_coordinates':[[scale,0],[scale,scale]],'radical_coordinates':[[0,0],[0,0]]})
            with self.assertRaises(ValueError):verify_artifact(a)
            a.data['rational_coordinates'][1]=[0,scale]
            self.assertEqual(verify_artifact(a).status,'exact')
        a.data['weights']=[1,0]
        with self.assertRaises(ValueError):verify_artifact(a)


if __name__=='__main__':unittest.main()
