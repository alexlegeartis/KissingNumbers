"""Actual small assembly generation, metric oracle witnesses and resource gates."""
from dataclasses import replace, asdict
from itertools import product
from pathlib import Path
import json
import time
import tempfile
import subprocess
import sys
from copy import deepcopy
from unittest.mock import patch

import numpy as np
import unittest

from discovery.core import RunContext
from discovery.engines.coset_search import (CosetAssemblyRequest, search_gray_coset_assemblies,
    coset_assembly, verify_artifact, _plans, _clique, _SearchCutoff)


def weight(v,q=8):
    return sum(0 if x==0 else q//2 if x==q//2 else q//4 for x in v)


class CosetSearchTests(unittest.TestCase):
    def fixture(self):
        return CosetAssemblyRequest(((1,0,1),(0,1,1)),((1,0),(0,1)),8,4,5,
            seed=4,offsets_per_family=64,information_sets=4,refinement_sets=8,
            max_vertices=32,max_cliques=16,families=('uniform',))

    def test_actual_five_coset_generation_and_independent_small_geometry(self):
        result=search_gray_coset_assemblies(self.fixture())
        self.assertEqual(result.forecast['total_cardinality_if_certified'],80)
        self.assertEqual(result.statistics['quotient_distinct_offsets'],31)
        self.assertTrue(result.candidates)
        candidate=result.candidates[0]
        self.assertEqual(candidate['distance_status'],'UNVERIFIED')
        self.assertNotIn('facts',result.to_dict())
        json.dumps(result.to_dict())
        # This independent finite check is a test fixture, not a production
        # promotion of sampled upper estimates to a distance certificate.
        points=set()
        for v in candidate['offsets']:
            for a,b in product(range(0,8,2),repeat=2):
                points.add(tuple((v[j]+x)%8 for j,x in enumerate((a,b,(a+b)%8))))
        self.assertEqual(len(points),80)
        points=sorted(points)
        self.assertEqual(min(weight([(a-b)%8 for a,b in zip(x,y)])
            for i,x in enumerate(points)for y in points[i+1:]),4)

    def test_computed_torsion_guard_avoids_insufficient_family(self):
        request=CosetAssemblyRequest(((1,1),),((1,),),8,5,3,
            seed=1,offsets_per_family=32,information_sets=2,refinement_sets=2,
            max_vertices=8,max_cliques=8,families=('torsion',))
        result=search_gray_coset_assemblies(request)
        self.assertEqual(result.torsion_guard['binary_covering_radius'],1)
        self.assertTrue(result.torsion_guard['active'])
        self.assertGreater(result.torsion_guard['prototypes_removed'],0)
        self.assertEqual(result.torsion_guard['maximum_compatible_cosets_in_that_family'],2)
        self.assertFalse(result.candidates)

    def test_invalid_rank_and_resource_requests_rejected_before_search(self):
        base=self.fixture()
        for request in (replace(base,parity_checks=((1,0),(1,0))),
                        replace(base,modulus=True),replace(base,target_distance=13),
                        replace(base,required_cosets=34),replace(base,seed=True),
                        replace(base,information_sets=10**8)):
            with self.assertRaises(ValueError):search_gray_coset_assemblies(request)
        with self.assertRaises(ValueError):
            search_gray_coset_assemblies(base,RunContext(Path('.'),{'max_pairs':1}))

    def test_forecast_companion_is_only_an_assumption(self):
        request=replace(self.fixture(),companion_cardinality_forecast=1000)
        result=search_gray_coset_assemblies(request)
        self.assertEqual(result.forecast['companion_cardinality_assumption'],1000)
        self.assertEqual(result.forecast['total_cardinality_if_certified'],1080)
        self.assertEqual(result.forecast['scope'],'cardinality forecast only; not verified spherical facts')
        self.assertTrue(result.status.startswith('unverified'))

    def test_mixed_requests_retain_torsion_vertices(self):
        request=CosetAssemblyRequest(((1,1),),((1,),),8,5,3,
            seed=1,offsets_per_family=32,information_sets=2,refinement_sets=2,
            max_vertices=8,max_cliques=8,families=('uniform','torsion'))
        result=search_gray_coset_assemblies(request)
        self.assertTrue(result.torsion_guard['active'])
        self.assertEqual(result.torsion_guard['prototypes_removed'],0)
        self.assertTrue(result.torsion_guard['mixed_requests_retain_torsion_vertices'])
        self.assertGreater(result.statistics['quotient_distinct_offsets'],1)

    def test_large_correction_space_is_sampled_without_population_array(self):
        G=np.eye(8,dtype=np.int64);original=np.arange
        def bounded_arange(*args,**kwargs):
            if args and isinstance(args[0],int)and args[0]>10000:
                raise AssertionError('large correction population allocated')
            return original(*args,**kwargs)
        with patch('numpy.arange',bounded_arange):
            plans=_plans(G,G,16,1,np.random.default_rng(1),64)
        self.assertEqual(len(plans[0][-1]),64)
        self.assertTrue(all(len(c)==8 for c in plans[0][-1]))

    def test_clique_deadline_and_branch_cutoff_are_explicit(self):
        adjacency=[((1<<12)-1)^(1<<i)for i in range(12)]
        with self.assertRaises(_SearchCutoff):
            _clique(adjacency,5,deadline=time.monotonic()-1)
        with self.assertRaises(_SearchCutoff):
            _clique(adjacency,5,max_branches=1)

    def test_registered_actual_generation_never_certifies_forecast(self):
        from discovery.runner import run_program,output_artifact,registry
        from discovery.certification import verify_artifact as dispatch_verify
        params=asdict(self.fixture())
        program={'schema_version':1,'name':'actual_coset_fixture','nodes':[
            {'id':'search','op':'search.coset_assembly','params':params}],
            'outputs':['search']}
        self.assertIn('search.coset_assembly',registry()[0])
        artifact=output_artifact(run_program(program,RunContext(Path('.'))))
        self.assertEqual(artifact.kind,'search_result')
        self.assertEqual(artifact.facts,{})
        self.assertTrue(artifact.data['result']['candidates'])
        self.assertEqual(dispatch_verify(artifact).status,'unverified')
        artifact.verification.status='exact'
        artifact.verification.details={'facts':{'dimension':12,'cardinality':80}}
        self.assertEqual(dispatch_verify(artifact).status,'unverified')
        damaged=deepcopy(artifact);damaged.facts={'dimension':12,'cardinality':80}
        with self.assertRaises(ValueError):dispatch_verify(damaged)
        damaged=deepcopy(artifact)
        damaged.data['result']['candidates'][0]['distance_status']='EXACT'
        with self.assertRaises(ValueError):dispatch_verify(damaged)
        with tempfile.TemporaryDirectory()as tmp:
            filename=Path(tmp)/'candidate.json';filename.write_text(json.dumps(artifact.to_dict()))
            outcome=subprocess.run([sys.executable,'-m','discovery','verify',str(filename)],
                capture_output=True,text=True,check=False)
        self.assertEqual(outcome.returncode,2,outcome.stdout+outcome.stderr)
        self.assertEqual(json.loads(outcome.stdout)['status'],'unverified')


if __name__=='__main__':unittest.main()
