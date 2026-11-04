from copy import deepcopy
from pathlib import Path
from unittest.mock import patch
import unittest

from discovery.core import Artifact,RunContext
from discovery.certification import verify_artifact
from discovery.engines.exact_orbits import integer_orbit
from discovery.engines.geometry import integer_code,section


class CertifiedSectionTests(unittest.TestCase):
    def source(self):
        cycle=[[int(j==(i+1)%4) for j in range(4)] for i in range(4)]
        flip=[[int(i==j)*(-1 if i==0 else 1) for j in range(4)] for i in range(4)]
        permutation=[1,0,2,3]
        swap=[[int(j==permutation[i]) for j in range(4)] for i in range(4)]
        ctx=RunContext(Path('.'),{'max_pairs':23})
        a=integer_orbit([],{'seed':[1,1,0,0],'generators':[{'matrix':cycle,'denominator':1},{'matrix':flip,'denominator':1},{'matrix':swap,'denominator':1}]},ctx)
        return a,ctx

    def test_inherited_separation_exceeds_explicit_pair_budget(self):
        a,ctx=self.source();b=section([a],{'normal':[1,0,0,0]},ctx)
        self.assertEqual((b.facts['cardinality'],b.facts['rank_exact']),(12,3))
        r=verify_artifact(b,ctx)
        self.assertEqual(r.details['unordered_output_pairs_covered'],66)
        self.assertTrue(r.details['source_separation_inherited'])

    def test_batch_and_chained_normals(self):
        a,ctx=self.source()
        b=section([a],{'normals':[[1,0,0,0],[0,1,0,0]]},ctx)
        c=section([section([a],{'normal':[1,0,0,0]},ctx)],{'normal':[0,1,0,0]},ctx)
        self.assertEqual((b.facts['cardinality'],b.facts['rank_exact']),(4,2))
        self.assertEqual(b.data['coordinates'],c.data['coordinates'])

    def test_weighted_metric_and_rational_normal(self):
        ctx=RunContext(Path('.'))
        a=integer_code([],{'coordinates':[[2,0],[1,1],[-1,1],[-2,0],[-1,-1],[1,-1]],'weights':[1,3],'norm_squared':4},ctx)
        b=section([a],{'normals':[[1,[1,3]]]},ctx)
        self.assertEqual(b.data['coordinates'],[[-1,1],[1,-1]])
        self.assertEqual(b.facts['rank_exact'],1)

    def test_tampered_membership_metric_norm_and_source_reject(self):
        a,ctx=self.source();b=section([a],{'normal':[1,0,0,0]},ctx)
        for change in ('point','weights','norm','source'):
            c=Artifact.from_dict(b.to_dict());c.facts={}
            if change=='point':c.data['coordinates'][0]=[1,1,0,0]
            elif change=='weights':c.data['weights'][0]=2
            elif change=='norm':c.data['norm_squared']=3
            else:c.data['source']['data']['coordinates'].pop()
            with self.assertRaises(ValueError):verify_artifact(c,ctx)

    def test_malformed_normal_rejects_before_source_replay(self):
        a,ctx=self.source()
        with patch('discovery.certification.verify_artifact',side_effect=AssertionError('source replayed')):
            with self.assertRaisesRegex(ValueError,'shaped normals'):
                section([a],{'normals':[[1,0]]},ctx)
            with self.assertRaisesRegex(ValueError,'exactly one'):
                section([a],{'normal':[1,0,0,0],'normals':[[1,0,0,0]]},ctx)

    def test_boolean_coordinate_or_metric_cannot_equal_integer_schema(self):
        a,ctx=self.source();b=section([a],{'normal':[1,0,0,0]},ctx)
        for field in ('coordinate','weight'):
            c=Artifact.from_dict(b.to_dict());c.facts={}
            if field=='weight':c.data['weights'][0]=True
            else:
                row=next(row for row in c.data['coordinates'] if 1 in row)
                row[row.index(1)]=True
            with self.assertRaises(ValueError):verify_artifact(c,ctx)


if __name__=='__main__':unittest.main()
