"""Exact Gram, weighted section, scalar budget and corruption controls."""
from copy import deepcopy
import contextlib
import io
import tempfile
from fractions import Fraction
from pathlib import Path
import json
import subprocess
import sys
import unittest
from unittest.mock import patch

from discovery.core import Artifact, RunContext, Verification
from discovery.certification import verify_artifact, verify_integer_spherical
from discovery.exact_gram import bounded_fraction, exact_integer_rank, verify_rational_gram
from discovery.engines import geometry
from discovery.runner import run_program, output_artifact, registry


class GramTests(unittest.TestCase):
    def setUp(self):self.ctx=RunContext(Path('.'))
    def triangle(self):return {'gram':[[1,[-1,2],[-1,2]],[[-1,2],1,[-1,2]],[[-1,2],[-1,2],1]]}
    def test_rational_singular_gram_exact_rank(self):
        result=verify_rational_gram(self.triangle(),self.ctx)
        self.assertEqual(result.status,'exact')
        self.assertEqual(result.details['facts']['dimension'],2)
        self.assertEqual(result.details['facts']['max_inner_product'],[-1,2])
        self.assertEqual(result.details['facts']['pairs_checked'],3)
        one=verify_rational_gram({'gram':[[[2,3]]],'norm_squared':[2,3]},self.ctx)
        self.assertEqual(one.details['facts']['dimension'],1)
    def test_reject_false_geometry(self):
        cases=[{'gram':[[1,-2],[-2,1]]},
               {'gram':[[1,-1,-1],[-1,1,[1,2]],[-1,[1,2],1]]},
               {'gram':[[1,1],[1,1]]}, {'gram':[[1,[3,4]],[[3,4],1]]},
               {'gram':[[1,0],[[-1,2],1]]}, {'gram':[[1,0],[0,2]]},
               {'gram':[[1,0],[0,1]],'norm_squared':2}, {'gram':[[0]]},
               {'gram':[[1.0]]}, {'gram':[[True]]}, {'gram':[[[1,0]]]},
               {'gram':[[[1,-2]]]}, {'gram':[[1,0]]}, {'gram':[]},
               {'gram':[[1]],'trust_saved_rank':1}]
        for data in cases:
            with self.subTest(data=data),self.assertRaises(ValueError):verify_rational_gram(data,self.ctx)
    def test_facts_and_status_are_not_authority(self):
        a=geometry.gram_code([],self.triangle(),self.ctx)
        self.assertEqual(verify_artifact(a,self.ctx).status,'exact')
        for key,value in [('dimension',3),('rank_exact',3),('cardinality',999),('norm_squared',2),('new_sota',True)]:
            b=deepcopy(a);b.facts[key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):verify_artifact(b,self.ctx)
        a.data['gram'][0][1]=1
        with self.assertRaises(ValueError):verify_artifact(a,self.ctx)
    def test_budget_precedes_parsing_and_arithmetic(self):
        for limit in [{'max_gram_entries':8},{'max_pairs':2},{'max_materialized_points':2}]:
            with patch('discovery.exact_gram.bounded_fraction') as parser:
                with self.assertRaises(ValueError):verify_rational_gram(self.triangle(),RunContext(Path('.'),limit))
                parser.assert_not_called()
        with patch('discovery.exact_gram.Fraction') as constructor:
            with self.assertRaises(ValueError):bounded_fraction([1,2**50],RunContext(Path('.'),{'max_rational_bits':8}))
            constructor.assert_not_called()
        with self.assertRaises(ValueError):verify_rational_gram(self.triangle(),RunContext(Path('.'),{'max_rank_work':1}))
        with self.assertRaises(ValueError):verify_rational_gram(self.triangle(),RunContext(Path('.'),{'max_dimension':1}))
        with self.assertRaises(ValueError):verify_integer_spherical({'coordinates':[[2**50]],'norm_squared':2**100},RunContext(Path('.'),{'max_rational_bits':8}))
    def test_fraction_free_rank_handles_out_of_order_pivots(self):
        self.assertEqual(exact_integer_rank([[0,2,0],[1,1,3],[0,0,1]],self.ctx)[0],3)
        self.assertEqual(exact_integer_rank([[0,2,1],[1,1,0],[2,0,-1]],self.ctx)[0],2)
        self.assertEqual(exact_integer_rank([[0,2,1],[1,1,0],[1,0,-1]],self.ctx)[0],3)
        # Independent rational reference on diverse rectangular integer inputs.
        import numpy as np
        import sympy as sp
        rng=np.random.default_rng(341)
        for height,width in [(3,5),(7,4),(1,8),(9,9)]:
            for _ in range(8):
                rows=rng.integers(-4,5,size=(height,width)).tolist()
                self.assertEqual(exact_integer_rank(rows,self.ctx)[0],int(sp.Matrix(rows).rank()))
        # Deliberately dependent inputs and shuffled nonunit pivots catch
        # incorrect partial-row scaling that generic full-rank samples miss.
        for height,width in [(7,5),(5,8),(8,8)]:
            for factor_rank in range(1,min(height,width)+1):
                left=rng.integers(-3,4,size=(height,factor_rank))
                right=rng.integers(-3,4,size=(factor_rank,width))
                matrix=left@right
                matrix*=rng.integers(1,5,size=(height,1))
                matrix=matrix[rng.permutation(height)][:,rng.permutation(width)]
                rows=matrix.tolist()
                self.assertEqual(exact_integer_rank(rows,self.ctx)[0],int(sp.Matrix(rows).rank()))
        with self.assertRaises(ValueError):exact_integer_rank([[1,1],[1,-1]],RunContext(Path('.'),{'max_rank_work':1}))
    def test_full_rank_still_validates_remaining_rows(self):
        for bad in ([True,0],[1],[2**100,0]):
            with self.assertRaises(ValueError):
                exact_integer_rank([[1,0],[0,1],bad],RunContext(Path('.'),{'max_rational_bits':32}))
        self.assertEqual(exact_integer_rank([[1,0],[0,1]]+[[2,3]]*100,self.ctx)[0],2)
    def test_optimized_python_retains_rejection(self):
        code="from discovery.exact_gram import verify_rational_gram; verify_rational_gram({'gram':[[1,-2],[-2,1]]})"
        result=subprocess.run([sys.executable,'-O','-c',code],cwd=Path(__file__).resolve().parents[2],capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0);self.assertIn('positive semidefinite',result.stderr)
    def test_registered_dag_and_symbolic_status(self):
        self.assertIn('geometry.gram_code',registry()[0])
        p={'schema_version':1,'nodes':[{'id':'triangle','op':'geometry.gram_code','params':self.triangle()}],'outputs':['triangle']}
        a=output_artifact(run_program(p,self.ctx))
        self.assertEqual(a.kind,'gram_spherical_code');self.assertEqual(a.facts['dimension'],2)
        self.assertIn('symbolic',a.facts['dimension_scope'])

    def test_cli_compare_and_scheduler_include_exact_grams(self):
        from discovery.__main__ import main
        from discovery.scheduler import search_grammar
        p={'schema_version':1,'nodes':[{'id':'triangle','op':'geometry.gram_code','params':self.triangle()}],'outputs':['triangle']}
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);receipt=root/'receipt.json';catalog=root/'catalog.json'
            receipt.write_text(json.dumps(run_program(p,self.ctx)))
            catalog.write_text(json.dumps({'comparisons':[{'dimension':2,'cardinality':6,'source':'test fixture'}]}))
            output=io.StringIO()
            with contextlib.redirect_stdout(output):
                code=main(['compare',str(receipt),'--catalog',str(catalog)])
            self.assertEqual(code,0);self.assertEqual(json.loads(output.getvalue())['cardinality'],3)
            grammar={'schema_version':1,'slots':[{'id':'code','choices':[{'op':'geometry.gram_code','params':self.triangle()}]}],'outputs':['code']}
            result=search_grammar(grammar,self.ctx,root/'search',max_programs=1,time_limit_seconds=15)
            self.assertEqual(result['best_exact_by_dimension']['2']['cardinality'],3)


class SectionTests(unittest.TestCase):
    def setUp(self):
        self.ctx=RunContext(Path('.'))
        self.source=geometry.integer_code([],{'coordinates':[[1,1],[-1,-1],[1,-1],[-1,1]],'weights':[1,3],'norm_squared':4},self.ctx)
    def test_weighted_rational_normal_and_exact_span(self):
        result=geometry.section([self.source],{'normal':[[3,2],[1,2]]},self.ctx)
        self.assertEqual(result.data['coordinates'],[[1,-1],[-1,1]])
        self.assertEqual(result.data['weights'],[1,3])
        self.assertEqual((result.facts['dimension'],result.facts['ambient_dimension'],result.facts['cardinality']),(1,2,2))
        self.assertEqual(verify_artifact(result,self.ctx).status,'exact')
        source=geometry.integer_code([],{'coordinates':[[1,0],[-1,0]],'weights':[1,3],'norm_squared':1},self.ctx)
        same=geometry.section([source],{'normal':[0,1]},self.ctx)
        self.assertEqual(same.facts['dimension'],source.facts['dimension'])
    def test_section_corruption_rejected(self):
        a=geometry.section([self.source],{'normal':[3,1]},self.ctx)
        for mutation in ('normal','rows','metric','source','facts','float'):
            b=deepcopy(a)
            if mutation=='normal':b.data['section']['normal']=[1,0]
            elif mutation=='rows':b.data['coordinates'].pop()
            elif mutation=='metric':b.data['weights']=[3,1]
            elif mutation=='source':b.data['section']['source']['coordinates'][0]=[9,9]
            elif mutation=='facts':b.facts['dimension']=2
            else:b.data['section']['normal']=[3.0,1]
            with self.subTest(mutation=mutation),self.assertRaises(ValueError):verify_artifact(b,self.ctx)
    def test_invalid_normal_source_and_budget(self):
        for normal in ([0,0],[1],[1,0],[1.0,1],[[1,0],1]):
            with self.subTest(normal=normal),self.assertRaises(ValueError):geometry.section([self.source],{'normal':normal},self.ctx)
        with patch('discovery.certification.verify_artifact') as verify:
            with self.assertRaises(ValueError):geometry.section([self.source],{'normal':[1]},self.ctx)
            verify.assert_not_called()
        fake=Artifact('spherical_code',{'real_coordinates':[[1.,0.]]},verification=Verification('exact'))
        with self.assertRaises(ValueError):geometry.section([fake],{'normal':[0,1]},self.ctx)
        bad=deepcopy(self.source);bad.facts['cardinality']=999
        with self.assertRaises(ValueError):geometry.section([bad],{'normal':[3,1]},self.ctx)
        with self.assertRaises(ValueError):geometry.section([self.source],{'normal':[3,1]},RunContext(Path('.'),{'max_pairs':1}))
    def test_section_budget_precedes_normal_fraction_allocation(self):
        oversized=Artifact('spherical_code',{'coordinates':[[1,0,0],[-1,0,0]],'weights':[1,1,1],'norm_squared':1})
        stored=Artifact('spherical_code',{'coordinates':[[1,0,0],[-1,0,0]],'weights':[1,1,1],'norm_squared':1,
                 'section':{'source':oversized.data,'normal':[0,1,0]}})
        ctx=RunContext(Path('.'),{'max_dimension':2})
        with patch('discovery.exact_gram.bounded_fraction') as parser:
            with self.assertRaises(ValueError):geometry.section([oversized],{'normal':[0,1,0]},ctx)
            parser.assert_not_called()
            with self.assertRaises(ValueError):verify_artifact(stored,ctx)
            parser.assert_not_called()
    def test_section_dag_replays_membership(self):
        p={'schema_version':1,'nodes':[{'id':'source','op':'geometry.integer_code','params':self.source.data},
           {'id':'cut','op':'geometry.section','inputs':['source'],'params':{'normal':[3,1]}}],'outputs':['cut']}
        a=output_artifact(run_program(p,self.ctx));self.assertEqual(a.facts['cardinality'],2)

if __name__=='__main__':unittest.main()
