import contextlib
from copy import deepcopy
import io
import json
from pathlib import Path
import tempfile
import unittest
from discovery.core import Artifact,Verification,RunContext,digest
from discovery.certification import verify_artifact,verify_integer_spherical
from discovery.runner import run_program,output_artifact,validate_program
from discovery.__main__ import main


def square():
    return {'id':'square','op':'geometry.integer_code','params':{'coordinates':[[1,0],[-1,0],[0,1],[0,-1]],'norm_squared':1}}

def program():return {'schema_version':1,'nodes':[square()],'outputs':['square']}

class CoreTests(unittest.TestCase):
    def setUp(self):self.ctx=RunContext(Path('.'))
    def test_dag_and_orthogonal_composition(self):
        p=program();p['nodes'].insert(0,{'id':'union','op':'geometry.orthogonal_union','inputs':['square','square']});p['outputs']=['union']
        r=run_program(p,self.ctx);a=output_artifact(r)
        self.assertEqual((a.facts['dimension'],a.facts['cardinality']),(4,8))
        self.assertEqual(verify_artifact(a,self.ctx).status,'exact')
    def test_bad_graphs_are_rejected(self):
        for change in ('cycle','duplicate','missing','unsupported','type'):
            p=program()
            if change=='cycle':p['nodes']=[{'id':'a','op':'geometry.orthogonal_union','inputs':['a','a']}];p['outputs']=['a']
            elif change=='duplicate':p['nodes'].append(square())
            elif change=='missing':p['outputs']=['absent']
            elif change=='unsupported':p['nodes'][0]['op']='os.system'
            elif change=='type':p['nodes'].append({'id':'bad','op':'geometry.sign_embed','inputs':['square']})
            with self.subTest(change=change),self.assertRaises(ValueError):run_program(p,self.ctx)
    def test_limits_precede_pair_work(self):
        with self.assertRaises(ValueError):run_program(program(),RunContext(Path('.'),{'max_pairs':1}))
        with self.assertRaises(ValueError):RunContext(Path('.'),{'unknown':5})
        with self.assertRaises(ValueError):RunContext(Path('.'),{'max_pairs':float('inf')})
    def test_invalid_coordinate_certificates(self):
        data=square()['params']
        for mutation in ('duplicate','norm','float','weights','pair'):
            bad=deepcopy(data)
            if mutation=='duplicate':bad['coordinates'][0]=bad['coordinates'][1]
            elif mutation=='norm':bad['norm_squared']=2
            elif mutation=='float':bad['coordinates'][0][0]=1.0
            elif mutation=='weights':bad['weights']=[1,-1]
            else:bad={'coordinates':[[3,4],[4,3]],'norm_squared':25}
            with self.subTest(mutation=mutation),self.assertRaises(ValueError):verify_integer_spherical(bad,self.ctx)
    def test_large_integer_backend_is_exact(self):
        scale=10**15
        data={'coordinates':[[scale,0],[-scale,0],[0,scale]],'norm_squared':scale**2}
        r=verify_integer_spherical(data,self.ctx)
        self.assertEqual(r.details['integer_backend'],'Python integers');self.assertEqual(r.details['dimension'],2)
    def test_bounded_binary32_and_int64_paths_match_scalar_geometry(self):
        for scale,backend in [(1,'bounded exact binary32'),(10000,'int64 proven overflow-safe')]:
            rows=[[3*scale,4*scale],[-5*scale,0],[4*scale,-3*scale]]
            result=verify_integer_spherical({'coordinates':rows,'norm_squared':25*scale**2},self.ctx)
            products=[sum(x*y for x,y in zip(a,b)) for i,a in enumerate(rows) for b in rows[i+1:]]
            self.assertEqual(result.details['integer_backend'],backend)
            self.assertEqual(result.details['max_inner_product'],max(products))
        with self.assertRaisesRegex(ValueError,'max_gram_entries'):
            verify_integer_spherical(square()['params'],RunContext(Path('.'),{'max_gram_entries':15}))
    def test_status_and_facts_cannot_forge_proof(self):
        a=output_artifact(run_program(program(),self.ctx));a.facts['cardinality']=1000
        with self.assertRaises(ValueError):verify_artifact(a,self.ctx)
        a.facts['cardinality']=4;a.data['coordinates'][0]=[2,0]
        with self.assertRaises(ValueError):verify_artifact(a,self.ctx)
    def test_forged_auxiliary_facts_and_receipt_hashes(self):
        r=run_program(program(),self.ctx)
        for key,value in [('norm_squared',999),('max_inner_product',999),('dimension',True)]:
            a=output_artifact(deepcopy(r));a.facts[key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):verify_artifact(a,self.ctx)
        for field in ('program_sha256','event'):
            bad=deepcopy(r)
            if field=='event':bad['events'][0]['artifact_sha256']='corrupt'
            else:bad[field]='corrupt'
            with self.subTest(field=field),self.assertRaises(ValueError):output_artifact(bad)
    def test_uncertified_geometry_and_novelty_claims_rejected(self):
        for key,value in [('max_inner_product_upper_bound',-1),
                          ('dimension_scope','exact rank 999'),('rank_exact',999),
                          ('global_optimality_claim',True),('new_sota',True),
                          ('worldwide_priority_claim',True),('new_sota',0)]:
            a=output_artifact(run_program(program(),self.ctx));a.facts[key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):verify_artifact(a,self.ctx)
        a=output_artifact(run_program(program(),self.ctx));a.facts['new_sota']=False
        self.assertEqual(verify_artifact(a,self.ctx).status,'exact')
    def test_serialization_and_digest(self):
        p=program();self.assertEqual(digest(p),digest(json.loads(json.dumps(p))))
        with self.assertRaises(ValueError):Artifact.from_dict({'kind':'bad','data':{'x':float('nan')}})
    def test_cli_roundtrip_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as d,contextlib.redirect_stdout(io.StringIO()),contextlib.redirect_stderr(io.StringIO()):
            source=Path(d)/'recipe.json';source.write_text(json.dumps(program()));out=Path(d)/'run.json'
            self.assertEqual(main(['run',str(source),'--out',str(out)]),0)
            self.assertEqual(main(['verify',str(out)]),0)
            self.assertEqual(main(['run',str(source),'--out',str(out)]),1)
            r=json.loads(out.read_text());r['artifacts']['square']['facts']['dimension']=100
            out.write_text(json.dumps(r));self.assertEqual(main(['verify',str(out)]),1)
    def test_implicit_gray_sign_embedding(self):
        p={'schema_version':1,'nodes':[
          {'id':'f','op':'rings.binary_factor','params':{'n':7,'index':0}},
          {'id':'h','op':'rings.hensel','inputs':['f'],'params':{'modulus':4}},
          {'id':'c','op':'rings.extend','inputs':['h']},
          {'id':'b','op':'rings.gray','inputs':['c']},
          {'id':'s','op':'geometry.sign_embed','inputs':['b']}],'outputs':['s']}
        a=output_artifact(run_program(p,self.ctx));self.assertEqual(a.kind,'implicit_spherical_code')
        self.assertEqual((a.facts['dimension'],a.facts['cardinality']),(16,256))
        self.assertFalse(a.verification.details['materialized'])
        self.assertIn('ambient',a.facts['dimension_scope'])
        a.data['binary_code']['facts']['min_distance']=99
        with self.assertRaises(ValueError):verify_artifact(a,self.ctx)

if __name__=='__main__':unittest.main()
