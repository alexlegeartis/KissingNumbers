"""Coverage claims must fail closed when source, inventory or provenance is altered."""
import sys
sys.dont_write_bytecode=True
import copy,importlib.util,json,pathlib,unittest
HERE=pathlib.Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('coverage_builder',HERE/'build_matrix.py')
builder=importlib.util.module_from_spec(spec);spec.loader.exec_module(builder)
class CoverageChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.matrix=json.loads((HERE/'matrix.json').read_text())
    def test_saved_matrix(self):
        checked=builder.check(copy.deepcopy(self.matrix))
        self.assertEqual([r['dimension'] for r in checked['rows']],list(range(1,25)))
    def test_missing_dimension_is_rejected(self):
        bad=copy.deepcopy(self.matrix);bad['rows'].pop()
        with self.assertRaisesRegex(ValueError,'exactly24'):builder.check(bad)
    def test_changed_target_is_rejected(self):
        bad=copy.deepcopy(self.matrix);bad['rows'][18]['sourced_target']=11948
        with self.assertRaisesRegex(ValueError,'Source mismatch'):builder.check(bad)
    def test_generated_control_cannot_be_labeled_search_discovered(self):
        bad=copy.deepcopy(self.matrix);row=next(r for r in bad['rows'] if r['integrated_runs']);row['integrated_runs'][0]['origin']='search-discovered'
        with self.assertRaisesRegex(ValueError,'provenance mismatch'):builder.check(bad)
    def test_artifact_cannot_be_moved_to_another_dimension(self):
        bad=copy.deepcopy(self.matrix);bad['rows'][1]['integrated_runs']=[copy.deepcopy(bad['rows'][0]['integrated_runs'][0])]
        with self.assertRaisesRegex(ValueError,'enclosing row'):builder.check(bad)
    def test_existing_integrated_recipe_cannot_be_hidden(self):
        bad=copy.deepcopy(self.matrix);row=bad['rows'][0];row.update(integrated_runs=[],implemented_recipe_or_none=None,coverage_status='none',imported_generated_search_discovered=[],independent_verification_status='no integrated exact target proof')
        with self.assertRaisesRegex(ValueError,'cover audited eligible inventory'):builder.check(bad)
    def test_imported_ancestry_survives_composition(self):
        program=json.loads((builder.ROOT/'discovery/recipes/orthogonal_union.json').read_text())
        self.assertEqual(builder.recipe_origin(program,program['outputs'][0]),'imported')
    def test_mixed_ancestry_is_distinguished(self):
        program={'nodes':[{'id':'coordinates','op':'geometry.integer_code'},{'id':'formula','op':'lattice.root_system'},{'id':'union','op':'geometry.orthogonal_union','inputs':['coordinates','formula']}]}
        self.assertEqual(builder.recipe_origin(program,'union'),'mixed')
    def test_unaudited_new_recipe_inventory_is_rejected(self):
        bad=copy.deepcopy(self.matrix);bad['audited_eligible_recipe_paths']=[]
        with self.assertRaisesRegex(ValueError,'inventory changed'):builder.check(bad)
if __name__=='__main__':unittest.main()
