"""Corruption controls for the arithmetic certificate checker."""
import copy
import itertools
import json
import unittest
import numpy as np
from verify import ROOT, load_source, verify_math


class CertificateControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.A, cls.B = load_source(ROOT/'data/source26.npz', 26)
        with np.load(ROOT/'data/original26.npz', allow_pickle=False) as z:
            cls.C = z['candidates']
            cls.ptr, cls.idx = z['blocker_indptr'], z['blocker_indices']
        cls.proof = json.loads((ROOT/'data/original26_dual.json').read_text())

    def check_rejected(self, proof, message, candidates=None):
        with self.assertRaisesRegex(ValueError, message):
            verify_math(self.C if candidates is None else candidates,
                        self.A, self.B, proof)

    def test_valid_certificate(self):
        self.assertEqual(verify_math(self.C, self.A, self.B, self.proof)
                         ['integer_gain_upper_bound'], 0)

    def test_changed_norm(self):
        points = self.C.copy()
        points[0, 0] += 1
        self.check_rejected(self.proof, 'candidate norms', points)

    def test_negative_weight(self):
        proof = copy.deepcopy(self.proof)
        proof['cliques'][0]['weight'] = [-1, 1]
        self.check_rejected(proof, 'nonnegative clique weight')

    def test_false_owner(self):
        proof = copy.deepcopy(self.proof)
        row = next(r for r in proof['cliques'] if r['owner'] >= 0)
        dots = self.A @ self.C[row['members'][0]].astype(np.int64)
        row['owner'] = int(np.flatnonzero((dots <= 64) & ~np.any(self.B, axis=1))[0])
        self.check_rejected(proof, 'source owner incompatibility')

    def test_compatible_points_are_not_a_conflict_clique(self):
        incident = {}
        for i in range(len(self.C)):
            for owner in self.idx[self.ptr[i]:self.ptr[i+1]]:
                incident.setdefault(int(owner), []).append(i)
        for owner, members in incident.items():
            for i, j in itertools.combinations(members, 2):
                if sum(int(x)*int(y) for x, y in zip(self.C[i], self.C[j])) <= 64:
                    proof = copy.deepcopy(self.proof)
                    proof['cliques'].insert(0, dict(owner=owner, members=[i, j], weight=[1, 1]))
                    self.check_rejected(proof, 'candidate conflict clique')
                    return
        self.fail('No compatible points with a shared source blocker found')

    def test_missing_coverage(self):
        proof = copy.deepcopy(self.proof)
        proof['cliques'] = []
        proof['variable_bounds'] = []
        self.check_rejected(proof, 'candidate indicators covered')

    def test_wrong_declared_bound(self):
        proof = copy.deepcopy(self.proof)
        proof['exact_gain_upper_bound'] = [0, 1]
        self.check_rejected(proof, 'exact weighted bound')


if __name__ == '__main__':
    unittest.main()
