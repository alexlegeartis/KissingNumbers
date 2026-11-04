from fractions import Fraction
from pathlib import Path
import tempfile
import unittest

from discovery.benchmarks.blind_rediscovery import read_table,growth_diagnostic,cycle_catalogue


class BlindBenchmarkTests(unittest.TestCase):
    def test_colspan_is_an_exact_upper_bound_not_ratio_one(self):
        rows=['<tr><td>1</td><td colspan="2">2</td><td>1</td><td></td></tr>']
        rows.extend(f'<tr><td>{n}</td><td>{n*3}</td><td>{n*6}</td><td>2</td></tr>' for n in range(2,49))
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'table.html';path.write_text('<table>'+''.join(rows)+'</table>')
            table=read_table(path)
        self.assertEqual(table[1],{'lower':2,'upper':2,'proved_exact':True})
        self.assertEqual(table[2]['upper'],12)
        diagnostic={r['dimension']:r for r in growth_diagnostic(table)}
        self.assertEqual(diagnostic[2]['gap_expansion']['decimal'],2.0)
        self.assertEqual(diagnostic[3]['gap_expansion']['decimal'],1.0)
        self.assertFalse(diagnostic[2]['new_record_evidence'])

    def test_catalogue_uses_only_partitions_of_the_dimension(self):
        for n in range(1,25):
            choices=cycle_catalogue(n)
            self.assertEqual(len(choices),len(set(map(tuple,choices))))
            self.assertTrue(all(sum(c)==n and min(c)>0 for c in choices))


if __name__=='__main__':unittest.main()
