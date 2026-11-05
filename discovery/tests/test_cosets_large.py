from pathlib import Path
import unittest
from discovery.core import RunContext
from discovery.engines.cosets import build

class LargeSyndromeTests(unittest.TestCase):
 def test_eighteen_bit_quotient(self):
  H=[[int(i==j) for j in range(18)] for i in range(18)]
  a=build([] ,{'parity_check':H},RunContext(Path('.')))
  self.assertEqual(a.facts['rank'],18);self.assertEqual(a.data['proof']['syndrome_states'],1<<18);self.assertEqual(a.facts['cardinality'],1)
 def test_nineteen_bit_quotient(self):
  H=[[int(i==j) for j in range(19)] for i in range(19)]
  a=build([] ,{'parity_check':H},RunContext(Path('.'),{'max_codewords':1<<19}));self.assertEqual(a.facts['rank'],19)
 def test_twenty_bit_rejected(self):
  H=[[int(i==j) for j in range(20)] for i in range(20)]
  with self.assertRaises(ValueError):build([] ,{'parity_check':H},RunContext(Path('.'),{'max_codewords':1<<20}))
if __name__=='__main__':unittest.main()
