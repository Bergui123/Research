"""Regression checks for the sharper benchmark certificate and new stencils."""
from fractions import Fraction
import json
from pathlib import Path
import subprocess
import sys
import unittest

ROOT=Path(__file__).resolve().parents[2]
VAL=ROOT/'review/validation'


class SharperB2Test(unittest.TestCase):
    def test_majorant_and_expansion_for_new_selection(self):
        majorant=json.loads((VAL/'b2-polynomial-majorant.json').read_text())
        self.assertLess(Fraction(majorant['B2_upper_rational']),56)
        output=subprocess.check_output([sys.executable,str(VAL/'check_expansion_enclosures.py'),
                                        '--case-set','b2-56'],text=True)
        cases=json.loads(output)['cases']
        self.assertEqual([(c['target'],Fraction(c['dx']),c['h']) for c in cases],
                         [('0.05',Fraction(1,2),8),('0.01',Fraction(1,4),16),
                          ('0.001',Fraction(2,15),16)])
        for c in cases:
            for node in c['nodes']:
                self.assertLess(Fraction(node['lower']),Fraction(node['upper']))

    def test_recost_rejects_a_bound_for_a_different_stencil(self):
        result=subprocess.run([sys.executable,str(VAL/'check_projected_recost.py'),
                               '--bounds',str(VAL/'projected-b2-56.json')],
                              capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('AssertionError',result.stderr)


if __name__=='__main__':unittest.main()
