"""Regression checks for sweep-specific expansion enclosures."""
import json
from pathlib import Path
import subprocess
import sys
import unittest
from fractions import Fraction as F

VAL=Path(__file__).resolve().parent


class SweepExpansionTest(unittest.TestCase):
    def run_case(self,node,cert,selection):
        return subprocess.run([sys.executable,str(VAL/'sweep_expansion.py'),'--nodes',str(VAL/node),
                               '--certificate',str(VAL/cert),'--selection',str(VAL/selection)],
                              capture_output=True,text=True)

    def test_selected_sweep_nodes_are_enclosed(self):
        result=self.run_case('sweep-nodes-rho0-lam1.json.gz','sweep-cert-rho0-lam1.json',
                             'sweep-reselection-rho0-lam1.json')
        self.assertEqual(result.returncode,0,result.stderr)
        data=json.loads(result.stdout)
        self.assertEqual((data['rho'],data['lambda'],data['n'],data['h']),('0','1',128,16))
        self.assertEqual(len(data['nodes']),4)
        for row in data['nodes']:self.assertLess(F(row['lower']),F(row['upper']))

    def test_mismatched_certificate_is_rejected(self):
        result=self.run_case('sweep-nodes-rho0-lam1.json.gz','sweep-cert-rho02-lam1.json',
                             'sweep-reselection-rho0-lam1.json')
        self.assertNotEqual(result.returncode,0)
        self.assertIn('AssertionError',result.stderr)

    def test_rational_grid_bound_checks_actual_selected_stencil(self):
        ident='rho0-lam1'
        result=subprocess.run([sys.executable,str(VAL/'sweep_grid_bound.py'),
            '--nodes',str(VAL/f'sweep-nodes-{ident}.json.gz'),
            '--certificate',str(VAL/f'sweep-cert-{ident}.json'),
            '--expansion',str(VAL/f'sweep-expansion-{ident}.json')],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        data=json.loads(result.stdout)
        self.assertEqual((data['rho'],data['lambda'],data['n'],data['h']),('0','1',128,16))
        self.assertTrue(data['E_below_S'])
        self.assertLess(F(data['relative_grid_bound_upper']),F(1,200))

    def test_recost_uses_matching_sweep_stencil(self):
        ident='rho0-lam1'
        result=subprocess.run([sys.executable,str(VAL/'sweep_recost.py'),
            '--bounds',str(VAL/f'sweep-grid-bound-{ident}.json'),
            '--selection',str(VAL/f'sweep-reselection-{ident}.json')],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        data=json.loads(result.stdout)
        self.assertIsNotNone(data['best'])
        self.assertEqual((data['best']['nu'],2**data['best']['n'],data['best']['h']),(8.,128,16))

    def test_recost_rejects_a_different_parameter_case(self):
        result=subprocess.run([sys.executable,str(VAL/'sweep_recost.py'),
            '--bounds',str(VAL/'sweep-grid-bound-rho0-lam1.json'),
            '--selection',str(VAL/'sweep-reselection-rho02-lam1.json')],capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('AssertionError',result.stderr)

    def test_grid_bound_rejects_mismatched_certificate(self):
        result=subprocess.run([sys.executable,str(VAL/'sweep_grid_bound.py'),
            '--nodes',str(VAL/'sweep-nodes-rho0-lam1.json.gz'),
            '--certificate',str(VAL/'sweep-cert-rho02-lam1.json'),
            '--expansion',str(VAL/'sweep-expansion-rho0-lam1.json')],capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('AssertionError',result.stderr)


if __name__=='__main__':unittest.main()
