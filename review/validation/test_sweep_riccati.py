"""Checks for the parameterized exact Riccati residual verifier."""
import copy
from fractions import Fraction as F
import unittest

from sweep_riccati import verify


class SweepRiccatiTest(unittest.TestCase):
    def fixture(self):
        # P(t)=2I solves P'=2I-P^2/2 exactly, for lambda=1 and rho=0.
        q=10**6
        return {'segments':1,'denominator':q,'rho':'0','lambda':'1',
                'components':['P11','P12','P22'],
                'nodes':[[2*q,0,2*q] for _ in range(3)]}

    def test_exact_constant_solution(self):
        result=verify(self.fixture(),min_hat=F(1),expected_initial=(F(2),F(0),F(2)))
        self.assertEqual(F(result['residual_integral_upper']),0)
        self.assertEqual(F(result['true_row_upper']),2)
        self.assertEqual(F(result['true_min_eigen_lower']),1)

    def test_invalid_initial_or_nonpositive_piece_rejected(self):
        data=self.fixture()
        with self.assertRaises(AssertionError):
            verify(data,min_hat=F(1))
        bad=copy.deepcopy(data)
        bad['nodes'][1][0]=-10**6
        with self.assertRaises(AssertionError):
            verify(bad,min_hat=F(1),expected_initial=(F(2),F(0),F(2)))


if __name__=='__main__':unittest.main()
