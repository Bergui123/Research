import math
import unittest
from fractions import Fraction as F

import numpy as np
from scipy.linalg import expm

from check_projected_box_bound import (ceil_scaled, log_interval, response_upper,
                                       sqrt_down, sqrt_up)


class ProjectedBoxTests(unittest.TestCase):
    def test_fixed_point_rounds_up(self):
        self.assertGreaterEqual(F(ceil_scaled(F(1, 3)), 10**18), F(1, 3))

    def test_log_and_sqrt_bounds(self):
        for q in (F(1, 3), F(1), F(2), F(5, 2)):
            lo, hi = log_interval(q)
            self.assertLessEqual(float(lo)-1e-15, math.log(float(q)))
            self.assertGreaterEqual(float(hi)+1e-15, math.log(float(q)))
            low, high = sqrt_down(q), sqrt_up(q)
            self.assertLessEqual(low**2, q)
            self.assertGreaterEqual(high**2, q)

    def test_small_response_against_independent_matrix_exponential(self):
        n, dx, nu, anchor, h = 8, F(1), F(8), F(1), 1
        result = response_upper(anchor, dx, n, nu, h)
        idx = result['anchor_index']
        coords = np.array([float(anchor + (i-idx)*dx) for i in range(n)])
        rate = float(nu/dx**2)
        potential = (21/25)*coords**2/(4*float(nu))
        generator = np.diag(-2*rate-potential) + np.diag(np.full(n-1, rate), 1) + np.diag(np.full(n-1, rate), -1)
        aug = np.zeros((n+1, n+1))
        aug[:n,:n] = generator
        aug[0,-1] = aug[n-1,-1] = rate
        actual = expm(aug)[:n,-1]
        for name, offset in [('center',0),('plus',1),('minus',-1)]:
            self.assertLessEqual(actual[idx+offset], float(result[name])+2e-14)
            self.assertLess(float(result[name])-actual[idx+offset], 1e-9)

    def test_parameterized_potential_against_matrix_exponential(self):
        n, dx, nu, anchor, h = 8, F(1), F(8), F(1), 1
        coeff=F(1,2)
        result=response_upper(anchor,dx,n,nu,h,potential_coeff=coeff)
        idx=result['anchor_index']
        coords=np.array([float(anchor+(i-idx)*dx) for i in range(n)])
        rate=float(nu/dx**2)
        potential=float(coeff)*coords**2/(4*float(nu))
        matrix=np.diag(-2*rate-potential)+np.diag(np.full(n-1,rate),1)+np.diag(np.full(n-1,rate),-1)
        aug=np.zeros((n+1,n+1));aug[:n,:n]=matrix;aug[0,-1]=aug[n-1,-1]=rate
        actual=expm(aug)[:n,-1]
        for name,offset in [('center',0),('plus',1),('minus',-1)]:
            self.assertLessEqual(actual[idx+offset],float(result[name])+2e-14)
            self.assertLess(float(result[name])-actual[idx+offset],1e-9)


if __name__ == '__main__':
    unittest.main()
