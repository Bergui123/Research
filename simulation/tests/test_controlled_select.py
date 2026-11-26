"""Exhaustive controlled SELECT check on a 512-dimensional block encoding."""
import unittest

import numpy as np

from model import X0, axis_grid
from select_circuit import Circuit, count, select, set_qubit
from step8_resources import W10, potential_scale


class ControlledSelectTest(unittest.TestCase):
    def test_controlled_select_and_inverse_on_every_basis_state(self):
        n, frac, m, dx, nu = 2, 2, 2, 0.5, 2.0
        npts = 2**n
        axes = [axis_grid(npts, npts*dx/2, X0[i], dx)[0] for i in range(2)]
        _, k_bits, q5 = potential_scale(nu, *axes, dx)
        offsets = [round(2**frac * axis[0]/dx) for axis in axes]
        control, j1, j2, t, k1, k0, r = np.indices((2,npts,npts,2,2,2,2**m)).reshape(7,-1)
        circuit = Circuit(control.size)
        load = lambda value, bits: [set_qubit(circuit, ((value >> i) & 1).astype(bool)) for i in range(bits)]
        regs = (load(j1,n), load(j2,n), *(set_qubit(circuit, x.astype(bool)) for x in (t,k1,k0)), load(r,m))
        ctrl = set_qubit(circuit, control.astype(bool))
        initial = [value.copy() for value in circuit.bits]
        base = circuit.live
        select(circuit,n,frac,offsets,k_bits,m,regs,W10,control=ctrl)
        read = lambda reg: sum(circuit.bits[q].astype(int) << i for i,q in enumerate(reg))
        kinetic = (control == 1) & (t == 0)
        step = np.where(k0, -1, 1) * kinetic
        self.assertTrue(np.array_equal(read(regs[0]), (j1+step*(1-k1))%npts))
        self.assertTrue(np.array_equal(read(regs[1]), (j2+step*k1)%npts))
        self.assertTrue(np.array_equal(circuit.bits[regs[4]], (k0 ^ kinetic).astype(bool)))
        q = np.floor(q5 * 2**(m-k_bits)).astype(int)[j1,j2]
        sign = kinetic | ((control == 1) & (t == 1) & ~(r < q))
        self.assertTrue(np.array_equal(circuit.sign, sign.astype(bool)))
        self.assertEqual(circuit.live, base)
        self.assertEqual(circuit.ands, count(n,frac,m,W10)[0] + 2)
        select(circuit,n,frac,offsets,k_bits,m,regs,W10,control=ctrl)
        self.assertEqual(circuit.live, base)
        self.assertFalse(circuit.sign.any())
        for before, after in zip(initial, circuit.bits):
            self.assertTrue(np.array_equal(before,after))


if __name__ == '__main__':
    unittest.main()
