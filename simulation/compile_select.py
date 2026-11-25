"""Exact simulation of the compiled SELECT at the costed sizes, and a comparison with Qualtran.

For every selected configuration the gate list of select_circuit.py is simulated on every grid point and on random
basis states: position, direction bit, sign and comparator integer must be right and every work qubit must return to
zero. Its T count, which step8_resources.py charges, is compared with Qualtran's block costs for the same layout.
"""

import json

import numpy as np
from qualtran import QInt, QUInt
from qualtran.bloqs.arithmetic import Add, AddK, LinearDepthHalfLessThan, Negate, Product, Square

from model import RESULTS, X0, axis_grid, save
from select_circuit import Circuit, count, select, set_qubit
from step8_resources import W10, potential_scale, tcount


def qualtran_select(n, frac, m):
	"""Block costs: two increments, compute and uncompute of the arithmetic, two half comparisons."""
	b = n + frac + 1
	kinetic = 2 * (tcount(AddK(QUInt(n), 1)) + 4)
	n_add = 4 + W10.bit_count()
	arith = (2 * tcount(AddK(QUInt(b), 1)) * (frac > 0) + 2 * tcount(Negate(QInt(b)))
		+ 2 * tcount(Square(b - 1)) + (W10 > 0) * (tcount(Product(b - 1, b - 1)) + tcount(Negate(QInt(2 * b))))
		+ n_add * tcount(Add(QUInt(2 * b + 3))))
	return kinetic + 2 * arith + 2 * tcount(LinearDepthHalfLessThan(QUInt(m)))


def check(cfg, readout, batch=200_000, seed=11):
	"""Simulate SELECT for one costed configuration on random basis states plus every grid point."""
	rng = np.random.default_rng(seed)
	n, m, dx = cfg["n"], cfg["m"], cfg["dx"]
	npts = 2**n
	if readout == "A":
		nu, frac = cfg["nu"], 2
		xs = [axis_grid(npts, cfg["L"], X0[a], dx)[0] for a in (0, 1)]
	else:
		nu, frac = 1.0, 0
		xs = [(np.arange(npts) - npts // 2) * dx] * 2
	_, k_bits, q5 = potential_scale(nu, xs[0], xs[1], dx)
	offsets = [int(round(2**frac * x[0] / dx)) for x in xs]
	grid = np.stack(np.meshgrid(np.arange(npts), np.arange(npts), indexing="ij"), -1).reshape(-1, 2)
	pos = np.concatenate([grid, rng.integers(0, npts, (batch, 2))])
	size = len(pos)
	anc = rng.integers(0, 2, (3, size)).astype(bool)
	anc[0, :len(grid)] = True  # every grid point through the potential branch
	rv = rng.integers(0, 2**m, size)
	c = Circuit(size)
	load = lambda v, k: [set_qubit(c, (v >> i & 1).astype(bool)) for i in range(k)]
	regs = (load(pos[:, 0], n), load(pos[:, 1], n), *[set_qubit(c, a) for a in anc], load(rv, m))
	base = c.live
	c.peak = c.live
	before = [b.copy() for b in c.bits]
	kinetic, q_circ = select(c, n, frac, offsets, k_bits, m, regs, W10)
	assert c.live == base
	q_ref = np.floor(q5 * 2.0 ** (m - k_bits)).astype(np.int64)[pos[:, 0], pos[:, 1]]
	t, k1, k0 = anc
	assert np.array_equal(q_circ[t], q_ref[t])  # t = 0 may have shifted the position first
	kin = ~t
	step = np.where(k0, -1, 1) * kin
	exp1 = (pos[:, 0] + step * ~k1) % npts
	exp2 = (pos[:, 1] + step * k1) % npts
	read = lambda reg: sum(c.bits[q].astype(np.int64) << i for i, q in enumerate(reg))
	assert np.array_equal(read(regs[0]), exp1) and np.array_equal(read(regs[1]), exp2)
	assert np.array_equal(c.bits[regs[4]], k0 ^ kin) and np.array_equal(read(regs[5]), rv)
	for i in (2, 3):
		assert np.array_equal(c.bits[regs[i]], before[regs[i]])
	assert np.array_equal(c.sign, kin | (t & ~(rv < q_ref)))
	ands, kin_ands, work = count(n, frac, m, W10)
	assert (c.ands, kinetic, c.peak - base) == (ands, kin_ands, work)
	assert 4 * ands == cfg["T_walk_parts"]["select_kinetic"] + cfg["T_walk_parts"]["potential_oracle"]
	assert work == cfg["qubits_parts"]["arith_work"]
	return {"readout": readout, "target": cfg["target"], "n": n, "m": m, "basis_states": size, "T_select": 4 * ands,
		"T_kinetic": 4 * kin_ands, "T_select_qualtran": qualtran_select(n, frac, m), "work": work}


if __name__ == "__main__":
	# The adders on every input of a small register.
	a, t = np.meshgrid(np.arange(32), np.arange(64), indexing="ij")
	for inv in (False, True):
		c = Circuit(a.size)
		ra = [set_qubit(c, (a.ravel() >> i & 1).astype(bool)) for i in range(5)]
		rt = [set_qubit(c, (t.ravel() >> i & 1).astype(bool)) for i in range(6)]
		c.add(ra, rt, inv)
		out = sum(c.bits[q].astype(int) << i for i, q in enumerate(rt))
		assert np.array_equal(out, (t.ravel() + (-1 if inv else 1) * a.ravel()) % 64) and c.live == 11 and c.ands == 5
	for k in range(2, 12):  # the squarer on every input
		v = np.arange(2**k)
		c = Circuit(v.size)
		rv, out = [set_qubit(c, (v >> i & 1).astype(bool)) for i in range(k)], c.alloc(2 * k)
		c.square(rv, out)
		assert np.array_equal(sum(c.bits[q].astype(int) << i for i, q in enumerate(out)), v * v)
		c.square(rv, out, inverse=True)
		c.release(out)
	rows = [check(b, "A") for b in json.loads((RESULTS / "step8.json").read_text())["best"].values()]
	rows += [check(b, "B") for b in json.loads((RESULTS / "step9_sampling.json").read_text())["best"].values()]
	for r in rows:
		print(RESULTS.name, r)
	save("compile_select", rows)
