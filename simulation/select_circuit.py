"""SELECT of the two-asset Hamiltonian as a list of X, CNOT, AND, AND^dagger, Z and CZ gates.

SELECT is a permutation of basis states with signs, so the gate list can be simulated exactly on computational
basis states (compile_select.py). AND costs 4 T; AND^dagger is measurement based and free (Gidney 2018). The T count
and the workspace charged in step8_resources.py are those of this circuit.
"""

import functools

import numpy as np


class Circuit:
	"""Basis-state simulator on a batch; a qubit is an index into `bits`."""

	def __init__(self, batch):
		self.bits, self.free, self.live, self.peak, self.ands = [], [], 0, 0, 0
		self.sign = np.zeros(batch, bool)
		self.batch = batch

	def alloc(self, k=1):
		out = []
		for _ in range(k):
			if self.free:
				out.append(self.free.pop())
			else:
				self.bits.append(np.zeros(self.batch, bool))
				out.append(len(self.bits) - 1)
		self.live += k
		self.peak = max(self.peak, self.live)
		return out

	def release(self, qs):
		for q in qs:
			assert not self.bits[q].any(), "work qubit not returned to zero"
			self.free.append(q)
		self.live -= len(qs)

	def x(self, q):
		self.bits[q] ^= True

	def cnot(self, c, t):
		self.bits[t] ^= self.bits[c]

	def and_(self, a, b):
		(t,) = self.alloc()
		self.bits[t] ^= self.bits[a] & self.bits[b]
		self.ands += 1
		return t

	def unand(self, a, b, t):
		self.bits[t] ^= self.bits[a] & self.bits[b]
		self.release([t])

	def z(self, q):
		self.sign ^= self.bits[q]

	def cz(self, a, b):
		self.sign ^= self.bits[a] & self.bits[b]

	def add(self, a, t, inverse=False):
		"""t += a (mod 2^len(t)), len(a) <= len(t); Gidney adder, len(t) - 1 ANDs at most."""
		if inverse:
			for q in t:
				self.x(q)
		carries = []  # (left operand, right operand, AND qubit)
		c = None
		for i in range(len(t) - 1):
			ai = a[i] if i < len(a) else None
			if ai is None and c is None:
				break
			if c is None:
				carries.append((ai, t[i], self.and_(ai, t[i])))
			elif ai is None:
				carries.append((c, t[i], self.and_(c, t[i])))
			else:  # maj(a, t, c) = ((a ^ c)(t ^ c)) ^ c
				self.cnot(c, ai), self.cnot(c, t[i])
				new = self.and_(ai, t[i])
				self.cnot(c, new)
				carries.append((ai, t[i], new))
			c = carries[-1][2]
		for i in range(len(t) - 1, -1, -1):
			ai = a[i] if i < len(a) else None
			cin = carries[i - 1][2] if 0 < i <= len(carries) else None
			if i < len(carries):
				left, right, cout = carries[i]
				if ai is not None and cin is not None:
					self.cnot(cin, cout)
					self.unand(left, right, cout)
					self.cnot(cin, ai), self.cnot(cin, t[i])
				else:
					self.unand(left, right, cout)
			if cin is not None:
				self.cnot(cin, t[i])
			if ai is not None:
				self.cnot(ai, t[i])
		if inverse:
			for q in t:
				self.x(q)

	def add_constant(self, k, t):
		a = self.alloc(len(t))
		for i, q in enumerate(a):
			if k >> i & 1:
				self.x(q)
		self.add(a, t)
		for i, q in enumerate(a):
			if k >> i & 1:
				self.x(q)
		self.release(a)

	def controlled_increment(self, c, t):
		self.add([c], t)

	def row(self, ctrl, a, t, inverse):
		"""t +-= ctrl * a with a temporary row of ANDs."""
		tmp = [self.and_(ctrl, q) for q in a]
		self.add(tmp, t, inverse)
		for q, w in zip(a, tmp):
			self.unand(ctrl, q, w)

	def square(self, a, out, inverse=False):
		"""out +-= a^2 = sum_i a_i 4^i + sum_{i<j} a_i a_j 2^{i+j+1}; the carries of row i end three bits
		above its last term (checked on every input in the self-test)."""
		n = len(a)
		if not inverse:
			for i in range(n):  # out = 0 before: copy the diagonal, no carries
				self.cnot(a[i], out[2 * i])
		for i in range(n - 1) if not inverse else range(n - 2, -1, -1):
			self.row(a[i], a[i + 1:], out[2 * i + 2:i + n + 4], inverse)
		if inverse:
			for i in range(n):
				self.cnot(a[i], out[2 * i])

	def multiply(self, a, b, out, inverse=False):
		n = len(a)
		for i in range(n) if not inverse else range(n - 1, -1, -1):
			if i == 0:  # out = 0: the first row is written directly
				for q, t in zip(b, out):
					if not inverse:
						self.ands += 1
					self.bits[t] ^= self.bits[a[0]] & self.bits[q]
			else:
				self.row(a[i], b, out[i:i + len(b) + 1], inverse)


def select(c, n, frac, offsets, k_bits, m, regs, w10, control=None):
	"""SELECT of the two-asset Hamiltonian. `offsets` are the classical constants of the coordinates
	2^frac x_a / dx = 2^frac j_a + offset_a; q is the top m bits of the integer 4^frac 5 x^T C x / dx^2.
	If `control` is supplied, apply SELECT only when that qubit is one."""
	j1, j2, t, k1, k0, r = regs
	# Kinetic part: -S_a^{+-1} and a flip of the direction bit, controlled on t = 0.
	c.x(t)
	on = t if control is None else c.and_(t, control)
	ax1 = c.and_(on, k1)
	(ax0,) = c.alloc()
	c.cnot(on, ax0), c.cnot(ax1, ax0)
	for ctrl, reg in ((ax0, j1), (ax1, j2)):
		for q in reg:
			c.cnot(k0, q)  # k0 = 1: conjugated increment = decrement
		c.controlled_increment(ctrl, reg)
		for q in reg:
			c.cnot(k0, q)
	c.z(on)
	c.cnot(on, k0)
	c.cnot(on, ax0), c.cnot(ax1, ax0)
	c.release([ax0])
	c.unand(on, k1, ax1)
	if control is not None:
		c.unand(t, control, on)
	c.x(t)
	kinetic = c.ands

	# Potential part: coordinates, squares and product, accumulator, comparison, phase, uncompute.
	b = n + frac + 1
	zs, signs = [], []
	for reg, off in ((j1, offsets[0]), (j2, offsets[1])):
		z = c.alloc(b)
		for i, q in enumerate(reg):
			c.cnot(q, z[frac + i])
		low = off % 2**frac
		for i in range(frac):
			if low >> i & 1:
				c.x(z[i])
		c.add_constant(((off - low) >> frac) % 2 ** (b - frac), z[frac:])
		for q in z[:-1]:  # |z|: conditional complement, then add the sign bit
			c.cnot(z[-1], q)
		c.controlled_increment(z[-1], z[:-1])
		zs.append(z[:-1]), signs.append(z[-1])
	w = b - 1
	s1, s2 = c.alloc(2 * w), c.alloc(2 * w)
	p = c.alloc(2 * w) if w10 else []
	c.square(zs[0], s1), c.square(zs[1], s2)
	if w10:
		c.multiply(zs[0], zs[1], p)
		c.cnot(signs[0], signs[1])  # sign of the product
	acc = c.alloc(2 * b + 3)
	terms = [(s1, 0), (s1, 2), (s2, 0), (s2, 2)] + [(p, i) for i in range(w10.bit_length()) if w10 >> i & 1]

	def accumulate(inverse):
		for reg, sh in terms if not inverse else terms[::-1]:
			if reg is p:
				for q in acc:
					c.cnot(signs[1], q)
			c.add(reg, acc[sh:], inverse)
			if reg is p:
				for q in acc:
					c.cnot(signs[1], q)
	accumulate(False)
	top = k_bits + 2 * frac
	q_reg = acc[top - m:top]
	# [r < q]: carry out of q + ~r, i.e. q - r - 1 >= 0.
	for q in r:
		c.x(q)
	chain, carry = [], None
	for rq, qq in zip(r, q_reg):
		if carry is None:
			new = c.and_(rq, qq)
		else:
			c.cnot(carry, rq), c.cnot(carry, qq)
			new = c.and_(rq, qq)
			c.cnot(carry, new)
		chain.append((rq, qq, carry, new))
		carry = new
	potential_on = t if control is None else c.and_(t, control)
	c.z(potential_on)
	c.cz(potential_on, carry)  # -(-1)^{[r < q]} on the enabled potential branch
	if control is not None:
		c.unand(t, control, potential_on)
	for rq, qq, cin, new in chain[::-1]:
		if cin is not None:
			c.cnot(cin, new)
		c.unand(rq, qq, new)
		if cin is not None:
			c.cnot(cin, rq), c.cnot(cin, qq)
	for q in r:
		c.x(q)
	value = sum(c.bits[q].astype(np.int64) << i for i, q in enumerate(q_reg))
	accumulate(True)
	c.release(acc)
	if w10:
		c.cnot(signs[0], signs[1])
		c.multiply(zs[0], zs[1], p, inverse=True)
		c.release(p)
	c.square(zs[1], s2, inverse=True), c.square(zs[0], s1, inverse=True)
	c.release(s1 + s2)
	for (reg, off), mag, sg in zip(((j2, offsets[1]), (j1, offsets[0])), zs[::-1], signs[::-1]):
		for q in mag:
			c.x(q)
		c.controlled_increment(sg, mag)  # x -> ~(~x + s) subtracts the sign bit
		for q in mag:
			c.x(q)
		for q in mag:
			c.cnot(sg, q)
		z = mag + [sg]
		low = off % 2**frac
		c.add_constant(-((off - low) >> frac) % 2 ** (b - frac), z[frac:])
		for i in range(frac):
			if low >> i & 1:
				c.x(z[i])
		for i, q in enumerate(reg):
			c.cnot(q, z[frac + i])
		c.release(z)
	return kinetic, value


def set_qubit(c, values):
	(q,) = c.alloc()
	c.bits[q] = values.copy()
	return q


@functools.lru_cache(maxsize=None)
def count(n, frac, m, w10):
	"""(ANDs, ANDs of the kinetic part, peak work qubits); independent of the input and of the constants."""
	c = Circuit(1)
	reg = lambda k: [set_qubit(c, np.zeros(1, bool)) for _ in range(k)]
	regs = (reg(n), reg(n), *reg(3), reg(m))
	base = c.peak = c.live
	b = n + frac + 1
	kinetic, _ = select(c, n, frac, (0, 0), max(m - 2 * frac, 2 * b + 3 - 2 * frac), m, regs, w10)
	return c.ands, kinetic, c.peak - base
