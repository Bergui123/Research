"""Derive a deterministic Hamilton-Jacobi-Bellman equation symbolically.

Supply a running-cost integrand, state variables, controls, and state dynamics.
The problem represented is

	V(x, t) = min_u integral_t^T L(x(s), u(s), s) ds + terminal_cost(x(T))
	x_dot = f(x, u, t).

The control is eliminated using its first-order optimality condition. This
returns the minimizing stationary control; convexity of the Hamiltonian in
the controls is needed for it to be a global minimum.
"""

from dataclasses import dataclass
from typing import Sequence

import sympy as sp


@dataclass(frozen=True)
class HJBResult:
	"""Symbolic value function, optimal control, and HJB equation."""

	value_function: sp.Expr
	hamiltonian: sp.Expr
	minimized_hamiltonian: sp.Expr
	optimal_control: dict[sp.Symbol, sp.Expr]
	pde: sp.Equality
	terminal_condition: sp.Equality


def derive_hjb(
	running_cost: sp.Expr,
	states: Sequence[sp.Symbol],
	controls: Sequence[sp.Symbol],
	dynamics: Sequence[sp.Expr],
	time: sp.Symbol,
	terminal_time: sp.Expr,
	terminal_cost: sp.Expr = sp.Integer(0),
) -> HJBResult:
	"""Build the deterministic HJB PDE from an integral's ingredients.

	``running_cost`` is the integrand L, and ``dynamics`` lists x_dot = f
	in the same order as ``states``. ``terminal_cost`` is the cost at T,
	expressed using the state symbols. Hard terminal constraints need an
	appropriate penalty or a restricted state domain.
	"""
	if not states:
		raise ValueError("Provide at least one state variable.")
	if len(dynamics) != len(states):
		raise ValueError("Provide one dynamics expression for each state.")
	if len(set(states)) != len(states) or len(set(controls)) != len(controls):
		raise ValueError("State and control variables must be unique.")
	if set(states) & set(controls):
		raise ValueError("State and control variables must be distinct.")

	value = sp.Function("V")(time, *states)
	value_gradient = [sp.diff(value, state) for state in states]
	hamiltonian = sp.sympify(running_cost) + sum(
		gradient * drift for gradient, drift in zip(value_gradient, dynamics)
	)

	if controls:
		stationarity = [sp.diff(hamiltonian, control) for control in controls]
		candidates = sp.solve(stationarity, controls, dict=True)
		complete_candidates = [
			candidate for candidate in candidates
			if all(control in candidate for control in controls)
		]
		if len(complete_candidates) != 1:
			raise ValueError(
				"Could not determine one unique stationary control. "
				"Check that the Hamiltonian has a unique minimizer, or "
				"supply a model with a solvable first-order condition."
			)
		optimal_control = complete_candidates[0]
	else:
		optimal_control = {}

	minimized_hamiltonian = sp.simplify(hamiltonian.subs(optimal_control))
	pde = sp.Eq(
		sp.diff(value, time) + minimized_hamiltonian,
		0,
	)
	terminal_condition = sp.Eq(
		value.subs(time, terminal_time),
		sp.sympify(terminal_cost),
	)

	return HJBResult(
		value_function=value,
		hamiltonian=sp.simplify(hamiltonian),
		minimized_hamiltonian=minimized_hamiltonian,
		optimal_control=optimal_control,
		pde=pde,
		terminal_condition=terminal_condition,
	)


def almgren_chriss_example() -> HJBResult:
	"""Derive the one-asset PDE for L = eta*v^2 + lambda*sigma^2*x^2."""
	x, selling_rate, time, horizon = sp.symbols(
		"x v t T", real=True
	)
	eta, risk_aversion, variance, penalty = sp.symbols(
		"eta lambda sigma2 Phi", positive=True, real=True
	)
	running_cost = eta * selling_rate**2 + risk_aversion * variance * x**2
	return derive_hjb(
		running_cost=running_cost,
		states=[x],
		controls=[selling_rate],
		dynamics=[-selling_rate],
		time=time,
		terminal_time=horizon,
		terminal_cost=sp.Rational(1, 2) * penalty * x**2,
	)


if __name__ == "__main__":
	result = almgren_chriss_example()
	print("Optimal selling rate:")
	sp.pprint(result.optimal_control)
	print("\nHJB PDE:")
	sp.pprint(result.pde)
	print("\nTerminal condition:")
	sp.pprint(result.terminal_condition)
