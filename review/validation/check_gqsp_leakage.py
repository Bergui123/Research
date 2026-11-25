"""Stress the two/three-segment leak-and-return estimate with exact 2x2 matrices.

This checks the block-matrix estimate in isolation. It does not synthesize Berry's
completion polynomials or validate the Jacobi--Anger truncation bound.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.special import jv


SEED = 20261124
TRIALS = 100_000
ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "review" / "validation" / "gqsp-leakage-results.json"


def tail(k: int, tau: float) -> float:
	return float(2 * np.sum(np.abs(jv(np.arange(k + 1, 300), tau))))


def run() -> dict:
	rng = np.random.default_rng(SEED)
	rows = []
	for k in (2, 4, 6):
		for tau in (0.1, 0.5, 1.0):
			eps = tail(k, tau)
			ell = 2 * np.sqrt(eps)
			if ell >= 1:
				continue
			c = np.sqrt(1 - ell**2)
			for nseg in (2, 3):
				phase = rng.uniform(-np.pi, np.pi, size=(TRIALS, nseg))
				mats = np.empty((TRIALS, nseg, 2, 2), dtype=np.complex128)
				mats[:, :, 0, 0] = c
				mats[:, :, 1, 1] = c
				mats[:, :, 0, 1] = ell * np.exp(1j * phase)
				mats[:, :, 1, 0] = -ell * np.exp(-1j * phase)
				prod = mats[:, 0]
				for j in range(1, nseg):
					prod = mats[:, j] @ prod
				# Diagonal-only paths never leave |+>; all remaining paths leak
				# out of |+> and return before the final segment.
				leak_return = np.abs(prod[:, 0, 0] - c**nseg)
				claimed = 2 * nseg**2 * eps
				pair_count_bound = (nseg * (nseg - 1) / 2) * ell**2
				rows.append({
					"K": k,
					"tau": tau,
					"segments": nseg,
					"tail": eps,
					"single_segment_leak_cap": ell,
					"claimed_bound_2_n2_tail": claimed,
					"analytic_path_bound_choose_n2_ell2": pair_count_bound,
					"sampled_max_leak_return": float(leak_return.max()),
					"sampled_max_over_claimed": float(leak_return.max() / claimed),
					"sampled_violations": int(np.sum(leak_return > claimed * (1 + 1e-12))),
				})
	return {
		"method": "Exact 2x2 unitary products with per-segment off-diagonal magnitude 2*sqrt(tail); random phases.",
		"seed": SEED,
		"trials_per_case": TRIALS,
		"cases": rows,
		"analytic_check": "For 2 or 3 segments, a path that leaves and returns to |+> has exactly two off-diagonal transitions. There are at most choose(n,2) such paths; each contributes at most ell^2=4*tail. This is <=2*n^2*tail.",
		"scope": "This verifies the claimed two/three-segment leakage estimate for the abstract block model. It does not verify the bound on Berry completion leakage, its normalization, or multi-segment direct approximation errors.",
	}


if __name__ == "__main__":
	data = run()
	OUT.write_text(json.dumps(data, indent=2) + "\n")
	print(json.dumps({"cases": len(data["cases"]), "max_ratio": max(r["sampled_max_over_claimed"] for r in data["cases"]),
		"violations": sum(r["sampled_violations"] for r in data["cases"]), "output": str(OUT)}, indent=2))
