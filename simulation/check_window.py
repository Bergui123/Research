"""Check the shifted-window condition c <= E_min on every grid used by a costed optimum, and
collect the wrap-around statistics quoted in the paper."""

import json

import numpy as np
from scipy.sparse.linalg import eigsh

from model import P_SHIFT, RESULTS, X0, axis_grid, grid_operator, save

best_a = json.loads((RESULTS / "step8.json").read_text())["best"].values()
best_b = json.loads((RESULTS / "step9_sampling.json").read_text())["best"].values()
ground = lambda hd: float(eigsh(hd.tocsc(), k=1, sigma=0, return_eigenvectors=False)[0])
out = []
for a in best_a:
	x1, _ = axis_grid(2 ** a["n"], 2 ** a["n"] * a["dx"] / 2, X0[0], a["dx"])
	x2, _ = axis_grid(2 ** a["n"], 2 ** a["n"] * a["dx"] / 2, X0[1], a["dx"])
	out.append(ground(grid_operator(x1, x2, a["dx"], a["nu"])[0]))
for b in best_b:
	x = (np.arange(2 ** b["n"]) - 2 ** b["n"] / 2) * b["dx"]
	out.append(ground(grid_operator(x, x, b["dx"], 1.0)[0]))
assert all(e > -P_SHIFT for e in out)
s7 = json.loads((RESULTS / "step7.json").read_text())
unsafe = [r["amp_err"] for r in s7["rows"] if not r["safe_Lp"] and r["n"] == 6 and r["p_a"] < 0 and r["kind"] in ("smooth_w1", "smooth_w2")]
res = {"c": -P_SHIFT, "E_min": out, "E_min_step7": min(s["E_min"] for s in s7["spectra"]),
	"unsafe_min": min(unsafe), "unsafe_max": max(unsafe), "weight40": max(s["weight_above"]["40"] for s in s7["spectra"])}
print(RESULTS.name, res)
save("window_check", res)
