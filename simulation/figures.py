"""Report figures from results/*.json (static PNG for LaTeX; values are also tabulated in the report)."""

import json

import matplotlib
import matplotlib.ticker
import numpy as np

from model import FIGURES, RESULTS

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Validated categorical order (dataviz reference palette); markers give the non-colour encoding.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
MARKERS = ["o", "s", "^", "D"]
plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#52514e", "axes.labelcolor": "#0b0b0b",
	"xtick.color": "#52514e", "ytick.color": "#52514e", "axes.spines.top": False,
	"axes.spines.right": False, "lines.linewidth": 1.5, "lines.markersize": 5,
	"figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb"})


def load(name):
	return json.loads((RESULTS / f"{name}.json").read_text())


def style(ax):
	for axis in (ax.xaxis, ax.yaxis):
		axis.set_minor_formatter(matplotlib.ticker.NullFormatter())  # minor log labels collide
	ax.grid(True, which="major", color="#d8d7d2", linewidth=0.6)
	ax.set_axisbelow(True)


def step4():
	rows = [r for r in load("step4")["grid"] if r["bc"] == "periodic"]
	fig, ax = plt.subplots(figsize=(4.6, 3.2))
	for i, L in enumerate((4.0, 5.0, 6.0, 8.0)):
		pts = sorted((r["dx"], r["rel_err_v"]) for r in rows if r["L"] == L)
		ax.loglog(*zip(*pts), marker=MARKERS[i], color=SERIES[i], label=f"$L={L:g}$")
	dx = np.array([0.03, 0.5])
	ax.loglog(dx, 0.4 * 0.222 * (dx / 0.5) ** 2, ":", color="#52514e", label=r"slope 2")
	ax.set_xlabel(r"grid spacing $\Delta x$")
	ax.set_ylabel(r"relative feedback error")
	ax.set_title(r"Step 4: $\nu=0.2$, periodic grid", fontsize=9)
	style(ax)
	ax.legend(frameon=False)
	fig.tight_layout()
	fig.savefig(FIGURES / "step4_convergence.png", dpi=200)


def step5():
	grid = load("step5")["grid"]
	rows = load("grid_error")
	fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
	pts = sorted((r["nu"], r["measured"], r["leading"], r["cancellation_free"]) for r in rows
		if r["n"] == 6 and r["h"] == 1 and r["nu"] >= 0.5)
	nu, meas, lead, cf = zip(*pts)
	axes[0].loglog(nu, cf, marker=MARKERS[1], color=SERIES[1], label="cancellation-free estimate")
	axes[0].loglog(nu, lead, color="#52514e", linestyle=":", label="leading order (Thm. 3)")
	axes[0].loglog(nu, meas, marker=MARKERS[0], color=SERIES[0], linestyle="none", label="measured")
	axes[0].set_ylabel("relative feedback error")
	axes[0].set_title("$N=64$, $h=1$", fontsize=9)
	axes[0].legend(frameon=False, fontsize=8)
	for i, n in enumerate((5, 6, 7, 8)):
		p = sorted((r["nu"], min(r["p_neighbours"])) for r in grid if r["n"] == n)
		axes[1].loglog(*zip(*p), marker=MARKERS[i], color=SERIES[i], label=f"$N={2**n}$")
	axes[1].set_ylabel(r"smallest neighbour probability $p_j$")
	axes[1].legend(frameon=False, fontsize=8, loc="lower right")
	for ax in axes:
		ax.set_xlabel(r"inventory diffusion $\nu$")
		style(ax)
	fig.tight_layout()
	fig.savefig(FIGURES / "step5_nu_tradeoff.png", dpi=200)


def step7():
	rows = load("step7")["rows"]
	fig, ax = plt.subplots(figsize=(4.6, 3.2))
	labels = {"kink": r"$e^{-|p|}$", "smooth_w0.5": r"$g_{0.5}$", "smooth_w1": r"$g_1$", "smooth_w2": r"$g_2$"}
	for i, kind in enumerate(labels):
		pts = sorted((r["eta_max"], r["amp_err"]) for r in rows if r["kind"] == kind and r["nu"] == 2.0
			and r["n"] == 6 and r["p_a"] < 0 and r["safe_Lp"])
		ax.loglog(*zip(*pts), marker=MARKERS[i], color=SERIES[i], label=labels[kind])
	ax.set_xlabel(r"frequency cut-off $\eta_{\max}$")
	ax.set_ylabel("largest relative amplitude error")
	style(ax)
	ax.legend(frameon=False)
	fig.tight_layout()
	fig.savefig(FIGURES / "step7_profiles.png", dpi=200)


def step8():
	a = load("step8")
	b = load("step9_sampling")["best"]
	fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
	eps = sorted(float(k) for k in a["best"])
	get = lambda d, e: next(v for k, v in d.items() if float(k) == e)
	for i, (label, y) in enumerate((("A: amplitudes, $r=11$", [get(a["best"], e)["T_total"] for e in eps]),
			("A: amplitudes, one set", [get(a["best"], e)["T_total_no_boost"] for e in eps]),
			("B: sampling", [get(b, e)["T_total"] for e in eps]))):
		axes[0].loglog(eps, y, marker=MARKERS[i], color=SERIES[i], label=label)
	axes[0].set_xlabel(r"target error $\varepsilon$")
	axes[0].set_ylabel("total T gates")
	for i, e in enumerate(eps[::-1]):
		pts = sorted((r["nu"], r["T_total"]) for r in a["per_nu"] if r["target"] == e)
		axes[1].loglog(*zip(*pts), marker=MARKERS[i], color=SERIES[i], label=fr"$\varepsilon={e:g}$")
	axes[1].set_xlabel(r"inventory diffusion $\nu$")
	axes[1].set_ylabel("total T gates, readout A")
	for ax in axes:
		style(ax)
	axes[0].legend(frameon=False, fontsize=8)
	axes[1].legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(0.42, 0.62))
	fig.tight_layout()
	fig.savefig(FIGURES / "tcount.png", dpi=200)


def step10():
	a = load("step10")["emulation"]
	b = load("step9_sampling")["best"]
	fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8), sharex=True)
	for i, r in enumerate(a):
		axes[0].hist(np.array(r["err_samples"]) / r["target"], bins=40, range=(0, 1.2), histtype="step",
			linewidth=1.5, color=SERIES[i], label=fr"$\varepsilon={r['target']:g}$")
		v = next(x for k, x in b.items() if float(k) == r["target"])["validation"]
		axes[1].hist(np.array(v["errors_at_X"]) / r["target"], bins=40, range=(0, 1.2), histtype="step",
			linewidth=1.5, color=SERIES[i], label=fr"$\varepsilon={r['target']:g}$")
	for ax, title in zip(axes, ("A: amplitude estimation, median of 11", "B: sampling")):
		ax.axvline(1.0, color="#52514e", linestyle=":", linewidth=1)
		ax.set_xlabel(r"feedback error at $X$ / target")
		ax.set_title(title, fontsize=9)
		style(ax)
	axes[0].legend(frameon=False, fontsize=8)  # same series in both panels
	axes[0].set_ylabel("trials (of 2000)")
	fig.tight_layout()
	fig.savefig(FIGURES / "step10_emulation.png", dpi=200)


def scaling():
	rows = load("scaling_m")["rows"]
	fig, ax = plt.subplots(figsize=(4.6, 3.2))
	m = [r["M"] for r in rows]
	for i, (key, label) in enumerate((("runs_A", "A: amplitudes"), ("runs_B_plain", "B: sampling, $p\\geq0$"),
			("runs_B", "B: sampling, shifted window"))):
		ax.semilogy(m, [r[key] for r in rows], marker=MARKERS[i], color=SERIES[i], label=label)
	ax.set_xlabel("number of assets $M$")
	ax.set_ylabel(r"coherent runs, $\varepsilon=1\%$")
	style(ax)
	ax.legend(frameon=False, fontsize=8)
	fig.tight_layout()
	fig.savefig(FIGURES / "scaling_m.png", dpi=200)


def nonquadratic():
	a = load("nonquadratic")
	fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9))
	nu = [r["nu"] for r in a["bias"]]
	for i, (key, label) in enumerate((("kappa0.5", r"$\kappa=0.5$"), ("kappa0.1", r"$\kappa=0.1$"), ("kappa0.02", r"$\kappa=0.02$"))):
		axes[0].loglog(nu, [r[key] for r in a["bias"]], marker=MARKERS[i], color=SERIES[i], label=label, markevery=3)
	axes[0].set_xlabel(r"inventory diffusion $\nu$")
	axes[0].set_ylabel(r"bias of the viscous rate at $X$")
	w = a["W0"]
	nus = np.geomspace(0.05, 8, 60)
	axes[1].loglog(nus, np.exp(-w / (4 * nus)), color="#52514e", linestyle=":", label=r"barrier $e^{-W(X)/4\nu}$")
	lev = [r for r in a["direct_amplitudes"] if 0.05 <= r["nu"] <= 8]
	axes[1].loglog([r["nu"] for r in lev], [r["direct"] for r in lev], marker=MARKERS[0], color=SERIES[0], label="direct", markevery=2)
	axes[1].loglog([r["nu"] for r in lev], [r["recentred"] for r in lev], marker=MARKERS[1], color=SERIES[1], label="recentred", markevery=2)
	axes[1].set_ylim(1e-26, 10)
	axes[1].set_xlabel(r"inventory diffusion $\nu$")
	axes[1].set_ylabel(r"smallest readout amplitude")
	for ax in axes:
		style(ax)
		ax.legend(frameon=False, fontsize=8, loc="lower right")
	fig.tight_layout()
	fig.savefig(FIGURES / "nonquadratic.png", dpi=200)


if __name__ == "__main__":
	FIGURES.mkdir(exist_ok=True)
	step4(), step5(), step7(), step8(), step10(), scaling(), nonquadratic()
