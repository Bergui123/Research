"""Write paper/generated.tex: every number and table body quoted in the paper, from results/*.json."""

import json
import math

from model import ROOT

R = ROOT / "results"
load = lambda n, d=R: json.loads((d / f"{n}.json").read_text())
out = []


def sci(x, d=1):
	if x == 0:
		return "0"
	e = math.floor(math.log10(abs(x)))
	mant = round(x / 10**e, d)
	if mant >= 10:
		mant, e = mant / 10, e + 1
	return f"{mant:.{d}f}\\times10^{{{e}}}"


def macro(name, value):
	out.append(f"\\newcommand{{\\{name}}}{{{value}}}")


def rows(name, lines):
	out.append(f"\\newcommand{{\\{name}}}{{%\n" + "\n".join(lines) + "}")


KIND = {"kink": "e^{-|p|}", "smooth_w0.5": "g_{1/2}", "smooth_w1": "g_1", "smooth_w2": "g_2"}
at = lambda d, e: next(v for k, v in d.items() if float(k) == e)
TAG = {0.05: "Five", 0.01: "One", 0.001: "Tenth"}
s8, s9, s10 = load("step8")["best"], load("step9_sampling"), load("step10")["emulation"]
EPS = (0.05, 0.01, 0.001)
A, B = [at(s8, e) for e in EPS], [at(s9["best"], e) for e in EPS]
E = [next(r for r in s10 if r["target"] == e) for e in EPS]
for e, a, b in zip(EPS, A, B):
	macro(f"TA{TAG[e]}", sci(a["T_total"]))
	macro(f"TAsingle{TAG[e]}", sci(a["T_total_no_boost"]))
	macro(f"QA{TAG[e]}", a["logical_qubits"])
	macro(f"TB{TAG[e]}", sci(b["T_total"]))
	macro(f"QB{TAG[e]}", b["logical_qubits"])
macro("RatioAB", f"{A[1]['T_total'] / B[1]['T_total']:.0f}")
macro("RatioABSingle", f"{A[1]['T_total_no_boost'] / B[1]['T_total']:.1f}")
macro("StepsA", A[1]["walk_steps_per_run"]); macro("StepsB", B[1]["walk_steps_per_run"])
macro("ReadoutRuns", sci(A[1]["calls_per_estimate_set"]))
macro("RunsBOne", sci(B[1]["runs"]))
macro("OracleShareA", f"{100 * A[1]['T_walk_parts']['potential_oracle'] / A[1]['T_walk_step']:.0f}")
macro("PrepShareA", f"{100 * A[1]['T_prep'] / A[1]['T_run']:.1f}")
macro("SerialDaysB", f"{B[1]['T_total'] * 1e-6 / 86400:.0f}")
macro("SerialYearsA", f"{A[1]['T_total'] * 1e-6 / 86400 / 365.25:.1f}")
Or = [r["T_walk_parts"]["potential_oracle"] / r["T_walk_step"] for r in A + B]
macro("OracleMin", f"{100 * min(Or):.0f}"); macro("OracleMax", f"{100 * max(Or):.0f}")
macro("KappaB", f"{s9['kappa']:.2f}")


def line(label, vals):
	return label + " & " + " & ".join(str(v) for v in vals) + r"\\"


rows("TableA", [
	line(r"$\nu$", [f"{a['nu']:g}" for a in A]), line(r"grid $N\times N$; $\Delta x$", [f"${2**a['n']}^2$; {a['dx']:.3f}" for a in A]),
	line(r"readout span $h$", [a["h"] for a in A]),
	line(r"$p$-profile; $\eta_{\max}$; $n_p$", [f"${KIND[a['kind']]}$; {a['eta_max']:.0f}; {a['n_p']}" for a in A]),
	line(r"grid error charged (measured)", [f"${sci(a['eps_grid'])}$ (${sci(a['eps_grid_measured'])}$)" for a in A]),
	line(r"LCU normalisation $\alpha$", [f"{a['alpha_walk']:.0f}" for a in A]), line(r"$P_{\rm succ}$", [f"{a['P_succ']:.2f}" for a in A]),
	line("walk steps per run (standard QSP)", [f"{a['walk_steps_per_run']} ({a['walk_steps_std_qsp']})" for a in A]),
	line("T per walk step (potential oracle)", [f"{a['T_walk_step']:.0f} ({a['T_walk_parts']['potential_oracle']:.0f})" for a in A]),
	line("T per coherent run", [f"${sci(a['T_run'])}$" for a in A]),
	line(r"$\log_2M_j$, four amplitudes", [", ".join(map(str, a["AE_log2_iterations"])) for a in A]),
	line("coherent runs per estimate set", [f"${sci(a['calls_per_estimate_set'])}$" for a in A]),
	line("total T, one estimate set", [f"${sci(a['T_total_no_boost'])}$" for a in A]),
	line(r"\textbf{total T, median of 11}", [f"$\\mathbf{{{sci(a['T_total'])}}}$" for a in A]),
	line(r"\textbf{logical qubits}", [f"\\textbf{{{a['logical_qubits']}}}" for a in A])])
rows("TableB", [
	line(r"grid $N\times N$; half-width $L/\sqrt\nu$", [f"${2**b['n']}^2$; {b['L']:g}" for b in B]),
	line(r"$p$-profile; $\eta_{\max}$; $n_p$", [f"${KIND[b['kind']]}$; {b['eta_max']:.0f}; {b['n_p']}" for b in B]),
	line(r"systematic error $r$", [f"${sci(b['rho_systematic'])}$" for b in B]),
	line(r"LCU normalisation $\alpha$", [f"{b['alpha_walk']:.0f}" for b in B]), line(r"$P_{\rm succ}$", [f"{b['P_succ']:.2f}" for b in B]),
	line("walk steps per run", [b["walk_steps_per_run"] for b in B]),
	line("T per walk step (potential oracle)", [f"{b['T_walk_step']:.0f} ({b['T_walk_parts']['potential_oracle']:.0f})" for b in B]),
	line("T per coherent run", [f"${sci(b['T_run'])}$" for b in B]),
	line("accepted samples", [f"${sci(b['accepted_samples'])}$" for b in B]),
	line("coherent runs", [f"${sci(b['runs'])}$" for b in B]),
	line(r"\textbf{total T, 95\% target}", [f"$\\mathbf{{{sci(b['T_total'])}}}$" for b in B]),
	line(r"\textbf{logical qubits}", [f"\\textbf{{{b['logical_qubits']}}}" for b in B])])
q = lambda d: [d["qubits_parts"][k] for k in ("x", "p", "prep_lcu", "gqsp", "arith_work", "phase_grad")]
rows("TableQ", [line(f"A, {100 * e:g}\\%", q(a) + [1, a["logical_qubits"]]) for e, a in zip(EPS, A)]
	+ [line(f"B, {100 * e:g}\\%", q(b) + [0, b["logical_qubits"]]) for e, b in zip(EPS, B)])
rows("TableE", [line(f"{100 * e:g}\\%", [f"${sci(x['amp_rel_err'])}$ (${sci(x['amp_budget'])}$)", f"${sci(x['deterministic_rel_err'])}$",
	f"${sci(x['q95_err'])}$", f"{100 * x['frac_within_target']:.1f}", f"{100 * x['median3_frac_within']:.1f}",
	f"{100 * x['single_shot_frac_within']:.1f}", f"{100 * b['validation']['frac_uniform_within']:.1f}",
	f"{100 * b['validation']['frac_at_X_within']:.1f}"]) for e, x, b in zip(EPS, E, B)])
s5 = load("step5")["grid"]
rows("TableG", [line(f"{nu:g}", [f"${sci(next(r for r in s5 if r['nu'] == nu and r['n'] == n)['rel_err_v'])}$" for n in (5, 6, 7, 8)]
	+ [f"${sci(min(next(r for r in s5 if r['nu'] == nu and r['n'] == 6)['p_neighbours']))}$"]) for nu in (0.2, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0)])
s7 = [r for r in load("step7")["rows"] if r["nu"] == 2.0 and r["n"] == 6 and r["p_a"] < 0 and r["safe_Lp"]]
pick = lambda k, e: next(r for r in s7 if r["kind"] == k and abs(r["eta_max"] - e) < 0.02 * e)
rows("TableP", [line(f"{e:g} ({pick('kink', e)['n_p']})", [f"${sci(pick(k, e)['amp_err'])}$ ({pick(k, e)['P_succ']:.2f})" for k in KIND]) for e in (5, 10, 20, 40)])
ge = load("grid_error")
sel = [(2.0, 5), (16.0, 6)]
rows("TableC", [line(f"{nu:g}; ${2**n}^2$; {r['h']}", [f"${sci(r['measured'])}$", f"${sci(r['leading'])}$", f"${sci(r['cancellation_free'])}$"])
	for nu, n in sel for r in ge if r["nu"] == nu and r["n"] == n])
big = [r for r in ge if r["nu"] >= 1]
macro("CfCount", len(big)); assert all(r["measured"] <= r["cancellation_free"] for r in big)
gt = load("grid_theorem")
for k, b in zip(("Zero", "One", "Two"), gt["B"]):
	macro(f"BK{k}", sci(b, 2))
rows("TableT", [line(f"{r['delta']:g}", [f"${sci(x)}$ (${sci(b)}$)" for x, b in zip(r["res"], r["bound"])]) for r in gt["orders"]])
o = gt["orders"]
rate = lambda k: math.log(o[-2]["res"][k] / o[-1]["res"][k]) / math.log(o[-2]["delta"] / o[-1]["delta"])
macro("OrderT", ", ".join(f"{rate(k):.2f}" for k in range(3)))
grid_estimates = load("grid_estimates")["rows"]
find = lambda a: next(c for c in grid_estimates if (c["nu"], c["n"], c["h"]) == (a["nu"], a["n"], a["h"]))
rows("TableD", [line(f"{100 * e:g}\\%", [f"{a['nu']:g}; ${2**a['n']}^2$; {a['h']}"] + [f"${sci(x)}$" for x in (c["measured"], c["cancellation_free"],
	c["remainder_lattice"], c["u_max"], c["numerical_error_estimate"])]) for e, a in zip(EPS, A) for c in [find(a)]])
for a in A:
	assert a["eps_grid"] >= find(a)["numerical_error_estimate"] and a["eps_grid"] >= find(a)["cancellation_free"]
available_estimates = [c for c in grid_estimates if c["numerical_error_estimate"] is not None]
macro("GridEstimateCount", len(available_estimates)); macro("GridEstimateTotal", len(grid_estimates))
macro("GridEstimateMedian", f"{sorted(c['numerical_error_estimate'] / c['measured'] for c in available_estimates)[len(available_estimates) // 2]:.1f}")
pts = [("0", "1"), ("0.2", "1"), ("0.4", "1"), ("0.6", "1"), ("0.8", "1"), ("0.4", "0.25"), ("0.4", "4")]
swp = []
for rho, lam in pts:
	d = R if (rho, lam) == ("0.4", "1") else R / f"rho{rho}_lam{lam}"
	a, b = at(load("step8", d)["best"], 0.01), at(load("step9_sampling", d)["best"], 0.01)
	c = next(c for c in load("grid_estimates", d)["rows"] if (c["nu"], c["n"], c["h"]) == (a["nu"], a["n"], a["h"]))
	assert a["eps_grid"] >= c["numerical_error_estimate"]  # sampled estimates are charged conditionally
	swp.append((a["T_total"], b["T_total"]))
	swp_line = line(f"{rho}; {lam}", [f"${sci(a['T_total'])}$", a["logical_qubits"], f"{a['nu']:g}; ${2**a['n']}^2$", f"${sci(b['T_total'])}$",
		b["logical_qubits"], f"{b['P_succ']:.2f}", f"{100 * b['validation']['frac_uniform_within']:.1f}"])
	ea = load("step10", d)["emulation"]
	swp[-1] += (swp_line, 100 * b["validation"]["frac_uniform_within"], 100 * next(e for e in ea if e["target"] == 0.01)["frac_within_target"])
rows("TableS", [x[2] for x in swp])
macro("SweepWithin", f"{min(x[3] for x in swp):.1f}"); macro("SweepWithinA", f"{min(x[4] for x in swp):.1f}")
macro("SweepAmin", sci(min(x[0] for x in swp))); macro("SweepAmax", sci(max(x[0] for x in swp)))
macro("SweepBmin", sci(min(x[1] for x in swp))); macro("SweepBmax", sci(max(x[1] for x in swp)))
sc = load("scaling_m")
m, fac = sc["rows"], sc["factors"]
macro("FacPlain", f"{fac['plain']:.2f}"); macro("FacMargin", f"{fac['margin']:.3f}"); macro("FacExact", f"{fac['exact']:.4f}")
macro("GrowthBTotal", f"{m[-1]['runs_B'] / m[0]['runs_B']:.0f}")
rows("TableM", [line(str(r["M"]), [f"${sci(r['P_succ_plain'])}$", f"{r['P_succ_shifted']:.2f}", f"{r['P_succ_exact']:.2f}",
	f"${sci(r['runs_A'])}$", f"${sci(r['runs_B_plain'])}$", f"${sci(r['runs_B'])}$"]) for r in m])
macro("GrowthA", f"{(m[-1]['runs_A'] / m[3]['runs_A']) ** (1 / (m[-1]['M'] - m[3]['M'])):.1f}")
c = load("compile_select") + [x for p in pts if p != ("0.4", "1") for x in load("compile_select", R / f"rho{p[0]}_lam{p[1]}")]
dq = [x["T_select_qualtran"] / x["T_select"] - 1 for x in c]  # Qualtran blocks relative to the compiled circuit
macro("CompileMin", f"{100 * max(dq):.0f}"); macro("CompileMax", f"{-100 * min(dq):.0f}"); macro("CompileCount", len(c))
g1 = lambda pa: next(x for x in s7all if x["kind"] == "smooth_w1" and x["nu"] == 2.0 and x["n"] == 6 and x["safe_Lp"] and abs(x["eta_max"] - 20) < 0.4 and abs(x["p_a"] - pa) < 1e-9)
s7all = load("step7")["rows"]
pa = min(x["p_a"] for x in s7all)
macro("WindowGain", f"{g1(pa)['P_succ'] / g1(0.0)['P_succ']:.1f}")
wc = [load("window_check", R if p == ("0.4", "1") else R / f"rho{p[0]}_lam{p[1]}") for p in pts]
assert all(min(w["E_min"]) > w["c"] for w in wc)
macro("UnsafeBench", f"{100 * wc[2]['unsafe_max']:.1f}"); macro("UnsafeMax", f"{100 * max(w['unsafe_max'] for w in wc):.0f}")
macro("WeightBench", sci(wc[2]["weight40"])); macro("WeightMax", sci(max(w["weight40"] for w in wc)))
one, three = [100 * x["single_shot_frac_within"] for x in s10], [100 * x["median3_frac_within"] for x in s10]
macro("SingleMin", f"{min(one):.0f}"); macro("SingleMax", f"{max(one):.0f}"); macro("ThreeMin", f"{min(three):.1f}")
low = [x["lower_bound"] for x in s10] + [b["validation"]["lower_bound"] for b in B]
macro("LowerBound", f"{100 * min(low):.1f}")
nq = {k: load(n) for k, n in ((0.1, "nonquadratic"), (0.5, "nonquadratic_kappa0.5"))}
d1 = lambda k, e, l: next((x for x in nq[k]["direct"] if x["target"] == e and x["levels"] == l), None)
rc = lambda k, e: next((x for x in nq[k]["recentred"] if x["target"] == e), None)
tab = []
for k in (0.1, 0.5):
	for e in (0.01, 0.001):
		a, a3, q = d1(k, e, 1), d1(k, e, 3), rc(k, e)
		tab.append(line(f"{k:g}; {100 * e:g}\\%", [f"${sci(a['nu'])}$", f"${sci(a['min_b'])}$", f"${sci(a['runs'])}$", f"${sci(a3['runs'])}$" if a3 else "--"]
			+ ([f"${sci(q['nu'])}$", f"${2**q['n']}^2$; {q['h']}", f"${sci(q['runs'])}$", f"${sci(q['T_total'])}$", q["logical_qubits"]] if q else ["--"] * 5)) if a else None)
rows("TableN", [t for t in tab if t])
macro("NqDirectOne", sci(d1(0.1, 0.01, 1)["runs"])); macro("NqRecOne", sci(rc(0.1, 0.01)["runs"]))
macro("NqDirectTenth", sci(d1(0.1, 0.001, 1)["runs"])); macro("NqRecTenth", sci(rc(0.1, 0.001)["runs"]))
macro("NqRichOne", sci(d1(0.1, 0.01, 3)["runs"])); macro("NqTOne", sci(rc(0.1, 0.01)["T_total"])); macro("NqQOne", rc(0.1, 0.01)["logical_qubits"])
macro("NqTTenth", sci(rc(0.1, 0.001)["T_total"])); macro("NqNuOne", f"{d1(0.1, 0.01, 1)['nu']:.2f}")
macro("NqRef", f"{100 * nq[0.1]['reference_error']:.1f}"); macro("NqRefHard", f"{100 * nq[0.5]['reference_error']:.1f}")
macro("NqSlope", f"{nq[0.1]['slope']:.3f}"); macro("NqW", f"{nq[0.1]['W0']:.1f}")
macro("NqLower", f"{100 * min(x['within_lower'] for v in nq.values() for x in v['recentred']):.1f}")
macro("NqBiasZero", sci(max(x["kappa0"] for x in nq[0.1]["bias"]), 0))
dp = load("direct_polynomial")
dpf = lambda ro, e: next(x for x in dp if x["readout"] == ro and x["target"] == e)
ratio = lambda x: f"{x:.0f}" if x >= 10 else f"{x:.1f}"
macro("DirA", ratio(dpf("A", 0.01)["ratio"])); macro("DirB", ratio(dpf("B", 0.01)["ratio"]))
macro("DirQA", dpf("A", 0.01)["qubits_schr"] - dpf("A", 0.01)["logical_qubits"]); macro("DirQB", dpf("B", 0.01)["qubits_schr"] - dpf("B", 0.01)["logical_qubits"])
rows("TableX", [line(f"{ro}, {100 * e:g}\\%", [x["steps"], f"({x['steps_schr']})", f"${sci(x['runs'])}$", f"(${sci(x['runs_schr'])}$)",
	f"${sci(x['T_total'])}$ ({ratio(x['ratio'])})", f"{x['logical_qubits']} ({x['qubits_schr']})"]) for ro in "AB" for e in EPS for x in [dpf(ro, e)]])
z = load("step10")["zero_risk"]
macro("ZeroRiskErr", sci(max(r["rel_err"] for r in z)))
bc = load("boundary_check")
macro("BoundaryDiff", sci(max(r["rel_diff"] for r in bc)))
(ROOT.parent / "paper" / "generated.tex").write_text("% Written by simulation/make_numbers.py; do not edit.\n" + "\n".join(out) + "\n")
print("\n".join(l for l in out if "TableA" not in l and l.startswith("\\newcommand") and "%\n" not in l))
