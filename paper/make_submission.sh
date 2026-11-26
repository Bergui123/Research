#!/bin/sh
# Build the arXiv source package (paper/arxiv-submission.zip) and compile it from a clean directory.
set -e
cd "$(dirname "$0")/.."
out=$(mktemp -d)
mkdir -p "$out/src/figures" "$out/src/anc"
cp paper/Paper.tex paper/generated.tex paper/grid_theorem.tex paper/grid_box.tex paper/grid_proof.tex "$out/src/"
for f in $(grep -o 'includegraphics[^{]*{[^}]*}' paper/Paper.tex | sed 's/.*{//; s/}//'); do
	cp "simulation/figures/$f" "$out/src/figures/"
done
zip -q -r "$out/src/anc/research-sources.zip" simulation lean \
 review/validation/check_amplitude_shares.py review/validation/amplitude-shares.json \
 review/validation/check_boundary_bound.py review/validation/check_boundary_bound_precision.py \
 review/validation/boundary-bound.md review/validation/boundary-bound.json review/validation/boundary-bound-precision.json \
 review/validation/check_sweep_readout.py review/validation/check_sweep_readout_all.py \
 review/validation/sweep-results.json review/validation/sweep-results-before-refresh.json \
 review/validation/check_alternative_boundary_search.py review/validation/alternative-boundary-selection.json \
 review/validation/check_alternative_boundary_emulation.py review/validation/alternative-boundary-emulation.json \
 review/validation/b2_global_majorant.py review/validation/b2-global-majorant.json \
 review/validation/check_b2_global_majorant.py review/validation/b2-global-majorant-check.json \
 review/validation/check_b2_symbolic.py review/validation/b2-global-majorant-symbolic.json \
 review/validation/check_b2_budget_impact.py review/validation/b2-budget-impact.json \
 review/validation/b2-global-majorant.md \
 review/validation/generate_riccati_nodes.py review/validation/riccati-quad-nodes.json \
 review/validation/verify_riccati_nodes.py review/validation/riccati-residual-bound.json \
 review/validation/riccati-residual-bound.md review/validation/test_riccati_certificate.py \
 review/validation/check_riccati_residual_symbolic.py review/validation/riccati-residual-symbolic.json \
 review/validation/check_riccati_interpolation.py review/validation/riccati-interpolation-check.json \
 review/validation/b2_riccati_majorant.py review/validation/b2-riccati-majorant.json \
 review/validation/check_b2_riccati_symbolic.py review/validation/b2-riccati-symbolic.json \
 review/validation/check_b2_riccati_majorant.py review/validation/b2-riccati-check.json \
 review/validation/check_b2_riccati_reselection.py review/validation/b2-riccati-reselection.json \
 review/validation/check_b2_riccati_emulation.py review/validation/b2-riccati-emulation.json \
 review/validation/generate_polynomial_nodes.py review/validation/polynomial-quad-nodes.json.gz \
 review/validation/verify_polynomial_nodes.py review/validation/polynomial-residual-bound.json \
 review/validation/polynomial-residual-bound.md review/validation/test_polynomial_certificate.py \
 review/validation/check_polynomial_residual_symbolic.py review/validation/check_expansion_enclosures.py \
 review/validation/expansion-enclosures.json review/validation/check_expansion_against_float.py \
 review/validation/check_projected_box_bound.py review/validation/projected-box-bound.json \
 review/validation/projected-box-bound.md review/validation/test_projected_box_bound.py \
 review/validation/check_projected_box_numerical.py review/validation/projected-box-numerical.json \
 review/validation/check_projected_recost.py review/validation/projected-recost.json \
 review/validation/check_projected_emulation.py review/validation/projected-emulation.json \
 -x '*/__pycache__/*' 'lean/.lake/*' 'simulation/figures/*'
(cd "$out/src" && for i in 1 2 3; do pdflatex -interaction=nonstopmode -halt-on-error Paper.tex > /dev/null; done
	if grep -E "undefined|Citation.*undefined" Paper.log; then exit 1; fi
	grep "Output written" Paper.log
	rm -f Paper.aux Paper.log Paper.out Paper.pdf)
rm -f paper/arxiv-submission.zip
(cd "$out/src" && zip -q -r "$OLDPWD/paper/arxiv-submission.zip" .)
unzip -l paper/arxiv-submission.zip
rm -r "$out"
