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
zip -q -r "$out/src/anc/research-sources.zip" simulation lean review/validation/check_amplitude_shares.py review/validation/amplitude-shares.json review/validation/check_boundary_bound.py review/validation/check_boundary_bound_precision.py review/validation/boundary-bound.md review/validation/boundary-bound.json review/validation/boundary-bound-precision.json review/validation/check_sweep_readout.py review/validation/check_sweep_readout_all.py review/validation/sweep-results.json review/validation/sweep-results-before-refresh.json -x '*/__pycache__/*' 'lean/.lake/*' 'simulation/figures/*'
(cd "$out/src" && for i in 1 2 3; do pdflatex -interaction=nonstopmode -halt-on-error Paper.tex > /dev/null; done
	if grep -E "undefined|Citation.*undefined" Paper.log; then exit 1; fi
	grep "Output written" Paper.log
	rm -f Paper.aux Paper.log Paper.out Paper.pdf)
rm -f paper/arxiv-submission.zip
(cd "$out/src" && zip -q -r "$OLDPWD/paper/arxiv-submission.zip" .)
unzip -l paper/arxiv-submission.zip
rm -r "$out"
