#!/bin/sh
# lambda / rho sweeps at the 1% target; each point writes results/rho<r>_lam<l>/.
set -e
py=../.venv/bin/python
for point in "0 1" "0.2 1" "0.6 1" "0.8 1" "0.4 0.25" "0.4 4"; do
	set -- $point
	export EXEC_RHO=$1 EXEC_LAM=$2 EXEC_TARGETS=0.01
	echo "== rho=$1 lambda=$2"
	$py step5_nu_sweep.py > /dev/null
	$py step7_schrodinger.py > /dev/null
	$py certify_grid.py
	$py step8_resources.py
	$py step9_sampling.py
	$py step10_validation.py
	$py compile_select.py
	$py check_window.py
done
