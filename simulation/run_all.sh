#!/bin/sh
# The whole numerical pipeline, in order (about an hour on 16 cores).
set -e
py=../.venv/bin/python
for s in step3_benchmark step4_pde step5_nu_sweep boundary_check grid_error grid_theorem chebyshev step7_schrodinger \
	certify_grid step8_resources step9_sampling step10_validation direct_polynomial check_window compile_select scaling_m; do
	echo "== $s"
	$py $s.py
done
echo "== nonquadratic"
$py nonquadratic.py
EXEC_KAPPA=0.5 $py nonquadratic.py
sh run_sweeps.sh
