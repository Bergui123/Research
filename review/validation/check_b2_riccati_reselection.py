"""Repeat the benchmark alternative search with the coarse B2=875 bound.

All modified inputs are temporary. This is a numerical resource experiment:
the box evolution, expansion and other circuit assumptions are not enclosed.
"""

from fractions import Fraction
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import check_alternative_boundary_search as alternative
from model import LAM, RHO, T

ROOT = Path(__file__).resolve().parents[2]


def main():
    assert (RHO, LAM, T) == (0.4, 1.0, 1.0)
    calculation = json.loads(subprocess.check_output(
        [sys.executable, str(ROOT / "review/validation/b2_riccati_majorant.py")], text=True
    ))
    saved = json.loads((ROOT / "review/validation/b2-riccati-majorant.json").read_text())
    assert calculation == saved
    assert Fraction(saved["B2_upper_rational"]) < 875
    source = alternative.RESULTS
    sampled = json.loads((source / "grid_estimates.json").read_text())["B2"]
    with tempfile.TemporaryDirectory(prefix="research-b2-reselection-") as directory:
        temporary = Path(directory)
        for name in ("step5.json", "step7.json", "step8.json"):
            shutil.copyfile(source / name, temporary / name)
        data = json.loads((source / "grid_estimates.json").read_text())
        data["B2"] = 875
        (temporary / "grid_estimates.json").write_text(json.dumps(data))
        alternative.RESULTS = temporary
        result = alternative.case(["0.05", "0.01", "0.001"], extension=True)
    # The reusable function labels its input as sampled; distinguish this run.
    result["B2_sampled_charge"] = sampled
    result["B2_applied"] = 875
    print(json.dumps({"scope": "Benchmark only; B2=875 and the boundary bound are analytic, other grid quantities remain numerical.",
                      "case": result}, indent=2))


if __name__ == "__main__":
    main()
