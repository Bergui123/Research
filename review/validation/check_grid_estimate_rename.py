"""Check that grid-estimate relabeling preserves all values and resource results.

The comparison baseline is the committed snapshot before the terminology fix.
This is a reproducibility check, not a numerical enclosure check.
"""

import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = "8e8ea0a"
SWEEPS = (
    "",
    "rho0_lam1",
    "rho0.2_lam1",
    "rho0.4_lam0.25",
    "rho0.4_lam4",
    "rho0.6_lam1",
    "rho0.8_lam1",
)


def git_json(path):
    raw = subprocess.check_output(["git", "show", f"{BASE}:{path}"], cwd=ROOT)
    return json.loads(raw)


def main():
    results = ROOT / "simulation" / "results"
    rows = {}
    for sweep in SWEEPS:
        prefix = f"simulation/results/{sweep}/" if sweep else "simulation/results/"
        old = git_json(prefix + "certified.json")
        current = json.loads((results / sweep / "grid_estimates.json").read_text())
        migrated = {
            **old,
            "rows": [
                {
                    ("numerical_error_estimate" if k == "certified_error" else k): v
                    for k, v in row.items()
                }
                for row in old["rows"]
            ],
        }
        assert current == migrated, f"grid estimates changed for {sweep or 'benchmark'}"

        old_step8 = git_json(prefix + "step8.json")
        current_step8 = json.loads((results / sweep / "step8.json").read_text())
        assert current_step8 == old_step8, f"resource selection changed for {sweep or 'benchmark'}"
        estimates = [r for r in current["rows"] if r["numerical_error_estimate"] is not None]
        rows[sweep or "benchmark"] = {
            "rows": len(current["rows"]),
            "available_estimates": len(estimates),
            "B2": current["B2"],
            "grid_values_unchanged": True,
            "step8_unchanged": True,
        }

    old_generated = subprocess.check_output(
        ["git", "show", f"{BASE}:paper/generated.tex"], cwd=ROOT, text=True
    )
    current_generated = (ROOT / "paper" / "generated.tex").read_text()
    old_generated = (
        old_generated.replace("\\CertCount", "\\GridEstimateCount")
        .replace("\\CertTotal", "\\GridEstimateTotal")
        .replace("\\CertMedian", "\\GridEstimateMedian")
    )
    assert current_generated == old_generated, "generated numerical macros or tables changed"

    report = {
        "baseline_commit": BASE,
        "datasets": rows,
        "generated_numbers_unchanged": True,
        "scope": "Terminology and schema only; numerical estimates are not validated enclosures.",
    }
    output = ROOT / "review" / "validation" / "grid-estimate-label-check.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
