"""Choose rational interpolant nodes; numerical integration is not the certificate."""
import gzip
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "simulation"))
from grid_theorem import solution  # noqa: E402

N = 4000
DENOMINATOR = 10**15
INDICES = [[(i, d - i) for d in range(0, 4 * k + 1, 2)
            for i in range(d + 1)] for k in (1, 2)]


def main():
    values = solution()(np.arange(2 * N + 1) / (2 * N)).reshape(2, 9, 9, -1)
    rows = [[[round(float(values[k, i, j, t]) * DENOMINATOR)
              for i, j in INDICES[k]] for k in range(2)]
            for t in range(2 * N + 1)]
    data = {"segments": N, "denominator": DENOMINATOR,
            "indices": INDICES, "nodes": rows}
    raw = json.dumps(data, separators=(",", ":")).encode()
    target = ROOT / "review/validation/polynomial-quad-nodes.json.gz"
    target.write_bytes(gzip.compress(raw, mtime=0))
    print(f"Wrote {target}: {len(raw)} uncompressed bytes")


if __name__ == "__main__":
    main()
