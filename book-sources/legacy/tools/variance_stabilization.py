#!/usr/bin/env python3
"""Compute exact Poisson sums for the variance-stabilization illustration."""
import math
from pathlib import Path

destination = Path(__file__).resolve().parents[1] / "figures/denoising/variance-stabilization-corrected.dat"
transforms = [lambda n: 2 * math.sqrt(n + 3 / 8),
              lambda n: math.sqrt(n) + math.sqrt(n + 1),
              lambda n: 2 * math.sqrt(n)]
rows = ["mean anscombe freeman square_root"]
for i in range(401):
    mean = i / 20
    probabilities = [math.exp(-mean)]
    for n in range(1, 241):
        probabilities.append(probabilities[-1] * mean / n)
    assert abs(math.fsum(probabilities) - 1) < 1e-13
    variances = []
    for transform in transforms:
        values = [transform(n) for n in range(len(probabilities))]
        expectation = math.fsum(p * v for p, v in zip(probabilities, values))
        variance = math.fsum(p * (v - expectation) ** 2 for p, v in zip(probabilities, values))
        variances.append(variance)
    rows.append(" ".join(f"{x:.12g}" for x in [mean, *variances]))
content = "\n".join(rows) + "\n"
if not destination.exists() or destination.read_text() != content:
    destination.write_text(content)
print(destination)
