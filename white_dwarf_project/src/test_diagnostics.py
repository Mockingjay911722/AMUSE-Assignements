"""
test_diagnostics.py
~~~~~~~~~~~~~~~~~~~
Test suite for sensitivity and convergence diagnostics module.
"""

import numpy as np
from white_dwarf_project.src.diagnostics import (
    decompose_pericentre_difference,
    analyze_convergence_status,
)


def test_pericentre_decomposition():
    # Synthetic time grid: 1000 points over 100,000 years
    t = np.linspace(0, 1e5, 1000)

    # Coarse run: expanding orbit + Lidov-Kozai oscillation
    a1 = 5.0 + 0.5 * (t / 1e5)
    e1 = 0.5 + 0.3 * np.sin(2 * np.pi * t / 2e4)

    # Fine run: identical physical envelope, but with a tiny Kozai period shift (phase drift)
    # and tiny mass interpolation difference in a
    a2 = a1 - 1e-4 * (t / 1e5)
    e2 = 0.5 + 0.3 * np.sin(2 * np.pi * t / (2e4 * 1.0005))

    decomp = decompose_pericentre_difference(t, a1, e1, a2, e2)

    # 1. Mathematical consistency: delta_q == dq_from_a + dq_from_e everywhere
    np.testing.assert_allclose(
        decomp["delta_q"],
        decomp["dq_from_a"] + decomp["dq_from_e"],
        atol=1e-15,
        err_msg="Pericentre decomposition identity violated!",
    )

    # 2. Check analysis report
    analysis = analyze_convergence_status(decomp, tolerance_au=1e-3)
    assert "dominant_discrepancy_source" in analysis
    assert analysis["dominant_discrepancy_source"] == "eccentricity_phase_shift"
    print("Test passed successfully!")
    print(f"Max delta_q: {analysis['max_abs_delta_q_au']:.6f} au")
    print(f"Dominant source: {analysis['dominant_discrepancy_source']}")
    print(f"Eccentricity vs semi-major axis ratio: {analysis['eccentricity_to_semimajor_ratio']:.1f}x")


if __name__ == "__main__":
    test_pericentre_decomposition()
