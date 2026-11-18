"""
diagnostics.py
~~~~~~~~~~~~~~
Numerical sensitivity and convergence diagnostics for coupled secular orbital
dynamics and stellar mass loss.

This module provides tools to decompose orbital pericentre discrepancies into:
1. Semi-major axis expansion deviations (induced by stellar mass interpolation)
2. Eccentricity excitation deviations (induced by Lidov-Kozai phase shifts)
"""

from typing import Dict, Any, Tuple
import numpy as np


def decompose_pericentre_difference(
    t: np.ndarray,
    a1: np.ndarray,
    e1: np.ndarray,
    a2: np.ndarray,
    e2: np.ndarray,
) -> Dict[str, np.ndarray]:
    """Decompose pericentre difference dq = q1 - q2 into contributions from a and e.

    Exact decomposition:
        q = a * (1 - e)
        q1 - q2 = (a1 - a2)*(1 - e2) - a1*(e1 - e2)
                = dq_from_a          + dq_from_e

    Parameters
    ----------
    t : np.ndarray
        Common time grid.
    a1, e1 : np.ndarray
        Semi-major axis (au) and eccentricity of run 1 (e.g. coarse refinement).
    a2, e2 : np.ndarray
        Semi-major axis (au) and eccentricity of run 2 (e.g. fine refinement).

    Returns
    -------
    dict
        Dictionary containing:
        - 'delta_q': total pericentre difference (q1 - q2)
        - 'dq_from_a': portion driven by semi-major axis difference
        - 'dq_from_e': portion driven by eccentricity difference
        - 'q1', 'q2': computed pericentre arrays
        - 'peak_index': index of maximum absolute difference
    """
    q1 = a1 * (1.0 - e1)
    q2 = a2 * (1.0 - e2)

    dq_from_a = (a1 - a2) * (1.0 - e2)
    dq_from_e = -a1 * (e1 - e2)
    delta_q = dq_from_a + dq_from_e

    peak_idx = int(np.argmax(np.abs(delta_q)))

    return {
        "time": t,
        "delta_q": delta_q,
        "dq_from_a": dq_from_a,
        "dq_from_e": dq_from_e,
        "q1": q1,
        "q2": q2,
        "peak_index": peak_idx,
        "max_abs_delta_q": float(np.abs(delta_q[peak_idx])),
        "max_abs_dq_from_a": float(np.max(np.abs(dq_from_a))),
        "max_abs_dq_from_e": float(np.max(np.abs(dq_from_e))),
    }


def analyze_convergence_status(
    decomp: Dict[str, Any],
    tolerance_au: float = 1e-3,
) -> Dict[str, Any]:
    """Evaluate whether a pericentre discrepancy is driven by secular phase drift
    or genuine secular envelope / mass divergence.

    Parameters
    ----------
    decomp : dict
        Output from `decompose_pericentre_difference`.
    tolerance_au : float
        Target tolerance limit in AU (default: 0.001 au).

    Returns
    -------
    dict
        Assessment results with status, dominant error source, and notes.
    """
    max_dq = decomp["max_abs_delta_q"]
    max_da = decomp["max_abs_dq_from_a"]
    max_de = decomp["max_abs_dq_from_e"]

    passed = max_dq <= tolerance_au
    dominant = "eccentricity_phase_shift" if max_de > max_da else "mass_interpolation"

    # Ratio of eccentricity contribution vs semi-major axis contribution
    ratio = max_de / max(max_da, 1e-15)

    return {
        "passed_tolerance": passed,
        "tolerance_au": tolerance_au,
        "max_abs_delta_q_au": max_dq,
        "dominant_discrepancy_source": dominant,
        "eccentricity_to_semimajor_ratio": float(ratio),
        "diagnosis": (
            "Discrepancy is within tolerance."
            if passed
            else (
                f"Discrepancy ({max_dq:.5f} au) exceeds tolerance ({tolerance_au} au), "
                f"primarily driven by {dominant} (ratio {ratio:.1f}x). "
                "Pointwise trajectory comparison is sensitive to slight secular phase shifts; "
                "verify whether envelope peaks and minimum clearances match before re-running fine tracks."
            )
        ),
    }
