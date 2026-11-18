"""
white_dwarf_project.src
~~~~~~~~~~~~~~~~~~~~~~~
Core utility and diagnostic modules for the white dwarf planetary system simulation.
"""

from .diagnostics import (
    decompose_pericentre_difference,
    analyze_convergence_status,
)

__all__ = [
    "decompose_pericentre_difference",
    "analyze_convergence_status",
]
