"""Independent neutral polymer-network research models."""

from .model import (
    EquilibriumResult,
    EquilibriumRoot,
    NetworkParameters,
    ScanPoint,
    chemical_potential_residual,
    equilibrium_scan,
    free_energy_components,
    free_energy_j_m3_dry,
    residual_derivative_phi,
    slab_release_fraction,
    solve_equilibrium,
    time_to_release_fraction,
)

__version__ = "0.3.0"
__all__ = [
    "EquilibriumResult", "EquilibriumRoot", "NetworkParameters", "ScanPoint",
    "chemical_potential_residual", "equilibrium_scan", "free_energy_components",
    "free_energy_j_m3_dry", "residual_derivative_phi", "slab_release_fraction",
    "solve_equilibrium", "time_to_release_fraction",
]
