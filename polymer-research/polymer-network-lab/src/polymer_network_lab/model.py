"""Neutral, dry-reference polymer-network thermodynamics and slab diffusion.

This module is independent of the workbench, GUI and other applications.  The
elastic convention is explicitly the Gaussian/neo-Hookean form used by Hong
et al. (2008, 2009): W_el = nu_e R T (I_1 - 3 - 2 log J) / 2.
See docs/SCIENCE.md for the reference state, derivation and limitations.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from math import isfinite, pi, sqrt
from numbers import Real
from typing import Iterable

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit, logit


GAS_CONSTANT_J_MOL_K = 8.31446261815324
PARAMETER_LIMITS = {
    "chi": (-5.0, 5.0),
    "nu_e_mol_m3": (1.0e-6, 1.0e6),
    "solvent_molar_volume_m3_mol": (1.0e-8, 1.0e-3),
    "temperature_K": (1.0, 3000.0),
}
PHI_MIN = 1.0e-14
PHI_MAX = 1.0 - 1.0e-14
MAX_SCAN_POINTS = 160
MAX_CURVE_POINTS = 1200


def finite_number(value: object, name: str, minimum: float, maximum: float) -> float:
    """Validate a scalar without accepting booleans, strings, NaN or infinity."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite number.")
    result = float(value)
    if not isfinite(result) or not minimum <= result <= maximum:
        raise ValueError(f"{name} must be between {minimum:g} and {maximum:g}.")
    return result


def integer_in_range(value: object, name: str, minimum: int, maximum: int) -> int:
    result = finite_number(value, name, minimum, maximum)
    if not result.is_integer():
        raise ValueError(f"{name} must be an integer.")
    return int(result)


@dataclass(frozen=True)
class NetworkParameters:
    """All densities refer to *dry polymer volume*, not swollen gel volume.

    nu_e_mol_m3 is a user-supplied effective elastic parameter in the stated
    model, never inferred from crosslinker feed, nominal conversion or MWD.
    chi is constant in composition and temperature during one calculation.
    """

    chi: float = 0.35
    nu_e_mol_m3: float = 100.0
    solvent_molar_volume_m3_mol: float = 18.0e-6
    temperature_K: float = 298.15

    def __post_init__(self) -> None:
        for name, bounds in PARAMETER_LIMITS.items():
            object.__setattr__(
                self, name, finite_number(getattr(self, name), name, *bounds)
            )

    @property
    def elastic_to_mixing_ratio(self) -> float:
        return self.nu_e_mol_m3 * self.solvent_molar_volume_m3_mol

    @property
    def energy_scale_j_m3_dry(self) -> float:
        return GAS_CONSTANT_J_MOL_K * self.temperature_K / self.solvent_molar_volume_m3_mol

    @property
    def dry_shear_modulus_pa(self) -> float:
        return self.nu_e_mol_m3 * GAS_CONSTANT_J_MOL_K * self.temperature_K


@dataclass(frozen=True)
class EquilibriumRoot:
    polymer_volume_fraction: float
    volume_swelling_ratio: float
    isotropic_stretch: float
    free_energy_j_m3_dry: float
    dimensionless_free_energy: float
    chemical_potential_residual_over_rt: float
    relative_residual: float
    curvature_dimensionless: float
    stability: str
    swollen_shear_modulus_pa: float


@dataclass(frozen=True)
class EquilibriumResult:
    parameters: NetworkParameters
    roots: tuple[EquilibriumRoot, ...]
    selected_index: int
    search_phi_min: float = PHI_MIN
    search_phi_max: float = PHI_MAX

    @property
    def selected(self) -> EquilibriumRoot:
        return self.roots[self.selected_index]


@dataclass(frozen=True)
class ScanPoint:
    value: float
    polymer_volume_fraction: float
    volume_swelling_ratio: float
    isotropic_stretch: float
    swollen_shear_modulus_pa: float
    stable_root_count: int


def _finite_array(values: object, name: str, low: float, high: float) -> np.ndarray:
    try:
        raw = np.asarray(values)
        if raw.dtype.kind not in "iuf":
            raise ValueError(f"{name} must contain numbers, not text or booleans.")
        array = raw.astype(float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must contain finite numbers.") from exc
    if array.size == 0:
        raise ValueError(f"{name} cannot be empty.")
    if array.size > MAX_CURVE_POINTS:
        raise ValueError(f"{name} is limited to {MAX_CURVE_POINTS} values.")
    if array.ndim > 1 or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be a scalar or a one-dimensional finite array.")
    if np.any(array < low) or np.any(array > high):
        raise ValueError(f"{name} must lie between {low:g} and {high:g}.")
    return array


def _return_shape(value: np.ndarray) -> float | np.ndarray:
    return float(value) if value.ndim == 0 else value


def _mixing_mu(phi: np.ndarray, chi: float) -> np.ndarray:
    """log(1-phi) + phi + chi phi^2 without dilute-limit cancellation."""
    dilute = (chi - 0.5) * phi**2
    for order in range(3, 12):
        dilute = dilute - phi**order / order
    return np.where(phi < 1.0e-4, dilute, np.log1p(-phi) + phi + chi * phi**2)


def _mu_unchecked(phi: np.ndarray, parameters: NetworkParameters) -> np.ndarray:
    a = parameters.elastic_to_mixing_ratio
    return _mixing_mu(phi, parameters.chi) + a * (np.cbrt(phi) - phi)


def chemical_potential_residual(
    polymer_volume_fraction: float | Iterable[float], parameters: NetworkParameters
) -> float | np.ndarray:
    """Return (mu_gel - mu_pure_solvent)/(RT), at free isotropic swelling.

    An equilibrium root is a stationary point; its stability is assessed from
    the *same* free energy by solve_equilibrium, not from residual size alone.
    """
    phi = _finite_array(polymer_volume_fraction, "polymer volume fraction", PHI_MIN, PHI_MAX)
    return _return_shape(_mu_unchecked(phi, parameters))


def residual_derivative_phi(
    polymer_volume_fraction: float | Iterable[float], parameters: NetworkParameters
) -> float | np.ndarray:
    phi = _finite_array(polymer_volume_fraction, "polymer volume fraction", PHI_MIN, PHI_MAX)
    derivative = phi * (2.0 * parameters.chi - 1.0 / (1.0 - phi))
    derivative += parameters.elastic_to_mixing_ratio * (phi ** (-2.0 / 3.0) / 3.0 - 1.0)
    return _return_shape(derivative)


def free_energy_components(
    volume_swelling_ratio: float | Iterable[float], parameters: NetworkParameters
) -> dict[str, float | np.ndarray]:
    """Free energies (J/m^3 of dry network), with the dry state set to zero.

    Q = V_swollen / V_dry >= 1.  The mixing energy is shifted by a constant
    relative to Hong et al.; derivatives and equilibrium are unaffected.
    """
    q = _finite_array(volume_swelling_ratio, "volume swelling ratio", 1.0, 1.0e14)
    flat = q.reshape(-1)
    mixing = np.zeros_like(flat)
    wet = flat > 1.0
    phi = 1.0 / flat[wet]
    mixing[wet] = (flat[wet] - 1.0) * np.log1p(-phi) + parameters.chi * (1.0 - phi)
    logq = np.log(flat)
    x = (2.0 / 3.0) * logq
    # Near Q=1, use exp(x)-1-x's series so the O(log(Q)^2) elastic energy is
    # not lost by subtracting two O(log(Q)) terms.
    series = x**2 / 2.0
    term = series.copy()
    for order in range(3, 12):
        term = term * x / order
        series = series + term
    elastic = 1.5 * parameters.elastic_to_mixing_ratio * np.where(
        x < 1.0e-4, series, np.expm1(x) - x
    )
    scale = parameters.energy_scale_j_m3_dry
    mixing = (scale * mixing).reshape(q.shape)
    elastic = (scale * elastic).reshape(q.shape)
    return {
        "mixing": _return_shape(mixing),
        "elastic": _return_shape(elastic),
        "total": _return_shape(mixing + elastic),
    }


def free_energy_j_m3_dry(
    volume_swelling_ratio: float | Iterable[float], parameters: NetworkParameters
) -> float | np.ndarray:
    return free_energy_components(volume_swelling_ratio, parameters)["total"]


def _stationary_residual_points(parameters: NetworkParameters) -> list[float]:
    """Candidate extrema of the residual from an exact degree-eight polynomial.

    With y = phi^(1/3), multiply dg/dphi by 3 y^2(1-y^3).  The resulting
    polynomial is a - 3a y^2 - a y^3 + (3a+6chi-3)y^5 - 6chi y^8.
    It supplies monotonic-interval boundaries in addition to the logit grid.
    """
    a = parameters.elastic_to_mixing_ratio
    chi = parameters.chi
    coefficients = np.array([a, 0.0, -3.0 * a, -a, 0.0, 3.0 * a + 6.0 * chi - 3.0,
                             0.0, 0.0, -6.0 * chi])
    roots = np.polynomial.polynomial.polyroots(coefficients)
    result = []
    for root in roots:
        if abs(float(root.imag)) < 1.0e-8 and 0.0 < float(root.real) < 1.0:
            phi = float(root.real) ** 3
            if PHI_MIN < phi < PHI_MAX:
                result.append(phi)
    return result


def solve_equilibrium(parameters: NetworkParameters | None = None) -> EquilibriumResult:
    """Enumerate stationary states and select the lowest stable free energy.

    The problem is homogeneous, neutral, free isotropic swelling in a pure
    solvent reservoir (mu_reservoir = 0).  This is not an inhomogeneous phase
    coexistence, constrained swelling, ionic, or dynamic swelling solver.
    """
    if parameters is None:
        parameters = NetworkParameters()
    if not isinstance(parameters, NetworkParameters):
        raise ValueError("parameters must be a NetworkParameters instance.")
    # Logit sampling resolves both very dilute and nearly dry states.  Residual
    # extrema supplement the grid so two crossings between grid nodes need not
    # be lost.  Bracketing is solved in logit coordinates for relative accuracy.
    critical = _stationary_residual_points(parameters)
    grid = np.unique(np.concatenate((
        expit(np.linspace(logit(PHI_MIN), logit(PHI_MAX), 1025)),
        np.asarray(critical), np.array([PHI_MIN, PHI_MAX]),
    )))
    values = _mu_unchecked(grid, parameters)
    candidates: list[float] = []
    for index in range(len(grid) - 1):
        left, right = float(grid[index]), float(grid[index + 1])
        fleft, fright = float(values[index]), float(values[index + 1])
        if fleft == 0.0:
            candidates.append(left)
        if np.signbit(fleft) != np.signbit(fright):
            root_logit = brentq(
                lambda z: float(_mu_unchecked(np.asarray(expit(z)), parameters)),
                float(logit(left)), float(logit(right)),
                xtol=5.0e-13, rtol=1.0e-14, maxiter=150,
            )
            candidates.append(float(expit(root_logit)))
    if values[-1] == 0.0:
        candidates.append(float(grid[-1]))
    # A tangent root has no sign change.  Include it only when the residual is
    # small relative to the two physical contributions, not to an arbitrary 1.
    for phi in critical:
        mixing = float(_mixing_mu(np.asarray(phi), parameters.chi))
        elastic = parameters.elastic_to_mixing_ratio * (np.cbrt(phi) - phi)
        scale = max(abs(mixing), abs(float(elastic)), np.finfo(float).tiny)
        if abs(mixing + elastic) <= 1.0e-10 * scale:
            candidates.append(phi)
    candidates.sort()
    unique: list[float] = []
    for phi in candidates:
        if not unique or abs(float(logit(phi) - logit(unique[-1]))) > 1.0e-8:
            unique.append(phi)
    if not unique:
        raise ValueError("No equilibrium root was resolved within the supported polymer-volume-fraction range.")

    roots: list[EquilibriumRoot] = []
    for phi in unique:
        q = 1.0 / phi
        stretch = q ** (1.0 / 3.0)
        mixing = float(_mixing_mu(np.asarray(phi), parameters.chi))
        elastic = parameters.elastic_to_mixing_ratio * (np.cbrt(phi) - phi)
        residual = mixing + float(elastic)
        # Normalize by the absolute terms of the equation, rather than by the
        # already-cancelled mixing sum.  In a poor solvent at vanishing elastic
        # density, that sum can be smaller than floating-point roundoff even
        # while phi and the absolute residual are well resolved.
        residual_scale = abs(float(np.log1p(-phi))) + phi + abs(parameters.chi * phi**2) + abs(float(elastic))
        relative = abs(residual) / max(residual_scale, np.finfo(float).tiny)
        dm = phi * (2.0 * parameters.chi - 1.0 / (1.0 - phi))
        de = parameters.elastic_to_mixing_ratio * (phi ** (-2.0 / 3.0) / 3.0 - 1.0)
        curvature = -phi**2 * (dm + de)  # d^2 [W/(RT/Vs)] / dQ^2
        tolerance = 1.0e-10 * phi**2 * (abs(dm) + abs(de))
        stability = "stable" if curvature > tolerance else "unstable" if curvature < -tolerance else "marginal"
        free_energy = float(free_energy_j_m3_dry(q, parameters))
        roots.append(EquilibriumRoot(
            polymer_volume_fraction=phi,
            volume_swelling_ratio=q,
            isotropic_stretch=stretch,
            free_energy_j_m3_dry=free_energy,
            dimensionless_free_energy=free_energy / parameters.energy_scale_j_m3_dry,
            chemical_potential_residual_over_rt=residual,
            relative_residual=relative,
            curvature_dimensionless=curvature,
            stability=stability,
            swollen_shear_modulus_pa=parameters.dry_shear_modulus_pa / stretch,
        ))
    stable = [index for index, root in enumerate(roots) if root.stability == "stable"]
    if not stable:
        raise ValueError("No locally stable homogeneous equilibrium was resolved for these parameters.")
    selected = min(stable, key=lambda index: roots[index].dimensionless_free_energy)
    return EquilibriumResult(parameters, tuple(roots), selected)


def equilibrium_scan(
    parameters: NetworkParameters, variable: str, values: Iterable[float]
) -> tuple[ScanPoint, ...]:
    if variable not in ("chi", "nu_e_mol_m3"):
        raise ValueError("scan variable must be chi or nu_e_mol_m3.")
    values_array = _finite_array(values, "scan values", *PARAMETER_LIMITS[variable])
    if values_array.ndim != 1 or not 2 <= len(values_array) <= MAX_SCAN_POINTS:
        raise ValueError(f"A scan must contain 2 to {MAX_SCAN_POINTS} values.")
    if np.any(np.diff(values_array) <= 0):
        raise ValueError("Scan values must be strictly increasing.")
    result = []
    for value in values_array:
        equilibrium = solve_equilibrium(replace(parameters, **{variable: float(value)}))
        root = equilibrium.selected
        result.append(ScanPoint(
            float(value), root.polymer_volume_fraction, root.volume_swelling_ratio,
            root.isotropic_stretch, root.swollen_shear_modulus_pa,
            sum(candidate.stability == "stable" for candidate in equilibrium.roots),
        ))
    return tuple(result)


def _release_from_fourier_number(fo: np.ndarray) -> np.ndarray:
    """Two-face slab solution with full thickness in Fo=D*t/L^2.

    Below Fo=0.005 the first boundary-overlap correction to the short-time
    expression is below 1e-20.  Above it, 128 odd Fourier modes are ample.
    This avoids cancellation at very small nonzero times.
    """
    shape = fo.shape
    flat = fo.reshape(-1)
    result = np.empty_like(flat)
    early = flat <= 0.005
    result[early] = 4.0 * np.sqrt(flat[early] / pi)
    odd = 2.0 * np.arange(128, dtype=float) + 1.0
    if np.any(~early):
        exponent = -pi**2 * flat[~early, None] * odd[None, :] ** 2
        remaining = (8.0 / pi**2) * np.sum(np.exp(exponent) / odd[None, :] ** 2, axis=1)
        result[~early] = 1.0 - remaining
    return np.clip(result, 0.0, 1.0).reshape(shape)


def slab_release_fraction(
    times_s: float | Iterable[float], full_thickness_m: float, diffusivity_m2_s: float
) -> float | np.ndarray:
    """Fraction released from an already equilibrated, uniformly loaded slab.

    Perfect sinks on BOTH faces; a constant user-supplied effective solute D;
    fixed full thickness L.  No dynamic swelling, binding, erosion or external
    mass-transfer resistance is included.  D is not calculated from nu_e.
    """
    thickness = finite_number(full_thickness_m, "full slab thickness (m)", 1.0e-12, 1.0e3)
    diffusivity = finite_number(diffusivity_m2_s, "diffusivity (m^2/s)", 1.0e-30, 1.0e-3)
    times = _finite_array(times_s, "times (s)", 0.0, 1.0e15)
    return _return_shape(_release_from_fourier_number(diffusivity * times / thickness**2))


def time_to_release_fraction(
    fraction: float, full_thickness_m: float, diffusivity_m2_s: float
) -> float:
    """Time (s) at a target released fraction; independent of plot duration."""
    fraction = finite_number(fraction, "target released fraction", 1.0e-6, 1.0 - 1.0e-6)
    thickness = finite_number(full_thickness_m, "full slab thickness (m)", 1.0e-12, 1.0e3)
    diffusivity = finite_number(diffusivity_m2_s, "diffusivity (m^2/s)", 1.0e-30, 1.0e-3)
    fourier_number = brentq(
        lambda fo: float(_release_from_fourier_number(np.asarray(fo))) - fraction,
        0.0, 2.0, xtol=1.0e-15, rtol=1.0e-14,
    )
    return fourier_number * thickness**2 / diffusivity
