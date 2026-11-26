"""Own mathematical implementations of constant-D transport in a full slab.

There is no dependency on an existing polymer modelling application. NumPy and
SciPy provide array operations, a generic stiff ODE integrator and optimization.
All lengths accepted by the scientific API are FULL slab thicknesses in metres.
"""

from __future__ import annotations

import csv
import io
import math
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, least_squares
from scipy.sparse import diags, lil_matrix
from scipy.special import erfc


MAX_POINTS = 1200
MAX_CELLS = 160
NUMERICAL_FO_CAP = 8.0
MODEL_ID = "constant-D-two-face-full-slab/1.0"


def finite_number(value, name: str, minimum=None, maximum=None) -> float:
    """Validate rather than silently accepting booleans, infinities or NaNs."""
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a finite number, not a boolean.")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite number.") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite.")
    if minimum is not None and result < minimum:
        raise ValueError(f"{name} must be at least {minimum:g}.")
    if maximum is not None and result > maximum:
        raise ValueError(f"{name} must be at most {maximum:g}.")
    return result


def integer_number(value, name: str, minimum: int, maximum: int) -> int:
    result = finite_number(value, name, minimum, maximum)
    if not result.is_integer():
        raise ValueError(f"{name} must be an integer.")
    return int(result)


def _physical_scales(diffusivity_m2_s, thickness_m) -> tuple[float, float]:
    d = finite_number(diffusivity_m2_s, "Diffusivity (m²/s)", 1e-24, 1e-4)
    length = finite_number(thickness_m, "Full slab thickness (m)", 1e-9, 1.0)
    return d, length


def _array(values, name: str, minimum=0.0, maximum=None) -> np.ndarray:
    try:
        raw = np.asarray(values)
        if raw.dtype.kind == "b":
            raise ValueError
        result = np.asarray(values, dtype=float)
    except (ValueError, TypeError, OverflowError) as exc:
        raise ValueError(f"{name} must contain finite numeric values.") from exc
    if result.ndim != 1 or not 1 <= result.size <= MAX_POINTS:
        raise ValueError(f"{name} must contain 1 to {MAX_POINTS} values.")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain finite values.")
    if np.any(result < minimum) or (maximum is not None and np.any(result > maximum)):
        raise ValueError(f"{name} contains values outside the supported range.")
    return result


def fraction_transferred_fourier(fourier_number) -> np.ndarray:
    """Fraction transferred with Fo=D*t/L² and FULL thickness L.

    An equivalent image series avoids subtracting two nearly equal numbers at
    short times. The Fourier series is rapidly convergent after Fo=0.05.
    The t=0 value is set from the initial condition, not a truncated series.
    """
    fo = _array(np.atleast_1d(fourier_number), "Fourier numbers", maximum=1e30)
    result = np.zeros_like(fo)
    short = (fo > 0) & (fo <= 0.05)
    if np.any(short):
        f = fo[short]
        root = np.sqrt(f / np.pi)
        answer = 4 * root
        # For Fo<=0.05, omitted k>=9 terms are smaller than exp(-405).
        for k in range(1, 9):
            sign = -1.0 if k % 2 else 1.0
            answer += sign * (
                8 * root * np.exp(-(k * k) / (4 * f))
                - 4 * k * erfc(k / (2 * np.sqrt(f)))
            )
        result[short] = answer
    long = fo > 0.05
    if np.any(long):
        odd = np.arange(1.0, 80.0, 2.0)
        remaining = (8 / np.pi**2) * np.sum(
            np.exp(-np.pi**2 * fo[long, None] * odd[None, :] ** 2)
            / odd[None, :] ** 2,
            axis=1,
        )
        result[long] = 1 - remaining
    return np.clip(result, 0, 1)


def fraction_transferred(time_s, diffusivity_m2_s, thickness_m) -> np.ndarray:
    """Normalized release or uptake; uniform initial state, two fixed faces."""
    d, length = _physical_scales(diffusivity_m2_s, thickness_m)
    times = _array(time_s, "Times (s)", maximum=1e12)
    return fraction_transferred_fourier(d * times / length**2)


def characteristic_times(diffusivity_m2_s, thickness_m) -> dict[str, float]:
    d, length = _physical_scales(diffusivity_m2_s, thickness_m)
    scale = length**2 / d
    result = {"diffusion_scale_s": scale, "first_mode_time_s": scale / np.pi**2}
    for fraction, label in [(0.5, "t50_s"), (0.9, "t90_s"), (0.99, "t99_s")]:
        root = brentq(
            lambda fo: float(fraction_transferred_fourier([fo])[0]) - fraction,
            0.0, 5.0, xtol=1e-14,
        )
        result[label] = float(root * scale)
    return result


@dataclass(frozen=True)
class FiniteVolumeResult:
    time_s: np.ndarray
    cell_centres_m: np.ndarray
    remaining_profiles: np.ndarray
    fraction_transferred: np.ndarray
    integrated_boundary_transfer: np.ndarray
    mass_balance_max_error: float
    time_tail_held_at_equilibrium: bool
    n_cells: int


def finite_volume_slab(
    time_s, diffusivity_m2_s, thickness_m, n_cells: int = 80
) -> FiniteVolumeResult:
    """Independent conservative spatial discretization plus BDF time solve.

    Cell centres lie half a cell from each Dirichlet boundary. A separate state
    integrates the flux leaving BOTH faces. It is not computed as 1-mean(c).
    At Fo>8 the negligible terminal state is held; no huge-time integration is
    attempted. The continuum remaining-mass bound there is exp(-8*pi²).
    """
    d, length = _physical_scales(diffusivity_m2_s, thickness_m)
    count = integer_number(n_cells, "Finite-volume cells", 8, MAX_CELLS)
    times = _array(time_s, "Times (s)", maximum=1e12)
    if np.any(np.diff(times) < 0):
        raise ValueError("Numerical output times must be in increasing order.")
    fo = d * times / length**2
    eval_fo, inverse = np.unique(np.minimum(fo, NUMERICAL_FO_CAP), return_inverse=True)
    initial = np.r_[np.ones(count), 0.0]
    if eval_fo[-1] == 0:
        states = np.repeat(initial[:, None], eval_fo.size, axis=1)
    else:
        diagonal = -2 * np.ones(count)
        diagonal[[0, -1]] = -3
        diffusion = count**2 * diags(
            [np.ones(count - 1), diagonal, np.ones(count - 1)], [-1, 0, 1],
            format="csc",
        )
        augmented = lil_matrix((count + 1, count + 1), dtype=float)
        augmented[:count, :count] = diffusion
        augmented[count, 0] = 2 * count
        augmented[count, count - 1] = 2 * count
        jacobian = augmented.tocsc()
        solved = solve_ivp(
            lambda _time, state: jacobian @ state,
            (0.0, float(eval_fo[-1])), initial, method="BDF", t_eval=eval_fo,
            jac=jacobian, rtol=2e-9, atol=2e-11,
        )
        if not solved.success:
            raise ValueError(f"Finite-volume integration failed: {solved.message}")
        states = solved.y
    profiles = states[:count, inverse].T
    outflux = states[count, inverse]
    balance = np.mean(profiles, axis=1) + outflux - 1
    if np.min(profiles) < -2e-8 or np.max(profiles) > 1 + 2e-8:
        raise ValueError("Numerical profiles exceeded the concentration bounds.")
    # Raw mass balance is calculated above before removing solver roundoff.
    display_profiles = np.clip(profiles, 0.0, 1.0)
    return FiniteVolumeResult(
        time_s=times.copy(),
        cell_centres_m=(np.arange(count) + 0.5) * length / count,
        remaining_profiles=display_profiles,
        fraction_transferred=np.clip(1 - np.mean(profiles, axis=1), 0, 1),
        integrated_boundary_transfer=outflux.copy(),
        mass_balance_max_error=float(np.max(np.abs(balance))),
        time_tail_held_at_equilibrium=bool(np.any(fo > NUMERICAL_FO_CAP)),
        n_cells=count,
    )


def parse_fraction_data(text: str) -> tuple[np.ndarray, np.ndarray]:
    """Read time_s,fraction CSV or two whitespace-separated columns."""
    if not isinstance(text, str) or not text.strip():
        raise ValueError("Fit data must contain time_s,fraction rows.")
    if len(text) > 100_000:
        raise ValueError("Fit data text is limited to 100,000 characters.")
    rows = []
    for line_number, line in enumerate(text.splitlines(), start=1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        cells = next(csv.reader(io.StringIO(line))) if "," in line else line.split()
        if len(cells) != 2:
            raise ValueError(f"Fit data line {line_number} must have two columns.")
        if not rows and [cell.strip().lower() for cell in cells] == ["time_s", "fraction"]:
            continue
        try:
            rows.append((float(cells[0]), float(cells[1])))
        except ValueError as exc:
            raise ValueError(f"Fit data line {line_number} is not numeric.") from exc
        if len(rows) > 1000:
            raise ValueError("Fit data is limited to 1000 observations.")
    if len(rows) < 5:
        raise ValueError("At least five observations are required for a D fit.")
    data = np.asarray(rows, dtype=float)
    times = _array(data[:, 0], "Fit times (s)", maximum=1e12)
    fractions = _array(data[:, 1], "Fit fractions", maximum=1.0)
    if np.unique(times[times > 0]).size < 3:
        raise ValueError("Fit data needs at least three distinct positive times.")
    order = np.argsort(times, kind="stable")
    return times[order], fractions[order]


@dataclass(frozen=True)
class DiffusivityFit:
    diffusivity_m2_s: float
    predicted: np.ndarray
    residuals: np.ndarray
    rmse: float
    log10_d_standard_error: float | None
    sensitivity_rms: float
    at_bound: bool
    warnings: tuple[str, ...]
    scan_log10_d: np.ndarray
    scan_sse: np.ndarray


def fit_diffusivity(
    time_s, fraction, thickness_m, lower_d: float = 1e-16, upper_d: float = 1e-8
) -> DiffusivityFit:
    """Fit log10(D) only. Thickness and normalization are held fixed.

    A bounded coarse scan initializes bounded least-squares, avoiding a flat
    plateau at a poor starting value. The objective scan is a sensitivity
    diagnostic, not a confidence interval. A local standard error is conditional
    on independent equal-variance fraction errors and an adequate physical model.
    """
    lower, length = _physical_scales(lower_d, thickness_m)
    upper = finite_number(upper_d, "Upper D bound (m²/s)", 1e-24, 1e-4)
    if upper <= lower:
        raise ValueError("Upper D bound must exceed lower D bound.")
    times = _array(time_s, "Fit times (s)", maximum=1e12)
    observed = _array(fraction, "Fit fractions", maximum=1.0)
    if times.size != observed.size or times.size < 5:
        raise ValueError("Fit requires at least five matching time/fraction pairs.")
    if np.unique(times[times > 0]).size < 3:
        raise ValueError("Fit needs at least three distinct positive times.")
    log_bounds = (math.log10(lower), math.log10(upper))

    def residual(log_d):
        return fraction_transferred(times, 10 ** float(log_d[0]), length) - observed

    grid = np.linspace(*log_bounds, 121)
    costs = np.array([np.sum(residual([v]) ** 2) for v in grid])
    initial = float(grid[int(np.argmin(costs))])
    # Bounded solver needs a feasible interior starting point.
    margin = (log_bounds[1] - log_bounds[0]) * 1e-6
    initial = float(np.clip(initial, log_bounds[0] + margin, log_bounds[1] - margin))
    solved = least_squares(
        residual, [initial], bounds=log_bounds, method="trf", jac="3-point",
        ftol=1e-12, xtol=1e-12, gtol=1e-12, max_nfev=200,
    )
    if not solved.success:
        raise ValueError(f"D fitting did not converge: {solved.message}")
    optimum = float(solved.x[0])
    d = 10**optimum
    residuals = residual([optimum])
    predicted = residuals + observed
    derivative = np.asarray(solved.jac)[:, 0]
    information = float(np.dot(derivative, derivative))
    sensitivity = float(np.sqrt(information / times.size))
    # A flat objective can stop at the deliberately interior initial point.
    # Include that feasibility margin in the boundary diagnostic, particularly
    # when the user searches many decades. Otherwise a lower-bound solution
    # can incorrectly look interior merely because its sensitivity is tiny.
    boundary_tolerance = max(1e-5, 2 * margin)
    at_bound = min(optimum - log_bounds[0], log_bounds[1] - optimum) <= boundary_tolerance
    warnings = []
    if at_bound:
        warnings.append("The optimum reaches a supplied D bound; D is not resolved within this search range.")
    if sensitivity < 0.01:
        warnings.append("Low sensitivity: changing D by a decade has little local effect on these data.")
    if np.all(predicted < 0.1) or np.all(predicted[times > 0] > 0.9):
        warnings.append("Only the early-time or near-equilibrium regime is sampled; extend the measurement window.")
    if np.any((times == 0) & (np.abs(observed) > 1e-6)):
        warnings.append("A nonzero fraction at t=0 conflicts with the assumed initial condition; no burst offset is fitted.")
    local_se = None
    if not at_bound and sensitivity >= 0.01:
        sse = float(np.sum(residuals**2))
        if sse > 1e-24:
            local_se = float(np.sqrt((sse / (times.size - 1)) / information))
        else:
            warnings.append("Near-zero synthetic or rounded residuals do not establish experimental uncertainty.")
    if local_se is not None and local_se > 0.3:
        warnings.append("Large local log10(D) standard error indicates a weakly constrained estimate.")
    return DiffusivityFit(
        diffusivity_m2_s=float(d), predicted=predicted, residuals=residuals,
        rmse=float(np.sqrt(np.mean(residuals**2))),
        log10_d_standard_error=local_se, sensitivity_rms=sensitivity,
        at_bound=at_bound, warnings=tuple(warnings), scan_log10_d=grid, scan_sse=costs,
    )


def steady_permeation(diffusivity_m2_s, thickness_m, solubility_mol_m3_Pa, delta_pressure_Pa):
    """Henry solution-diffusion relation with an explicit molar concentration basis."""
    d, length = _physical_scales(diffusivity_m2_s, thickness_m)
    solubility = finite_number(solubility_mol_m3_Pa, "Henry solubility (mol m⁻³ Pa⁻¹)", 1e-20, 1e6)
    delta_p = finite_number(delta_pressure_Pa, "Pressure difference (Pa)", -1e12, 1e12)
    permeability = d * solubility
    return {
        "permeability_mol_m_s_Pa": float(permeability),
        "steady_flux_mol_m2_s": float(permeability * delta_p / length),
    }
