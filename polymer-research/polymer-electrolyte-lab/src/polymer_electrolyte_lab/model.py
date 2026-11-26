"""Binary 1:1 electrolyte linear response in a polymer-fixed frame.

J_i = -sum_j L_ij grad(mu_tilde_j), with J in mol m^-2 s^-1 and
mu_tilde in J mol^-1. L is in mol^2 J^-1 m^-1 s^-1.
The neutral polymer is the velocity reference. This module neither
converts reference frames nor infers collective coefficients from chemistry.
"""
from __future__ import annotations
import math
import numpy as np

AVOGADRO = 6.02214076e23
ELEMENTARY_CHARGE = 1.602176634e-19
BOLTZMANN = 1.380649e-23
FARADAY = AVOGADRO * ELEMENTARY_CHARGE
GAS_CONSTANT = AVOGADRO * BOLTZMANN
PSD_RELATIVE_TOLERANCE = 8 * np.finfo(float).eps


def finite(value, name, *, lower=None, upper=None, strict=False):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a finite number.")
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError):
        raise ValueError(f"{name} must be a finite number.") from None
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite.")
    if lower is not None and (number <= lower if strict else number < lower):
        raise ValueError(f"{name} must be {'greater than' if strict else 'at least'} {lower}.")
    if upper is not None and number > upper:
        raise ValueError(f"{name} must not exceed {upper}.")
    return number


def nernst_einstein(concentration_mol_m3, d_plus, d_minus, temperature_K):
    c = finite(concentration_mol_m3, "Salt concentration (mol/m³)", lower=1e-6, upper=1e5)
    dp = finite(d_plus, "Cation self-diffusivity (m²/s)", lower=1e-30, upper=1e-4)
    dm = finite(d_minus, "Anion self-diffusivity (m²/s)", lower=1e-30, upper=1e-4)
    T = finite(temperature_K, "Temperature (K)", lower=1, upper=5000)
    lp = c * dp / (GAS_CONSTANT * T)
    lm = c * dm / (GAS_CONSTANT * T)
    return {
        "conductivity_S_m": FARADAY**2 * (lp + lm),
        "cation_current_fraction_ideal": dp / (dp + dm),
        "Lpp_ideal": lp, "Lmm_ideal": lm,
    }


def onsager(Lpp, Lmm, Lpm):
    pp = finite(Lpp, "L++", lower=0, upper=1)
    mm = finite(Lmm, "L−−", lower=0, upper=1)
    pm = finite(Lpm, "L+−", lower=-1, upper=1)
    bound = math.sqrt(pp) * math.sqrt(mm)
    # A zero diagonal forces the corresponding row/column to vanish exactly.
    # For positive diagonals allow only a few relative floating-point ULPs at
    # the PSD boundary. A dimensional absolute tolerance would admit genuinely
    # indefinite matrices when the physical coefficients are small.
    invalid_zero = (pp == 0 or mm == 0) and pm != 0
    invalid_bound = abs(pm) > bound and not math.isclose(
        abs(pm), bound, rel_tol=PSD_RELATIVE_TOLERANCE, abs_tol=0.0
    )
    if invalid_zero or invalid_bound:
        raise ValueError("The Onsager matrix must be positive semidefinite: |L+−| ≤ √(L++ L−−). Check units and reference frame.")
    scale = max(pp, mm)
    if scale == 0:
        return {"conductivity_S_m": 0.0, "t_plus_polymer_fixed": None,
                "minimum_eigenvalue": 0.0, "cross_correlation_ratio": None}
    # Scaling keeps the conditioning decision relative to the matrix itself,
    # including very small coefficients; it introduces no absolute L floor.
    p_scaled, m_scaled, cross_scaled = pp / scale, mm / scale, pm / scale
    charge_mobility_scaled = max(0.0, p_scaled + m_scaled - 2 * cross_scaled)
    sigma = (FARADAY**2 * scale) * charge_mobility_scaled
    t_plus = ((p_scaled - cross_scaled) / charge_mobility_scaled
              if charge_mobility_scaled > 1e-14 * (p_scaled + m_scaled) else None)
    eig = np.linalg.eigvalsh(np.array([[p_scaled, cross_scaled], [cross_scaled, m_scaled]], dtype=float)) * scale
    return {"conductivity_S_m": sigma, "t_plus_polymer_fixed": t_plus,
            "minimum_eigenvalue": float(eig[0]),
            "cross_correlation_ratio": max(-1.0, min(1.0, pm / bound)) if bound else None}

