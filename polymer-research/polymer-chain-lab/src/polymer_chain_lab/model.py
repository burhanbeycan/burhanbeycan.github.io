"""Scientific kernels. No GUI, file, network, or workbench dependencies.

All lengths in these APIs are nm, energies kJ/mol, molar masses g/mol,
temperature K, and scattering vectors nm^-1. See docs/SCIENCE.md.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
import io
import math

import numpy as np
from scipy.special import logsumexp

R_KJ_MOL_K = 0.00831446261815324
MAX_COMPONENTS = 512
MAX_RIS_BEADS = 200_000
BASES = ("number", "mass", "number_log_density", "mass_log_density")


def finite_number(value, name, minimum=None, maximum=None):
    """Validate a finite scalar without accepting boolean or complex values."""
    if isinstance(value, (bool, np.bool_, complex)):
        raise ValueError(f"{name} must be a finite number.")
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite number.") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite.")
    if minimum is not None and number < minimum:
        raise ValueError(f"{name} must be at least {minimum:g}.")
    if maximum is not None and number > maximum:
        raise ValueError(f"{name} must be at most {maximum:g}.")
    return number


def integer(value, name, minimum, maximum):
    result = finite_number(value, name, minimum, maximum)
    if result != int(result):
        raise ValueError(f"{name} must be an integer.")
    return int(result)


def _array(values, name, minimum=0.0, maximum=1e100, strict=False):
    try:
        result = np.asarray(values, dtype=float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must contain numeric values.") from exc
    if result.ndim != 1 or result.size == 0 or result.size > MAX_COMPONENTS:
        raise ValueError(f"{name} needs 1 to {MAX_COMPONENTS} entries.")
    if not np.isfinite(result).all():
        raise ValueError(f"{name} must contain finite values.")
    if np.any(result <= minimum if strict else result < minimum) or np.any(result > maximum):
        condition = "greater than" if strict else "at least"
        raise ValueError(f"{name} must be {condition} {minimum:g} and at most {maximum:g}.")
    return result


@dataclass(frozen=True)
class MolarMassDistribution:
    """Discrete number distribution, including a quadrature-discretized curve."""

    molar_mass_g_mol: np.ndarray
    number_fraction: np.ndarray
    input_basis: str

    @property
    def mass_fraction(self):
        scaled = self.molar_mass_g_mol / np.max(self.molar_mass_g_mol)
        weights = self.number_fraction * scaled
        return weights / weights.sum()

    @property
    def scattering_fraction(self):
        """Zero-angle coherent intensity weights at equal contrast per mass."""
        scaled = self.molar_mass_g_mol / np.max(self.molar_mass_g_mol)
        weights = self.number_fraction * scaled**2
        return weights / weights.sum()

    def moments(self):
        mass_scale = float(np.max(self.molar_mass_g_mol))
        scaled = self.molar_mass_g_mol / mass_scale
        m1 = float(np.dot(self.number_fraction, scaled))
        m2 = float(np.dot(self.number_fraction, scaled**2))
        m3 = float(np.dot(self.number_fraction, scaled**3))
        return {"Mn": mass_scale * m1, "Mw": mass_scale * m2 / m1,
                "Mz": mass_scale * m3 / m2, "dispersity": m2 / m1**2}


def make_distribution(molar_mass_g_mol, amount, basis="number"):
    """Normalize discrete amounts or trapezoidal densities per log10(M).

    Discrete duplicate masses are combined. Density abscissae must be unique;
    integration covers only the supplied mass range, without tail extrapolation.
    """
    if basis not in BASES:
        raise ValueError(f"Distribution basis must be one of {', '.join(BASES)}.")
    mass = _array(molar_mass_g_mol, "Molar masses", 0, 1e12, strict=True)
    if np.any(mass < 1e-6):
        raise ValueError("Molar masses must be at least 1e-6 g/mol.")
    amounts = _array(amount, "Distribution amounts")
    if mass.shape != amounts.shape:
        raise ValueError("Molar masses and distribution amounts must have equal lengths.")
    if not np.any(amounts > 0):
        raise ValueError("At least one distribution amount must be positive.")
    order = np.argsort(mass)
    mass, amounts = mass[order], amounts[order] / np.max(amounts)
    if basis.endswith("log_density"):
        if mass.size < 2 or np.any(np.diff(mass) <= 0):
            raise ValueError("A log-density curve needs at least two distinct molar masses.")
        increments = np.diff(np.log10(mass))
        widths = np.empty_like(mass)
        widths[0], widths[-1] = increments[0] / 2, increments[-1] / 2
        if mass.size > 2:
            widths[1:-1] = (increments[:-1] + increments[1:]) / 2
        amounts = amounts * widths
    else:
        unique, inverse = np.unique(mass, return_inverse=True)
        combined = np.zeros_like(unique)
        np.add.at(combined, inverse, amounts)
        mass, amounts = unique, combined
    if basis.startswith("mass"):
        amounts = amounts / mass
    positive = amounts > 0
    mass, amounts = mass[positive], amounts[positive]
    number_fraction = amounts / amounts.sum()
    mass.setflags(write=False)
    number_fraction.setflags(write=False)
    return MolarMassDistribution(mass, number_fraction, basis)


def parse_distribution(text, basis="number"):
    """Read a bounded two-column CSV (or tab-separated table)."""
    if not isinstance(text, str) or not text.strip():
        raise ValueError("Enter a two-column molar-mass distribution.")
    if len(text) > 200_000:
        raise ValueError("Distribution input exceeds 200,000 characters.")
    lines = [line for line in text.splitlines() if line.strip() and not line.lstrip().startswith("#")]
    if not lines:
        raise ValueError("Distribution contains no data rows.")
    delimiter = "\t" if "\t" in lines[0] and "," not in lines[0] else ","
    rows = list(csv.reader(io.StringIO("\n".join(lines)), delimiter=delimiter))
    headers = {"m", "molar_mass", "molar_mass_g_mol", "m_g_mol", "molar mass", "mass", "mw"}
    if rows and rows[0] and rows[0][0].strip().lower() in headers:
        if len(rows[0]) != 2:
            raise ValueError("The CSV header must have exactly two columns.")
        rows = rows[1:]
    if not 1 <= len(rows) <= MAX_COMPONENTS:
        raise ValueError(f"Enter 1 to {MAX_COMPONENTS} distribution rows.")
    mass, amount = [], []
    for index, row in enumerate(rows, start=1):
        if len(row) != 2:
            raise ValueError(f"Distribution row {index} must have exactly two columns.")
        mass.append(finite_number(row[0].strip(), f"Molar mass in row {index}"))
        amount.append(finite_number(row[1].strip(), f"Amount in row {index}"))
    return make_distribution(mass, amount, basis)


def degrees_of_polymerization(molar_mass_g_mol, repeat_mass_g_mol,
                              end_group_mass_g_mol=0.0):
    """Effective DP from M=DP*M0+Mend for one common chain chemistry.

    Mend is the net terminal mass correction relative to DP repetitions of the
    selected constitutional repeat contribution. It is not a second copy of
    the atoms already included in M0. Fractional values are retained for
    continuous/distribution bins; no integer chain reconstruction is implied.
    """
    mass = _array(np.atleast_1d(molar_mass_g_mol), "Molar masses", 0, 1e12, strict=True)
    repeat = finite_number(repeat_mass_g_mol, "Repeat-unit molar mass", 0.001, 1e6)
    terminal = finite_number(end_group_mass_g_mol, "Net end-group mass correction", 0, 1e6)
    if np.any(mass <= terminal):
        raise ValueError("Every chain molar mass must exceed its net end-group mass correction.")
    return (mass - terminal) / repeat


def gaussian_dimensions(contour_length_nm, kuhn_length_nm):
    """Return ensemble mean-square end-to-end and gyration dimensions."""
    length = _array(np.atleast_1d(contour_length_nm), "Contour lengths", 0, 1e14, strict=True)
    kuhn = finite_number(kuhn_length_nm, "Kuhn length", 1e-8, 1e8)
    re2 = length * kuhn
    return re2, re2 / 6


def wormlike_dimensions(contour_length_nm, persistence_length_nm):
    """Continuum WLC moments, including stable stiff-chain limits."""
    length = _array(np.atleast_1d(contour_length_nm), "Contour lengths", 0, 1e14, strict=True)
    persistence = finite_number(persistence_length_nm, "Persistence length", 1e-8, 1e8)
    x = length / persistence
    re2, rg2 = np.empty_like(x), np.empty_like(x)
    small = x < 0.1
    # Taylor series eliminate catastrophic cancellation through the rod limit.
    re_series = np.zeros_like(x[small])
    rg_series = np.zeros_like(x[small])
    for k in range(10, -1, -1):
        re_series = re_series * x[small] + 2 * (-1)**k / math.factorial(k + 2)
        rg_series = rg_series * x[small] + 2 * (-1)**k / math.factorial(k + 4)
    re2[small] = length[small]**2 * re_series
    rg2[small] = length[small]**2 * rg_series
    y = x[~small]
    re2[~small] = 2 * persistence**2 * (y + np.expm1(-y))
    rg2[~small] = persistence**2 * (y / 3 - 1 + 2 / y + 2 * np.expm1(-y) / y**2)
    return re2, rg2


def debye_form_factor(q_nm_inverse, rg2_nm2):
    """Gaussian continuum P(q); rows correspond to q, columns to components."""
    q = _array(np.atleast_1d(q_nm_inverse), "Scattering vectors", 0, 1e8)
    rg2 = _array(np.atleast_1d(rg2_nm2), "Mean-square radii of gyration", 0, 1e24, strict=True)
    x = q[:, None]**2 * rg2[None, :]
    result = np.empty_like(x)
    small = x < 1e-3
    y = x[small]
    result[small] = 1 + y * (-1 / 3 + y * (1 / 12 + y * (-1 / 60 + y / 360)))
    y = x[~small]
    result[~small] = 2 * (y + np.expm1(-y)) / y**2
    return result


def population_form_factor(q_nm_inverse, rg2_nm2, distribution):
    rg2 = np.atleast_1d(rg2_nm2)
    if rg2.size != distribution.number_fraction.size:
        raise ValueError("Each distribution component needs one radius of gyration.")
    return debye_form_factor(q_nm_inverse, rg2) @ distribution.scattering_fraction


def ris_statistics(n_torsions, temperature_K, gauche_energy_kj_mol, opposite_gauche_penalty_kj_mol):
    """Exact open finite-chain RIS messages and marginals in log space.

    State order: trans, gauche+, gauche-. Each torsion contributes a state
    energy. Each consecutive unlike-gauche pair contributes the pair penalty.
    """
    count = integer(n_torsions, "Number of torsions", 1, 498)
    temperature = finite_number(temperature_K, "Temperature", 1, 5000)
    gauche = finite_number(gauche_energy_kj_mol, "Gauche energy", -50, 50)
    penalty = finite_number(opposite_gauche_penalty_kj_mol, "Opposite-gauche penalty", -50, 100)
    beta = 1 / (R_KJ_MOL_K * temperature)
    energies = np.array([0., gauche, gauche])
    pair = np.zeros((3, 3))
    pair[1, 2] = pair[2, 1] = penalty
    transfer = -beta * (energies[None, :] + pair)
    forward = np.empty((count, 3))
    backward = np.zeros((count, 3))
    forward[0] = -beta * energies
    for i in range(1, count):
        forward[i] = logsumexp(forward[i - 1, :, None] + transfer, axis=0)
    for i in range(count - 2, -1, -1):
        backward[i] = logsumexp(transfer + backward[i + 1, None, :], axis=1)
    log_partition = float(logsumexp(forward[-1]))
    log_marginals = forward + backward
    marginals = np.exp(log_marginals - logsumexp(log_marginals, axis=1, keepdims=True))
    if count > 1:
        pair_log = forward[:-1, :, None] + transfer[None, :, :] + backward[1:, None, :]
        pair_marginals = np.exp(pair_log - logsumexp(pair_log, axis=(1, 2), keepdims=True))
        opposite_fraction = float(np.mean(pair_marginals[:, 1, 2] + pair_marginals[:, 2, 1]))
    else:
        pair_marginals = np.empty((0, 3, 3))
        opposite_fraction = None
    return {"log_partition": log_partition, "marginals": marginals,
            "pair_marginals": pair_marginals, "opposite_gauche_fraction": opposite_fraction,
            "backward": backward, "transfer": transfer, "state_energies": energies,
            "beta": beta}


def sample_ris(n_bonds=80, bond_length_nm=0.154, bond_angle_deg=112.,
               temperature_K=298.15, gauche_energy_kj_mol=1.5,
               opposite_gauche_penalty_kj_mol=4., n_chains=200, seed=1729):
    """Independent exact finite-chain torsion samples, reconstructed in 3D.

    This is an equilibrium local torsion model, not molecular dynamics.
    There is no nonlocal excluded volume, force field, or chemistry inference.
    """
    bonds = integer(n_bonds, "Number of bonds", 3, 500)
    chains = integer(n_chains, "Number of independent chains", 2, 2000)
    random_seed = integer(seed, "Random seed", 0, 2**32 - 1)
    if chains * (bonds + 1) > MAX_RIS_BEADS:
        raise ValueError(f"RIS run exceeds {MAX_RIS_BEADS:,} total beads; reduce bonds or chains.")
    bond = finite_number(bond_length_nm, "Bond length", 1e-6, 100)
    angle = finite_number(bond_angle_deg, "Bond angle", 1, 179)
    stats = ris_statistics(bonds - 2, temperature_K, gauche_energy_kj_mol, opposite_gauche_penalty_kj_mol)
    rng = np.random.default_rng(random_seed)
    states = np.empty((chains, bonds - 2), dtype=np.int8)
    initial = -stats["beta"] * stats["state_energies"] + stats["backward"][0]
    initial = np.exp(initial - logsumexp(initial))
    states[:, 0] = np.minimum(np.sum(rng.random(chains)[:, None] > np.cumsum(initial), axis=1), 2)
    for i in range(1, bonds - 2):
        logits = stats["transfer"] + stats["backward"][i, None, :]
        probabilities = np.exp(logits - logsumexp(logits, axis=1, keepdims=True))
        conditional = probabilities[states[:, i - 1]]
        states[:, i] = np.minimum(np.sum(rng.random(chains)[:, None] > np.cumsum(conditional, axis=1), axis=1), 2)
    theta = math.radians(angle)
    coordinates = np.zeros((chains, bonds + 1, 3))
    coordinates[:, 1, 0] = bond
    coordinates[:, 2] = coordinates[:, 1] + bond * np.array([-math.cos(theta), math.sin(theta), 0.])
    torsions = np.array([math.pi, math.pi / 3, -math.pi / 3])
    for i in range(bonds - 2):
        a, b, c = coordinates[:, i], coordinates[:, i + 1], coordinates[:, i + 2]
        tangent = (c - b) / bond
        normal = np.cross(b - a, tangent)
        normal /= np.linalg.norm(normal, axis=1)[:, None]
        transverse = np.cross(normal, tangent)
        phi = torsions[states[:, i]]
        direction = (-math.cos(theta) * tangent + math.sin(theta) *
                     (np.cos(phi)[:, None] * transverse + np.sin(phi)[:, None] * normal))
        coordinates[:, i + 3] = c + bond * direction
    end = coordinates[:, -1] - coordinates[:, 0]
    re2 = np.einsum("ij,ij->i", end, end)
    centered = coordinates - coordinates.mean(axis=1, keepdims=True)
    rg2 = np.mean(np.sum(centered**2, axis=2), axis=1)
    observed = np.bincount(states.ravel(), minlength=3) / states.size
    opposite = ((states[:, :-1] == 1) & (states[:, 1:] == 2)) | ((states[:, :-1] == 2) & (states[:, 1:] == 1))
    return {"coordinates_nm": coordinates, "states": states, "re2_nm2": re2, "rg2_nm2": rg2,
            "mean_re2_nm2": float(re2.mean()), "mean_rg2_nm2": float(rg2.mean()),
            "sem_re2_nm2": float(re2.std(ddof=1) / math.sqrt(chains)),
            "sem_rg2_nm2": float(rg2.std(ddof=1) / math.sqrt(chains)),
            "finite_characteristic_ratio": float(re2.mean() / (bonds * bond**2)),
            "sampled_state_fraction": observed, "exact_state_fraction": stats["marginals"].mean(axis=0),
            "sampled_opposite_gauche_fraction": float(opposite.mean()) if opposite.size else None,
            "exact_opposite_gauche_fraction": stats["opposite_gauche_fraction"],
            "log_partition": stats["log_partition"]}
