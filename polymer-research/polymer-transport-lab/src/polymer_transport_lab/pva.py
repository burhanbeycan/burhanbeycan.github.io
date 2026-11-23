"""Water activity steps in fixed PVA films on a dry-polymer mass basis."""
from __future__ import annotations

import csv
import io
from dataclasses import dataclass
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

from .model import finite_number, fraction_transferred_fourier, integer_number

ISO_HEADER = "record_id,material_id,temperature_K,branch,aw,w_g_g"
KIN_HEADER = "curve_id,material_id,temperature_K,branch,thickness_um,faces,aw_start,aw_end,time_s,w_g_g"


def csv_rows(body, header, maximum):
    if not isinstance(body, str) or len(body) > 120000:
        raise ValueError("Enter CSV text of at most 120000 characters.")
    reader = csv.DictReader(io.StringIO(body.strip()))
    if reader.fieldnames != header.split(","):
        raise ValueError(f"Required CSV header: {header}")
    rows = list(reader)
    if not 1 <= len(rows) <= maximum:
        raise ValueError(f"Supply 1 to {maximum} rows.")
    if any(None in row or any(not isinstance(v, str) or not v.strip() for v in row.values()) for row in rows):
        raise ValueError("Each CSV row must contain exactly the listed columns with nonempty values.")
    return [{k: v.strip() for k, v in row.items()} for row in rows]


def condition(row):
    material = row["material_id"]
    temp = finite_number(row["temperature_K"], "Temperature", 250, 400)
    branch = row["branch"]
    if branch not in {"sorption", "desorption"}:
        raise ValueError("Use sorption or desorption as the isotherm branch.")
    return material, temp, branch


@dataclass
class Isotherm:
    condition: tuple
    record_ids: set
    aw: np.ndarray
    water: np.ndarray

    def step(self, start, end):
        start = finite_number(start, "Initial water activity", self.aw[0], self.aw[-1])
        end = finite_number(end, "Final water activity", self.aw[0], self.aw[-1])
        if (end <= start if self.condition[2] == "sorption" else end >= start):
            raise ValueError("The humidity step must follow the selected isotherm branch.")
        wi, wf = np.interp([start, end], self.aw, self.water)
        if abs(wf - wi) < 1e-8:
            raise ValueError("The selected activities must give distinct initial and final loadings.")
        return wi, wf


def read_isotherm(body):
    rows = csv_rows(body, ISO_HEADER, 200)
    conditions = {condition(row) for row in rows}
    if len(conditions) != 1:
        raise ValueError("Provide a single material grade/history, temperature and sorption branch.")
    points = sorted((finite_number(r["aw"], "Water activity", 0, 1),
                     finite_number(r["w_g_g"], "Equilibrium water, g/g dry polymer", 0, 10)) for r in rows)
    aw, water = np.array(points).T
    if len(aw) < 2 or np.any(np.diff(aw) <= 0) or np.any(np.diff(water) < 0) or water[-1] <= water[0]:
        raise ValueError("Supply at least two distinct activities, monotone loading and positive loading span.")
    return Isotherm(conditions.pop(), {r["record_id"] for r in rows}, aw, water)


@dataclass
class Curve:
    key: str
    thickness: float
    faces: int
    start: float
    end: float
    times: np.ndarray
    observed: np.ndarray


def read_curves(body, iso):
    groups = {}
    for row in csv_rows(body, KIN_HEADER, 1000):
        if condition(row) != iso.condition or row["curve_id"] in iso.record_ids:
            raise ValueError("Kinetics must match the isotherm condition and have independent record IDs.")
        key = row["curve_id"]
        thickness = finite_number(row["thickness_um"], "Full film thickness (µm)", 0.01, 10000)
        faces = integer_number(row["faces"], "Exposed faces", 1, 2)
        start = finite_number(row["aw_start"], "Initial activity", 0, 1)
        end = finite_number(row["aw_end"], "Final activity", 0, 1)
        iso.step(start, end)
        meta = (thickness, faces, start, end)
        previous, observations = groups.setdefault(key, (meta, []))
        if meta != previous:
            raise ValueError(f"Curve {key} changes thickness, geometry or humidity step.")
        observations.append((finite_number(row["time_s"], "Elapsed seconds", 0, 1e9),
                             finite_number(row["w_g_g"], "Observed water, g/g dry polymer", 0, 10)))
    if len(groups) > 20:
        raise ValueError("At most 20 curves per comparison.")
    result = []
    for key, (meta, observations) in groups.items():
        times, observed = np.array(sorted(observations)).T
        if np.sum(times > 0) < 3 or np.any(np.diff(times) <= 0):
            raise ValueError(f"Curve {key} needs three or more positive, distinct times.")
        result.append(Curve(key, *meta, times, observed))
    return result


def loading_curve(iso, curve, dref, contrast, cells=32):
    """Cell-centred conservative flux of the Kirchhoff integral of D(w).

    Model coordinate extends from the exposed surface to a reflecting plane.
    Its length is L for a backed film and L/2 for a film with equal exposed faces.
    """
    dref = finite_number(dref, "Mid-isotherm diffusivity", 1e-19, 1e-7)
    contrast = finite_number(contrast, "Log D contrast", -6, 6)
    cells = integer_number(cells, "Finite-volume cells", 12, 100)
    wi, wf = iso.step(curve.start, curve.end)
    path = finite_number(curve.thickness, "Full film thickness (µm)", 0.01, 10000) * 1e-6
    path /= integer_number(curve.faces, "Exposed faces", 1, 2)
    times = np.asarray(curve.times, dtype=float)
    if times.ndim != 1 or times.size == 0 or not np.isfinite(times).all() or (times < 0).any() or (np.diff(times) <= 0).any():
        raise ValueError("Times must be finite, nonnegative and strictly increasing.")
    fo = times * dref / path**2
    if contrast == 0:
        return wi + (wf - wi) * fraction_transferred_fourier(fo / 4)
    # w is kg water / kg dry polymer. Within a fixed dry coordinate and
    # constant dry mass per volume this is a local effective diffusivity.
    span = iso.water[-1] - iso.water[0]
    mid = 0.5 * (iso.water[-1] + iso.water[0])
    slope = contrast / span
    prefactor = math.exp(slope * (wi - mid))
    boundary_potential = prefactor * math.expm1(slope * (wf - wi)) / slope

    def rhs(_, water):
        phi = prefactor * np.expm1(slope * (water - wi)) / slope
        # q>0 points into the film. Dirichlet face is half a cell away.
        q = np.r_[2 * (boundary_potential - phi[0]), -np.diff(phi), 0.0]
        return -cells**2 * np.diff(q)

    if fo[-1] == 0:
        return np.array([wi])
    solution = solve_ivp(rhs, (0, fo[-1]), np.full(cells, wi), method="BDF", t_eval=fo,
                         rtol=3e-6, atol=max(1e-11, 1e-8 * abs(wf - wi)))
    if not solution.success:
        raise ValueError("Water-transport solver failed: " + solution.message)
    return solution.y.mean(axis=0)


def compare_models(iso, curves, holdout):
    holdout = finite_number(holdout, "Nominal held-out full thickness (µm)", 0.01, 10000)
    train = [c for c in curves if c.thickness != holdout]
    test = [c for c in curves if c.thickness == holdout]
    if not train or not test:
        raise ValueError("Use a measured nominal thickness for the holdout and provide other thicknesses for calibration.")
    if len({(c.start, c.end) for c in train}) < 2:
        raise ValueError("Use at least two independent humidity steps in training to resolve a D(w) contrast.")
    bounds = ([-19, -6], [-7, 6])
    # Equal weight per curve, in absolute dry-basis loading units.
    def residual(p):
        beta = p[1] if len(p) > 1 else 0
        return np.concatenate([(loading_curve(iso, c, 10**p[0], beta) - c.observed) / np.sqrt(c.times.size) for c in train])

    constant = least_squares(residual, [-12.0], bounds=([-19], [-7]), xtol=2e-7, gtol=1e-12)
    # Independent analytical constant-D and FV variable-D implementations are
    # deliberately distinct; avoid switching discretization inside the optimizer.
    def variable_residual(p):
        beta = p[1] if abs(p[1]) > 1e-10 else 1e-10
        return residual([p[0], beta])
    variable = least_squares(variable_residual, [constant.x[0], 0.2], bounds=bounds, xtol=2e-6,
                             ftol=2e-7, gtol=1e-11, diff_step=2e-4, max_nfev=45)
    if not constant.success or not variable.success:
        raise ValueError("The calibrations did not converge. Inspect sampling times and humidity-step coverage.")
    models = []
    for fit, beta in ((constant, 0.0), (variable, float(variable.x[1]))):
        dref = 10**float(fit.x[0])
        predicted = {c.key: loading_curve(iso, c, dref, beta).tolist() for c in curves}
        def rmse(curveset):
            return float(np.sqrt(np.mean([np.mean((predicted[c.key] - c.observed)**2) for c in curveset])))
        singular = np.linalg.svd(fit.jac, compute_uv=False)
        models.append(dict(dref=dref, contrast=beta, prediction=predicted, train_rmse=rmse(train), test_rmse=rmse(test),
                           weak_sensitivity=bool(singular[-1] < max(1e-10, 1e-4 * singular[0])),
                           near_bound=bool(any(abs(x-l)<.005 or abs(x-h)<.005 for x,l,h in zip(fit.x, *bounds)))))
    return models, train, test


def literature_rows():
    """Published apparent D and equilibrium loadings evaluated from published GAB coefficients."""
    rows = []
    for aw, tabulated in ((0.38, 0.94), (0.61, 1.40), (0.74, 2.40)):
        k_aw = 0.83 * aw
        w = 0.062 * 1.70 * k_aw / ((1.0 - k_aw) * (1.0 + 0.70 * k_aw))
        rows.append([aw, w, tabulated, tabulated * 1e-11])
    return rows


def apparent_d_ratio(log_contrast, outer_eta=8.0, tolerance=2e-5):
    """D_app/D_dry from the semi-infinite humidity-step similarity solution.

    With y=w/w_surface and eta=x/(2 sqrt(D_surface t)), j=(D/D_surface)y'.
    At early times the apparent constant-D slope corresponds to
    D_app / D_dry = pi exp(log_contrast) j(0)^2/4.
    """
    from scipy.integrate import solve_bvp
    from scipy.special import erfc

    if log_contrast == 0:
        return 1.0
    def ode(eta, state):
        y, flux = state
        # Newton iterates may leave [0,1] temporarily; bounded exponent avoids
        # overflow in trial steps. The converged physical solution is monotone.
        derivative = flux * np.exp(np.clip(log_contrast * (1 - y), -50, 50))
        return np.vstack((derivative, -2 * eta * derivative))
    def boundary(left, right):
        return [left[0] - 1, right[0]]
    eta = np.linspace(0, outer_eta, 121)
    guess = np.array([erfc(eta), -2 * np.exp(-eta**2) / np.sqrt(np.pi)])
    solution = solve_bvp(ode, boundary, eta, guess, tol=tolerance, max_nodes=5000)
    if not solution.success:
        raise ValueError("Could not resolve the semi-infinite similarity solution for apparent D.")
    return float(np.pi * np.exp(log_contrast) * solution.y[1,0]**2 / 4)


def paper_humidity_check():
    """Predefined endpoint fit with central activity reserved; data are published apparent D."""
    from scipy.optimize import brentq

    low, mid, high = literature_rows()
    ratio = high[3] / low[3]
    # Equilibrium endpoints come from the independently fitted GAB isotherm.
    # B = log[D(w_high)/D(0)], where each experiment starts in a dry film.
    def residual(contrast):
        return apparent_d_ratio(contrast) / apparent_d_ratio(contrast * low[1] / high[1]) - ratio
    contrast = brentq(residual, 0, 4)
    dry_d = high[3] / apparent_d_ratio(contrast)
    predicted = dry_d * apparent_d_ratio(contrast * mid[1] / high[1])
    # Both models fit D_app itself, in m²/s, with equal endpoint weights.
    return {"dry_d": dry_d, "contrast": contrast, "constant": (low[3] + high[3]) / 2,
            "predicted": predicted, "observed": mid[3], "aw": mid[0]}
