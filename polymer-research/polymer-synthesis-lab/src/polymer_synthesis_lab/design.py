"""Prospective single-step design and a locally fitted residual for one PEtOx protocol."""
import csv
import io
import math
from collections import Counter
from dataclasses import dataclass

import numpy as np

from .crop import HEADER, Recipe, parse_cohort

POOL = """M0_M,feed
2,30
2,60
2,120
3,30
3,60
3,120
4,30
4,60
4,120"""


def key(m, feed):
    return format(m, ".6g"), format(feed, ".6g")


def features(m, feed):
    return np.array([1.0, math.log2(m / 3.0), math.log2(feed / 60.0)])


def synthetic_surface():
    """Known log-linear response for a deterministic software control in the GUI."""
    output = io.StringIO()
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(HEADER)
    for study in ("A", "B"):
        for m in (1, 2, 4):
            for feed in (20, 40, 80):
                r = Recipe(m, m/feed, 100, 10, "Me/OH")
                factor = math.exp(features(m, feed) @ [0.10, 0.20, -0.30])
                writer.writerow([f"synthetic-control-{study}", f"M{m}-F{feed}", "control", "mock reactor",
                                 "synthetic SEC response", r.ends, m, r.I, r.T, r.minutes,
                                 factor*r.predict()[1], ""])
    return output.getvalue()


def candidates(raw, temperature, minutes, ends):
    if not isinstance(raw, str) or len(raw) > 20000:
        raise ValueError("The design pool requires CSV text under 20000 characters.")
    reader = csv.reader(io.StringIO(raw.strip()), strict=True)
    pool = {}
    try:
        if next(reader, None) != ["M0_M", "feed"]:
            raise ValueError("Candidate CSV must begin with M0_M,feed.")
        for index, row in enumerate(reader, 2):
            try:
                if len(row) != 2:
                    raise ValueError("expected two fields")
                m, feed = map(float, row)
                if not math.isfinite(feed) or feed <= 0:
                    raise ValueError("feed must be positive and finite")
                recipe = Recipe(m, m / feed, temperature, minutes, ends)
                if recipe.minutes <= 0:
                    raise ValueError("the chosen residence/hold time must be positive")
                coordinates = key(m, feed)
                if coordinates in pool:
                    raise ValueError("duplicate candidate concentration and feed")
                pool[coordinates] = recipe
            except (OverflowError, ValueError) as exc:
                raise ValueError(f"Candidate CSV line {index}: {exc}") from exc
            if len(pool) > 200:
                raise ValueError("Use up to 200 candidate recipes.")
    except csv.Error as exc:
        raise ValueError("Malformed candidate CSV.") from exc
    if not pool:
        raise ValueError("At least one candidate recipe is required.")
    return pool


def series_id(row):
    return row.paper, row.series


def matrices(rows):
    counts = Counter(series_id(row) for row in rows)
    weights = np.array([1.0 / counts[series_id(row)] for row in rows])
    x = np.array([features(row.recipe.M, row.recipe.feed) for row in rows]).reshape(-1, 3)
    return x, weights


@dataclass(frozen=True)
class Residual:
    coefficients: np.ndarray
    design_points: np.ndarray

    def predict(self, recipe):
        phi = features(recipe.M, recipe.feed)
        log_mn = math.log(recipe.predict()[1]) + phi @ self.coefficients
        # Extremely ill-scaled custom assays can give non-representable extrapolations.
        return math.exp(log_mn) if -700 < log_mn < 700 else None

    def interpolates(self, recipe):
        # All observations and predictions here share T and hold time. Test the
        # 2D convex hull of actual combinations, not independent min/max ranges.
        return inside_polygon(self.design_points, features(recipe.M, recipe.feed)[1:])


def hull(points):
    points = sorted(set(map(tuple, points)))
    def cross(a, b, c):
        return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])
    def half(sequence):
        result = []
        for p in sequence:
            while len(result) >= 2 and cross(result[-2], result[-1], p) <= 0:
                result.pop()
            result.append(p)
        return result
    return np.array(half(points)[:-1] + half(reversed(points))[:-1])


def inside_polygon(polygon, point):
    edge = np.roll(polygon, -1, axis=0) - polygon
    offset = point - polygon
    return bool(np.all(edge[:,0]*offset[:,1] - edge[:,1]*offset[:,0] >= -1e-10))


def fit(rows):
    if len({series_id(r) for r in rows}) < 4:
        return None, "At least four independent reaction series are required by this release."
    coordinates = [key(r.recipe.M, r.recipe.feed) for r in rows]
    if len(set(coordinates)) < 4 or min(len(set(values)) for values in zip(*coordinates)) < 2:
        return None, "Need varied concentration and feed with at least four distinct paired recipes."
    x, weights = matrices(rows)
    scale = np.sqrt(weights)
    xw = x * scale[:,None]
    if np.linalg.matrix_rank(xw) < 3 or np.linalg.cond(xw) > 1000:
        return None, "Concentration and feed effects are confounded or poorly conditioned."
    y = np.array([math.log(r.Mn/r.recipe.predict()[1]) for r in rows])
    coefficients = np.linalg.lstsq(xw, y*scale, rcond=None)[0]
    return Residual(coefficients, hull(x[:,1:])), "Local empirical fit to apparent SEC Mn."


def folds(rows, split):
    if split not in ("paper", "joint"):
        raise ValueError("Choose paper or joint holdout.")
    groups = {}
    for r in rows:
        group = r.paper if split == "paper" else key(r.recipe.M, r.recipe.feed)
        groups.setdefault(group, []).append(r)
    for group, test in groups.items():
        if split == "paper":
            train = [r for r in rows if r.paper != group]
        else:
            series = {series_id(r) for r in test}
            train = [r for r in rows if series_id(r) not in series and
                     all(a != b for a,b in zip(key(r.recipe.M,r.recipe.feed), group))]
        model, reason = fit(train)
        tested = []
        for r in test:
            baseline = r.recipe.predict()[1]
            estimate = model.predict(r.recipe) if model else None
            tested.append((r,baseline,estimate))
        yield group,len(train),tested,reason


def design(observations, requested, temperature, minutes):
    all_rows = parse_cohort(observations)
    rows = [r for r in all_rows if math.isclose(r.recipe.T, temperature, abs_tol=1e-7)
            and math.isclose(r.recipe.minutes, minutes, abs_tol=1e-7)]
    if not rows:
        raise ValueError("No measured observations match the selected T and hold time. Select a represented protocol.")
    pool = candidates(requested, temperature, minutes, rows[0].recipe.ends)
    observed = {key(r.recipe.M, r.recipe.feed) for r in rows}
    x, w = matrices(rows)
    information = np.eye(3) + (x.T*w) @ x
    model, reason = fit(rows)
    ranked = []
    for coordinates, recipe in pool.items():
        if coordinates in observed:
            continue
        phi = features(recipe.M, recipe.feed)
        gain = 0.5 * math.log1p(phi @ np.linalg.solve(information, phi))
        prediction = model.predict(recipe) if model and model.interpolates(recipe) else None
        ranked.append((recipe, gain, prediction))
    ranked.sort(key=lambda r:(-r[1],r[0].M,r[0].feed))
    return rows, ranked, model, reason
