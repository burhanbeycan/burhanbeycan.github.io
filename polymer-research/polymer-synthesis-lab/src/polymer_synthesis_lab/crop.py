"""Recipe-only PEtOx prediction and a grouped audit against a single Mn assay."""
import csv
from dataclasses import dataclass
from importlib.resources import files
import io
import math
from statistics import median

import numpy as np

from .model import finite

# Table 2, DOI 10.1021/acs.macromol.2c01930 (MeOTs, EtOx, MeCN).
A_L_MOL_S = 1.06e8
EA_J_MOL = 71600.0
R_J_MOL_K = 8.31446261815324
REPEAT_G_MOL = 99.133  # C5H9NO, average atomic masses
ENDS = {"Me/OH": 32.042, "Me/N3": 57.056}
HEADER = ("paper", "series", "sample", "reactor", "assay", "end_groups", "M0_M", "I0_M",
          "T_C", "time_min", "Mn_g_mol", "X_observed")
SPLITS = ("paper", "series", "concentration", "feed", "concentration_and_feed")
COEFFICIENT_SOURCES = ("10.1021/acs.macromol.2c01930", "10.1021/acsmacrolett.5b00392")


def load_cohort(name):
    if name not in {"batch_azide", "flow_feed"}:
        raise ValueError("Unknown literature cohort.")
    return files(__package__).joinpath("data", name + ".csv").read_text(encoding="utf-8")


@dataclass(frozen=True)
class Recipe:
    M: float
    I: float
    T: float
    minutes: float
    ends: str

    def __post_init__(self):
        for key, lower, upper in (("M", 1e-6, 12), ("I", 1e-8, 12),
                                  ("T", 0, 200), ("minutes", 0, 100000)):
            object.__setattr__(self, key, finite(getattr(self, key), key, lower, upper))
        if self.ends not in ENDS:
            raise ValueError("Choose Me/OH or Me/N3 end groups.")

    @property
    def feed(self):
        return self.M / self.I

    @property
    def kp(self):
        return A_L_MOL_S * math.exp(-EA_J_MOL / (R_J_MOL_K * (273.15 + self.T)))

    @property
    def exposure(self):
        return self.kp * self.I * 60 * self.minutes

    def predict(self):
        x = -math.expm1(-self.exposure)
        residue = self.feed * REPEAT_G_MOL
        return x, ENDS[self.ends] + x * residue, ENDS[self.ends] + residue


@dataclass(frozen=True)
class Observation:
    paper: str
    series: str
    sample: str
    reactor: str
    assay: str
    recipe: Recipe
    Mn: float
    conversion: float | None


def parse_cohort(raw):
    if not isinstance(raw, str) or len(raw) > 150000:
        raise ValueError("Submit CSV text up to 150000 characters.")
    reader = csv.DictReader(io.StringIO(raw.strip()), strict=True)
    try:
        if not reader.fieldnames or sorted(reader.fieldnames) != sorted(HEADER):
            raise ValueError("CSV needs exactly these named columns: " + ",".join(HEADER))
        rows = []
        identifiers = set()
        for lineno, fields in enumerate(reader, 2):
            try:
                if None in fields or any(v is None for v in fields.values()):
                    raise ValueError("unexpected number of fields")
                fields = {k: v.strip() for k, v in fields.items()}
                if any(not fields[k] or len(fields[k]) > 160 for k in HEADER[:6]):
                    raise ValueError("missing or overlong sample, method or provenance identifier")
                paper = fields["paper"].lower().removeprefix("https://doi.org/")
                if any(source in paper for source in COEFFICIENT_SOURCES):
                    raise ValueError("kinetic coefficient sources cannot be counted as independent test papers")
                key = paper, fields["series"], fields["sample"]
                if key in identifiers:
                    raise ValueError("duplicate paper/series/sample identifier")
                identifiers.add(key)
                recipe = Recipe(fields["M0_M"], fields["I0_M"], fields["T_C"],
                                fields["time_min"], fields["end_groups"])
                if recipe.minutes == 0:
                    raise ValueError("a scored polymer sample requires a positive reaction time")
                m = finite(fields["Mn_g_mol"], "Measured Mn", 1e-6, 1e9)
                observed = fields["X_observed"]
                x = finite(observed, "Conversion", 0, 1) if observed else None
                rows.append(Observation(*key, fields["reactor"], fields["assay"], recipe, m, x))
            except ValueError as exc:
                raise ValueError(f"CSV row {lineno}: {exc}") from exc
            if len(rows) > 500:
                raise ValueError("Maximum 500 observations per audit.")
    except csv.Error as exc:
        raise ValueError("Invalid CSV: " + str(exc)) from exc
    if not rows:
        raise ValueError("Provide at least one measured Mn.")
    if len({(r.reactor, r.assay, r.recipe.ends) for r in rows}) != 1:
        raise ValueError("Audit one reactor type, comparable Mn assay and terminal chemistry per cohort.")
    return rows


def group_key(row, split):
    # Six significant figures tolerate round-off in I0 back-calculated from a nominal feed.
    m, feed = (format(v, ".6g") for v in (row.recipe.M, row.recipe.feed))
    return {"paper": row.paper, "series": (row.paper, row.series), "concentration": m,
            "feed": feed, "concentration_and_feed": (m, feed)}[split]


def audit(raw, split):
    """Published kinetic model is fixed. A secondary assay-scale model learns only from training Mn."""
    if split not in SPLITS:
        raise ValueError("Unknown grouping for validation.")
    rows = parse_cohort(raw)
    groups = {}
    for row in rows:
        groups.setdefault(group_key(row, split), []).append(row)
    results, folds = [], []
    for key, test in groups.items():
        train = [r for r in rows if group_key(r, split) != key]
        if split == "concentration_and_feed":
            train = [r for r in train if all(a != b for a, b in zip(group_key(r, split), key))]
        # Keep each reaction series at equal weight, even if sampled at many time points.
        ratios = {}
        for r in train:
            ratios.setdefault((r.paper, r.series), []).append(r.Mn / r.recipe.predict()[1])
        scale = median([median(values) for values in ratios.values()]) if ratios else None
        folds.append([str(key), len(train), len(test), scale])
        for r in test:
            x, mn, full = r.recipe.predict()
            calibrated = scale * mn if scale is not None else None
            def ape(prediction):
                return 100 * abs(prediction - r.Mn) / r.Mn if prediction is not None else None
            results.append({"sample": r.sample, "paper": r.paper, "series": r.series, "X": x,
                            "X_observed": r.conversion, "Mn": r.Mn, "kinetic": mn, "full": full,
                            "scaled": calibrated, "ape_kinetic": ape(mn), "ape_full": ape(full),
                            "ape_scaled": ape(calibrated)})
    summaries = {}
    for key in ("kinetic", "full", "scaled"):
        errors = [r["ape_" + key] for r in results if r["ape_" + key] is not None]
        summaries[key] = (len(errors), median(errors), float(np.percentile(errors, 90))) if errors else (0, None, None)
    return rows, folds, results, summaries
