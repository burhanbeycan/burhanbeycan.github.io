"""Stationary two-state sequence ensemble for partially hydrolysed PVAc."""
from collections import Counter
import math
import random


def integer(value, label, minimum, maximum):
    if (isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value)
            or int(value) != value or not minimum <= value <= maximum):
        raise ValueError(f"{label} must be an integer from {minimum} to {maximum}.")
    return int(value)


def fraction(value, label):
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f"{label} must be a finite mole fraction from 0 to 1.")
    return value


def pva_ensemble(mode, oh_fraction, dyads, dp, n_chains, seed):
    dp = integer(dp, "Model sequence length", 2, 64)
    n_chains = integer(n_chains, "Number of chains", 1, 4096)
    seed = integer(seed, "Random seed", 0, 2**32 - 1)
    if mode == "pva_random":
        h = fraction(oh_fraction, "Hydrolysed repeat fraction")
        mix = 2 * h * (1 - h)
        target = [h*h, mix, (1-h)**2]
    elif mode == "pva_nmr":
        target = [fraction(v, "NMR dyad fraction") for v in dyads]
        if len(target) != 3 or not math.isclose(sum(target), 1, rel_tol=0, abs_tol=1e-6):
            raise ValueError("The three normalized NMR dyad fractions must sum to one (within 0.000001).")
        target = [v/sum(target) for v in target]
        mix = target[1]
        h = target[0] + mix/2
    else:
        raise ValueError("Unsupported PVA ensemble mode.")
    # Stationarity requires equal OH→OAc and OAc→OH flux. Each conditional
    # probability is bounded by construction even for a pure homopolymer.
    acetate_fraction = target[2] + mix / 2
    switch_oh = mix/(2*h) if h else 0
    switch_oac = mix/(2*acetate_fraction) if acetate_fraction else 0
    rng = random.Random(seed)
    pairs = Counter()
    hydrolysis_hist = Counter()
    samples = []
    for idx in range(n_chains):
        state = "A" if rng.random() < h else "B"  # A = vinyl alcohol
        sequence = [state]
        for _ in range(1, dp):
            previous = state
            if rng.random() < (switch_oh if state == "A" else switch_oac):
                state = "B" if state == "A" else "A"
            pairs["".join(sorted((previous, state)))] += 1
            sequence.append(state)
        hydrolysis_hist[sequence.count("A")] += 1
        if idx < 5:
            samples.append("".join(sequence))
    return {"target": target, "oh_fraction": h,
            "observed_dyads": [pairs[k] / ((dp-1)*n_chains) for k in ("AA", "AB", "BB")],
            "sample_oh_fraction": sum(k*v for k, v in hydrolysis_hist.items()) / (dp*n_chains),
            "histogram": sorted(hydrolysis_hist.items()),
            "switch_oh_to_oac": switch_oh if h else None,
            "switch_oac_to_oh": switch_oac if acetate_fraction else None,
            "eta": mix/(2*h*acetate_fraction) if h and acetate_fraction else None,
            "mean_run_oh": 1/switch_oh if switch_oh else None,
            "mean_run_oac": 1/switch_oac if switch_oac else None,
            "examples": samples}
