"""Analytical chemical balances. These functions do not predict conversion or kinetics."""
from __future__ import annotations
import math


def finite(value, name, minimum=None, maximum=None):
    if isinstance(value, (bool, complex)):
        raise ValueError(f"{name} must be a finite scalar.")
    try:
        v = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite scalar.") from exc
    if not math.isfinite(v) or minimum is not None and v < minimum or maximum is not None and v > maximum:
        raise ValueError(f"{name} is outside its finite permitted range.")
    return v


def step_growth(n_aa_mol, n_bb_mol, conversion, mass_aa_g_mol, mass_bb_g_mol, byproduct_g_mol=0):
    """Ideal AA+BB acyclic population, including unreacted monomers.

    p refers to the initially limiting functional groups. Each intermolecular
    bond consumes one A and one B and removes at most one specified byproduct.
    DPn counts original AA/BB molecules, not paired chemical repeat units.
    """
    a = finite(n_aa_mol, "AA amount", 1e-12, 1e6)
    b = finite(n_bb_mol, "BB amount", 1e-12, 1e6)
    p = finite(conversion, "Limiting-group conversion", 0, 1)
    ma = finite(mass_aa_g_mol, "AA molar mass", 1e-6, 1e7)
    mb = finite(mass_bb_g_mol, "BB molar mass", 1e-6, 1e7)
    loss = finite(byproduct_g_mol, "Byproduct molar mass per bond", 0, 1e7)
    # At most two condensation losses per AA+BB pair. The remaining pair
    # residue must have strictly positive mass, even at full conversion.
    if 2 * loss >= ma + mb:
        raise ValueError("The specified byproduct removes all or more of the AA+BB pair mass.")
    ratio = min(a, b) / max(a, b)
    bonds = 2 * min(a, b) * p
    # This form preserves the small nonzero molecule count near completion.
    molecules = abs(a - b) + 2 * min(a, b) * (1 - p)
    initial_mass = a * ma + b * mb
    byproduct_mass = bonds * loss
    retained_mass = initial_mass - byproduct_mass
    return {
        "r": ratio, "p_limiting": p, "limiting_group": "A" if a < b else "B" if b < a else "balanced",
        "initial_molecules_mol": a + b, "bonds_mol": bonds,
        "remaining_molecules_mol": molecules,
        "remaining_A_mol": 2 * a - bonds, "remaining_B_mol": 2 * b - bonds,
        "conversion_A": bonds / (2 * a), "conversion_B": bonds / (2 * b),
        "initial_mass_g": initial_mass, "byproduct_mass_g": byproduct_mass,
        "retained_population_mass_g": retained_mass,
        "DPn_original_monomer_units": (a + b) / molecules if molecules else None,
        "Mn_population_g_mol": retained_mass / molecules if molecules else None,
        "balanced_Flory_length_dispersity": 1 + p if a == b else None,
    }


def living_chain_growth(monomer_mol, initiator_mol, conversion, efficiency, chains_per_initiator,
                        repeat_mass_g_mol, end_mass_g_mol):
    """Ideal number-average residue count and Mn conditional on supplied p and f.

    All initiated chains are counted, including zero-addition chains. No transfer,
    termination, combination, initiator mass inference or kinetic fitting.
    """
    mon = finite(monomer_mol, "Initial monomer", 1e-12, 1e6)
    ini = finite(initiator_mol, "Initial initiator", 1e-12, 1e6)
    p = finite(conversion, "Monomer conversion", 0, 1)
    eff = finite(efficiency, "Initiation efficiency", 1e-9, 1)
    sites = finite(chains_per_initiator, "Chains per initiator", 1, 100)
    if sites != int(sites):
        raise ValueError("Chains per initiator must be a positive integer.")
    m0 = finite(repeat_mass_g_mol, "Incorporated repeat mass", 1e-6, 1e7)
    mend = finite(end_mass_g_mol, "Combined end-group mass", 0, 1e7)
    chains = ini * eff * sites
    reacted = mon * p
    dp = reacted / chains
    mass = reacted * m0 + chains * mend
    return {"initiated_chains_mol": chains, "reacted_monomer_mol": reacted,
            "unreacted_monomer_mol": mon * (1 - p), "DPn": dp,
            "Mn_g_mol": dp * m0 + mend, "chain_population_mass_g": mass,
            "end_mass_fraction": (chains * mend / mass) if mass else None}


def terminal_copolymer(f_a, r_a, r_b):
    """Stationary terminal-model composition at one fixed monomer feed.

    Strict interior feeds avoid undefined transition probabilities for absent
    monomers. Ordered dyads come from stationary chain-end transitions.
    """
    f = finite(f_a, "Monomer A feed fraction", 1e-9, 1 - 1e-9)
    ra = finite(r_a, "A reactivity ratio", 0, 1e6)
    rb = finite(r_b, "B reactivity ratio", 0, 1e6)
    g = 1 - f
    paa = ra * f / (ra * f + g)
    pbb = rb * g / (f + rb * g)
    pab = g / (ra * f + g)
    pba = f / (f + rb * g)
    denominator = ra * f * f + 2 * f * g + rb * g * g
    fa = (ra * f * f + f * g) / denominator
    fb = 1 - fa
    return {"F_A": fa, "F_B": fb, "P_AA": paa, "P_AB": pab, "P_BA": pba, "P_BB": pbb,
            "dyad_AA": fa * paa, "dyad_AB": fa * pab, "dyad_BA": fb * pba, "dyad_BB": fb * pbb,
            "mean_A_run": 1 / pab, "mean_B_run": 1 / pba}
