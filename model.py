"""History-dependent environmental surface colonization, version 1.0.

Dimensionless, species-agnostic stochastic model of attachment, detachment and
division on a patchy surface. Six generative hypotheses share one simulator:

    passive    heterogeneous attachment and growth, no transfer effect
    carryover  passive dynamics plus residual cells left on a visited surface
    selection  passive dynamics plus an adhesive subpopulation among exposed cells
    cell       cell-carried history state (acquired on the surface, inherited)
    surface    surface-carried history state (deposited by contacts and residents)
    both       cell and surface history states together

Every hypothesis contains heritable cell-to-cell heterogeneity in adhesiveness
and quenched patch heterogeneity, so the passive null has heterogeneous rates.
Carryover mimics surface history and selection mimics cell history without any
history state; they are the confounders a transfer experiment must reject.

Parameters are per-replicate arrays, so one call can simulate many experiments
with different parameter draws. No organism is fitted.
"""
from dataclasses import dataclass, asdict, fields, replace
import numpy as np

VERSION = "1.0.0"
MECHANISMS = ("passive", "carryover", "selection", "cell", "surface", "both")
ARMS = ("fresh_fresh", "exposed_fresh", "fresh_conditioned", "exposed_conditioned")
LABELS = {"passive": "Heterogeneous passive", "carryover": "Residual-cell carryover",
          "selection": "Adhesive selection", "cell": "Cell history",
          "surface": "Surface history", "both": "Cell + surface history"}
SHORT = {"passive": "Passive", "carryover": "Carryover", "selection": "Selection",
         "cell": "Cell", "surface": "Surface", "both": "Both"}
ARM_LABELS = {"fresh_fresh": "Fresh cells / fresh surface",
              "exposed_fresh": "Exposed cells / fresh surface",
              "fresh_conditioned": "Fresh cells / visited surface",
              "exposed_conditioned": "Exposed cells / visited surface"}
# Couplings switched on by each hypothesis. Everything else is shared.
ACTIVE = {"passive": (), "carryover": ("carryover",), "selection": ("selection",),
          "cell": ("cell",), "surface": ("surface",), "both": ("cell", "surface")}
OBS_TIMES = tuple(float(t) for t in range(1, 19))


@dataclass(frozen=True)
class Config:
    # Geometry, dose and observation window
    patches: int = 6
    capacity_per_patch: int = 12
    arrivals: int = 40
    arrival_duration: float = 4.0
    horizon: float = 18.0
    dt: float = 0.05
    # Shared passive dynamics
    attachment_probability: float = 0.28
    division_rate: float = 0.70
    detachment_rate: float = 0.85
    patch_heterogeneity: float = 0.35
    cell_heterogeneity: float = 0.40
    # Cell history
    memory_acquisition_rate: float = 1.0
    cell_memory_decay_time: float = 4.0
    cell_attachment_effect: float = 1.8
    cell_retention_effect: float = 1.6
    # Surface history
    surface_decay_time: float = 12.0
    surface_attachment_effect: float = 2.2
    surface_retention_effect: float = 2.0
    contact_trace: float = 0.08
    resident_deposition_rate: float = 0.15
    # Confounders
    carryover_cells: float = 2.0          # mean residual cells per visited surface
    storage_loss_rate: float = 0.05       # loss of residual cells during the gap
    selection_shift: float = 1.0          # adhesiveness shift of exposed cells (SD units of logit)
    selection_relaxation_time: float = 40.0
    # Preparation
    conditioning_visits: int = 4
    conditioning_residence: float = 0.15
    conditioning_interval: float = 0.35
    # Operational endpoint
    establishment_threshold: int = 12
    establishment_persistence: float = 1.0

    def __post_init__(self):
        if self.patches < 1 or self.capacity_per_patch < 1 or self.arrivals < 0:
            raise ValueError("Invalid count")
        if not 0 < self.attachment_probability < 1:
            raise ValueError("attachment_probability must lie in (0,1)")
        if min(self.dt, self.horizon, self.cell_memory_decay_time,
               self.surface_decay_time, self.selection_relaxation_time) <= 0:
            raise ValueError("Times must be positive")
        if self.arrival_duration <= 0 or self.arrival_duration > self.horizon:
            raise ValueError("Invalid arrival duration")
        if self.establishment_threshold > self.patches * self.capacity_per_patch:
            raise ValueError("Threshold exceeds total capacity")
        for f in fields(self):
            if f.name in ("patches", "capacity_per_patch", "arrivals", "horizon", "dt",
                          "attachment_probability", "arrival_duration",
                          "establishment_threshold"):
                continue
            if getattr(self, f.name) < 0:
                raise ValueError(f.name)


# Parameters that may vary between replicates inside one batch.
VARYING = ("attachment_probability", "division_rate", "detachment_rate",
           "patch_heterogeneity", "cell_heterogeneity", "memory_acquisition_rate",
           "cell_memory_decay_time", "cell_attachment_effect", "cell_retention_effect",
           "surface_decay_time", "surface_attachment_effect", "surface_retention_effect",
           "contact_trace", "resident_deposition_rate", "carryover_cells",
           "storage_loss_rate", "selection_shift", "selection_relaxation_time")


def mechanism_params(cfg, mechanism):
    """Return a parameter dictionary with inactive couplings set to zero."""
    if mechanism not in MECHANISMS:
        raise ValueError(mechanism)
    p = asdict(cfg)
    on = ACTIVE[mechanism]
    if "cell" not in on:
        p["cell_attachment_effect"] = p["cell_retention_effect"] = 0.0
    if "surface" not in on:
        p["surface_attachment_effect"] = p["surface_retention_effect"] = 0.0
    if "carryover" not in on:
        p["carryover_cells"] = 0.0
    if "selection" not in on:
        p["selection_shift"] = 0.0
    return p


def conditioning_history(params):
    """History left by the prescribed transient preparation.

    Each visit ends in complete detachment and no reproduction occurs, so the
    history models leave zero attached cells. Values broadcast over arrays.
    The carryover hypothesis instead leaves residual cells (handled in simulate).
    """
    km = np.asarray(params["memory_acquisition_rate"], float)
    tm = np.asarray(params["cell_memory_decay_time"], float)
    ts = np.asarray(params["surface_decay_time"], float)
    q = np.asarray(params["contact_trace"], float)
    dep = np.asarray(params["resident_deposition_rate"], float)
    visits = int(params["conditioning_visits"])
    res, gap = float(params["conditioning_residence"]), float(params["conditioning_interval"])
    K = int(params["capacity_per_patch"])
    m = np.zeros(np.broadcast(km, tm, ts, q, dep).shape)
    s = np.zeros_like(m)
    duration = 0.0
    for j in range(visits):
        m = 1 - (1 - m) * np.exp(-km * res)
        s = 1 - (1 - s) * np.exp(-q)
        a = dep / K
        r = a + 1 / ts
        eq = a / r
        s = eq + (s - eq) * np.exp(-r * res)
        duration += res
        if j + 1 < visits:
            m = m * np.exp(-gap / tm)
            s = s * np.exp(-gap / ts)
            duration += gap
    return {"cell_memory": m, "surface_trace": s, "attached_cells": 0, "duration": duration}


def _column(params, name, R):
    v = np.asarray(params[name], float)
    if v.ndim == 0:
        v = np.full(R, float(v))
    if v.shape != (R,):
        raise ValueError(f"{name} must be scalar or length {R}")
    return v


def simulate_batch(cfg, params, exposed, conditioned, gap, seed=1234,
                   obs_times=OBS_TIMES, record=False):
    """Simulate R independent challenge experiments in one vectorized pass.

    cfg       Config holding structural settings (geometry, dose, step, endpoint).
    params    dict of per-replicate arrays (length R) or scalars for VARYING keys.
    exposed   boolean array: incoming cells previously made transient contacts.
    conditioned boolean array: surface previously received transient contacts.
    gap       array of times between preparation and challenge.

    Cells arrive on a fixed schedule, attach with a logistic probability, divide
    into vacant slots, and detach with a hazard reduced by patch quality, cell
    adhesiveness, and the two history states. Newborns inherit adhesiveness,
    history state and founder identity. Detached cells leave the domain.
    """
    exposed = np.asarray(exposed, bool)
    R = exposed.size
    conditioned = np.broadcast_to(np.asarray(conditioned, bool), (R,))
    gap = np.broadcast_to(np.asarray(gap, float), (R,))
    if np.any(gap < 0):
        raise ValueError("gap must be non-negative")
    P, K = cfg.patches, cfg.capacity_per_patch
    pr = {k: _column(params, k, R) for k in VARYING}
    rng = np.random.default_rng(seed)
    rows = np.arange(R)
    shape = (R, P, K)

    prep = conditioning_history({**asdict(cfg), **pr})
    m0 = np.where(exposed, prep["cell_memory"] * np.exp(-gap / pr["cell_memory_decay_time"]), 0.)
    s0 = np.where(conditioned, prep["surface_trace"] * np.exp(-gap / pr["surface_decay_time"]), 0.)
    sel = np.where(exposed, pr["selection_shift"] * np.exp(-gap / pr["selection_relaxation_time"]), 0.)

    patch_effect = rng.normal(0, 1, (R, P)) * pr["patch_heterogeneity"][:, None]
    arrival_patch = rng.integers(P, size=(R, cfg.arrivals))
    arrival_u = rng.random((R, cfg.arrivals))
    vacancy_u = rng.random((R, cfg.arrivals))
    arrival_times = (np.arange(cfg.arrivals) + .5) * cfg.arrival_duration / max(cfg.arrivals, 1)
    # Selected cells relax toward the population distribution while in suspension,
    # during the gap and during the challenge, exactly as cell history decays.
    arrival_z = (rng.normal(0, 1, (R, cfg.arrivals)) * pr["cell_heterogeneity"][:, None]
                 + sel[:, None] * np.exp(-arrival_times[None, :] / pr["selection_relaxation_time"][:, None]))
    carry_u = rng.random((R, P))
    carry_z = rng.normal(0, 1, (R, P, K)) * pr["cell_heterogeneity"][:, None, None]

    occupied = np.zeros(shape, bool)
    memory = np.zeros(shape, np.float64)
    adhes = np.zeros(shape, np.float64)
    founder = np.zeros(shape, np.int32)       # 0 = none, -1 = residual carryover, j+1 = arrival j
    generation = np.zeros(shape, np.int16)
    born = np.zeros(shape, np.float32)        # attachment or birth time
    surface = np.repeat(s0[:, None], P, axis=1)

    # Residual cells: Poisson number per patch, thinned by storage loss over the gap.
    lam = (pr["carryover_cells"] / P)[:, None] * np.exp(
        -pr["storage_loss_rate"][:, None] * np.exp(-patch_effect) * gap[:, None])
    lam = np.where(conditioned[:, None], lam, 0.)
    # Inverse-CDF Poisson draw keeps random streams aligned across arms.
    pmf = np.exp(-lam)
    cdf = pmf.copy()
    n_carry = np.zeros((R, P), np.int32)
    for k in range(1, K + 1):
        n_carry += (carry_u >= cdf)
        pmf = pmf * lam / k
        cdf = cdf + pmf
    place = np.arange(K)[None, None, :] < n_carry[:, :, None]
    occupied |= place
    adhes[place] = carry_z[place]
    founder[place] = -1
    initial_residual = place.sum(axis=(1, 2)).astype(np.int32)

    arrival_steps = np.floor(arrival_times / cfg.dt).astype(int)
    by_step = {}
    for j, st in enumerate(arrival_steps):
        by_step.setdefault(int(st), []).append(j)
    half = cfg.arrivals // 2

    z = lambda n: np.zeros(R, np.int32)
    arrived, rejected, detached, divisions = z(0), z(0), z(0), z(0)
    attached_early, attached_late = z(0), z(0)
    exposed_attach_m = np.zeros(R)
    founder_departures, founder_quick = z(0), z(0)
    desc_departures, desc_quick = z(0), z(0)
    residual_departures = z(0)
    streak = np.zeros(R)
    established_time = np.full(R, np.nan)
    integrated = np.zeros(R)
    obs_steps = {int(round(t / cfg.dt)) - 1: i for i, t in enumerate(obs_times)}
    counts_obs = np.zeros((R, len(obs_times)), np.int16)
    residual_obs = np.zeros((R, len(obs_times)), np.int16)
    trace = []

    base_logit = np.log(pr["attachment_probability"] / (1 - pr["attachment_probability"]))
    b = pr["division_rate"][:, None]
    d0 = pr["detachment_rate"][:, None, None]
    gm = pr["cell_retention_effect"][:, None, None]
    gs = pr["surface_retention_effect"][:, None, None]
    km = pr["memory_acquisition_rate"][:, None, None]
    tm = pr["cell_memory_decay_time"]
    bm, bs = pr["cell_attachment_effect"], pr["surface_attachment_effect"]
    dep = pr["resident_deposition_rate"][:, None]
    inv_ts = 1 / pr["surface_decay_time"][:, None]
    q = pr["contact_trace"]
    cell_on = bool(np.any(bm > 0) or np.any(gm > 0))
    steps = int(np.ceil(cfg.horizon / cfg.dt - 1e-9))

    for step in range(steps):
        t = step * cfg.dt
        dt = min(cfg.dt, cfg.horizon - t)
        old_count = occupied.sum(axis=2)
        a = dep * old_count / K
        rate = a + inv_ts
        eq = a / rate
        surface = eq + (surface - eq) * np.exp(-rate * dt)

        hazard = d0 * np.exp(-patch_effect[:, :, None] - adhes - gm * memory - gs * surface[:, :, None])
        departures = occupied & (rng.random(shape) < -np.expm1(-hazard * dt))
        if departures.any():
            r_i, p_i, k_i = np.nonzero(departures)
            res = (t + dt) - born[r_i, p_i, k_i]
            quick = res < 1.0
            g0 = generation[r_i, p_i, k_i] == 0
            fid = founder[r_i, p_i, k_i]
            arr = g0 & (fid > 0)
            founder_departures += np.bincount(r_i[arr], minlength=R).astype(np.int32)
            founder_quick += np.bincount(r_i[arr & quick], minlength=R).astype(np.int32)
            dsc = ~g0 & (fid > 0)
            desc_departures += np.bincount(r_i[dsc], minlength=R).astype(np.int32)
            desc_quick += np.bincount(r_i[dsc & quick], minlength=R).astype(np.int32)
            residual_departures += np.bincount(r_i[fid < 0], minlength=R).astype(np.int32)
            detached += np.bincount(r_i, minlength=R).astype(np.int32)
            occupied[departures] = False
            memory[departures] = 0
            adhes[departures] = 0
            founder[departures] = 0
            generation[departures] = 0
        if cell_on:
            memory = np.where(occupied, 1 - (1 - memory) * np.exp(-km * dt), memory)

        # Each vacant slot fills at rate b N / K, giving logistic growth b N (1 - N/K).
        live = occupied.sum(axis=2)
        birth_p = -np.expm1(-b * live / K * dt)
        births = (~occupied) & (rng.random(shape) < birth_p[:, :, None])
        if births.any():
            r_i, p_i, k_i = np.nonzero(births)
            occ = occupied[r_i, p_i]                      # parents present before births
            rank = (rng.random(r_i.size) * live[r_i, p_i]).astype(int)
            parent = np.argmax(np.cumsum(occ, axis=1) > rank[:, None], axis=1)
            memory[r_i, p_i, k_i] = memory[r_i, p_i, parent]
            adhes[r_i, p_i, k_i] = adhes[r_i, p_i, parent]
            founder[r_i, p_i, k_i] = founder[r_i, p_i, parent]
            generation[r_i, p_i, k_i] = generation[r_i, p_i, parent] + 1
            born[r_i, p_i, k_i] = t + dt
            occupied[r_i, p_i, k_i] = True
            divisions += np.bincount(r_i, minlength=R).astype(np.int32)

        for j in by_step.get(step, []):
            p = arrival_patch[:, j]
            arrived += 1
            surface[rows, p] = 1 - (1 - surface[rows, p]) * np.exp(-q)
            incoming_m = m0 * np.exp(-t / tm)
            logit = (base_logit + patch_effect[rows, p] + arrival_z[:, j]
                     + bm * incoming_m + bs * surface[rows, p])
            prob = 1 / (1 + np.exp(-logit))
            vac = ~occupied[rows, p]
            nfree = vac.sum(axis=1)
            choose = (vacancy_u[:, j] * np.maximum(nfree, 1)).astype(int)
            where = np.argmax(np.cumsum(vac, axis=1) > choose[:, None], axis=1)
            accept = (arrival_u[:, j] < prob) & (nfree > 0)
            rr, pp, ss = rows[accept], p[accept], where[accept]
            occupied[rr, pp, ss] = True
            memory[rr, pp, ss] = incoming_m[accept]
            adhes[rr, pp, ss] = arrival_z[accept, j]
            founder[rr, pp, ss] = j + 1
            generation[rr, pp, ss] = 0
            born[rr, pp, ss] = t + dt
            rejected += (~accept).astype(np.int32)
            if j < half:
                attached_early += accept.astype(np.int32)
            else:
                attached_late += accept.astype(np.int32)

        count = occupied.sum(axis=(1, 2))
        balance = initial_residual + arrived + divisions - detached - rejected - count
        if np.any(balance != 0):
            raise AssertionError("Cell accounting failed")
        streak = np.where(count >= cfg.establishment_threshold, streak + dt, 0)
        newly = (np.isnan(established_time) & (count >= cfg.establishment_threshold)
                 & (streak >= cfg.establishment_persistence - 1e-10))
        established_time[newly] = t + dt
        integrated += count * dt
        if step in obs_steps:
            i = obs_steps[step]
            counts_obs[:, i] = count
            residual_obs[:, i] = (occupied & (founder < 0)).sum(axis=(1, 2))
        if record:
            trace.append({"time": t + dt, "mean_attached": float(count.mean()),
                          "mean_surface_trace": float(surface.mean()),
                          "established_fraction": float(np.isfinite(established_time).mean())})

    if not (np.all(memory >= -1e-12) and np.all(memory <= 1 + 1e-12)
            and np.all(surface >= -1e-12) and np.all(surface <= 1 + 1e-12)):
        raise AssertionError("History state outside [0,1]")
    final = occupied.sum(axis=(1, 2))
    founders_ret = (occupied & (generation == 0) & (founder > 0)).sum(axis=(1, 2))
    residual_ret = (occupied & (founder < 0)).sum(axis=(1, 2))
    lineages = np.array([np.unique(founder[r][occupied[r]]).size for r in range(R)])
    return {"established": np.isfinite(established_time), "establishment_time": established_time,
            "final_attached": final, "arrivals": arrived, "rejected": rejected,
            "departures": detached, "divisions": divisions,
            "retained_founders": founders_ret,
            "retained_descendants": final - founders_ret - residual_ret,
            "retained_residual": residual_ret, "initial_residual": initial_residual,
            "residual_departures": residual_departures,
            "retained_lineages": lineages, "integrated_occupancy": integrated,
            "attached_early": attached_early, "attached_late": attached_late,
            "founder_departures": founder_departures, "founder_quick": founder_quick,
            "descendant_departures": desc_departures, "descendant_quick": desc_quick,
            "counts": counts_obs, "residual_counts": residual_obs, "trace": trace,
            "initial_memory": m0, "initial_surface_trace": s0, "initial_selection": sel}


def simulate(cfg, mechanism, arm, gap=0.0, n=200, seed=1234, record=False, overrides=None):
    """Simulate n replicate surfaces of one mechanism, arm and gap."""
    if mechanism not in MECHANISMS or arm not in ARMS or gap < 0 or n < 1:
        raise ValueError("Invalid mechanism, arm, gap, or ensemble size")
    params = mechanism_params(cfg, mechanism)
    if overrides:
        params.update(overrides)
        # Overrides may rescale active couplings but never switch on inactive ones.
        off = mechanism_params(replace(cfg, **{k: v for k, v in overrides.items()
                                               if k in asdict(cfg) and np.ndim(v) == 0}), mechanism)
        for k in ("cell_attachment_effect", "cell_retention_effect", "surface_attachment_effect",
                  "surface_retention_effect", "carryover_cells", "selection_shift"):
            if off[k] == 0:
                params[k] = 0.0
    exposed = np.full(n, arm in ("exposed_fresh", "exposed_conditioned"))
    conditioned = np.full(n, arm in ("fresh_conditioned", "exposed_conditioned"))
    return simulate_batch(cfg, params, exposed, conditioned, np.full(n, float(gap)),
                          seed=seed, record=record)
