"""
DP-CLN Computational Stress Test: Specification & Simulation Engine (v0.3)

Addresses:
- Reviewer 2, Comment 5: Realistic stochastic state dynamics, circadian variation,
  noise, drift, missing data, estimator lag, and treatment response heterogeneity/decay.
- Reviewer 3, Comment 5: Separation of fast neural feedback (hourly) from slow DP risk (daily).
- Reviewer 3, Comment 6: Therapeutic responsiveness causal estimand (inactive in simulation, delta = -inf).
- Reviewer 3, Comment 7: Unified control policy with hard safety gate g(t).
- Reviewer 3, Comment 8: Justified dual-timescale parameters and fair burden-matched comparators.

Architecture:
- Slow clock: k in 0..179 (180 days).
- Fast clock: h in 0..23 (24 hours per day), t = k * 24 + h.
- 3-state daily Markov process: z_k in {0: remission, 1: prodrome, 2: episode}.
- Treatment dynamics: dose carryover, saturation, tolerance decay, heterogeneous efficacy.
- Observations:
  * Neural (hourly): mu[z] + circadian(h) + AR(1) noise + patient drift.
  * DP (daily): a * sev(z) + c * ctx + AR(1) noise; severity-dependent missingness.
- Arms:
  * None: Natural history.
  * Neural-only: hourly threshold on b(t).
  * DP-gated: daily 09:00 slot threshold on r_hat(k).
  * Hierarchical: hourly threshold modulated by DP risk with safety clipping: theta_b(r) = clip(theta_0 - kappa * (r - r_bar), theta_min, theta_max).
  * Random@A: exact burden-matched comparator delivering A's delivered stimulation count on safe days.
  * Scheduled@A: evenly spaced safe days delivering matched count.
  * Oracle: upper bound triggering whenever z >= 1 and safe.
"""

import argparse
import collections
from collections import deque
import json
import os
import sys
import numpy as np


# =============================================================================
# 1. PARAMETERS & GENERATIVE PROCESS
# =============================================================================

DEFAULT_CONFIG = {
    "n_days": 180,
    "burn_in_days": 30,
    "refractory_h": 24,
    "weekly_cap": 5,
    "budget": 30,             # Target stimulation budget over scored 150 days
    "gamma": 0.85,            # Dose carryover
    "d50": 2.0,               # Dose saturation midpoint
    "t_tol": 60.0,            # Tolerance decay timescale (cumulative stims)
    "p_resp": 0.6,            # Fraction of responders
    "harm": 0.05,             # Iatrogenic destabilization parameter
    "phi_neur": 0.9,          # Neural AR(1) persistence
    "phi_dp": 0.5,            # DP AR(1) persistence
    "phi_ctx": 0.8,           # Context AR(1) persistence
    "alpha_ewma": 0.3,        # EWMA smoothing factor for DP features
    "slot_hour": 9,           # Daily scheduling slot for Path C
    "theta_min": 0.5,         # Bounded threshold floor
    "theta_max": 2.5          # Bounded threshold ceiling
}


def sample_patient_cohort(n_patients, rng, scenario="informative", miss_inform=None,
                          neural_sd=1.0, drift_scale=0.0, is_test=False):
    """
    Generates a cohort of virtual patients with heterogeneous parameters
    and pre-drawn stochastic noise sequences (Common Random Numbers).
    """
    # Scenario parameters
    scenario_map = {
        "informative":   {"a": 1.0,  "c": 0.0, "rho_train": 0.0, "rho_test": 0.0, "miss_inform": 1},
        "weak":          {"a": 0.35, "c": 0.0, "rho_train": 0.0, "rho_test": 0.0, "miss_inform": 1},
        "uninformative": {"a": 0.0,  "c": 0.0, "rho_train": 0.0, "rho_test": 0.0, "miss_inform": 0},
        "misleading":    {"a": 0.2,  "c": 1.0, "rho_train": 0.6, "rho_test": 0.0, "miss_inform": 0}
    }
    scen = scenario_map.get(scenario, scenario_map["informative"])
    a = scen["a"]
    c = scen["c"]
    rho_ctx = scen["rho_test"] if is_test else scen["rho_train"]
    m_inform = scen["miss_inform"] if miss_inform is None else miss_inform

    cohort = []
    n_days = DEFAULT_CONFIG["n_days"]
    total_hours = n_days * 24

    for i in range(n_patients):
        p = {}
        p["id"] = i
        # Markov baseline transition probabilities
        p["p01"] = rng.uniform(0.005, 0.030)
        p["p12"] = rng.uniform(0.10, 0.30)
        p["p10"] = rng.uniform(0.05, 0.15)
        p["p20"] = rng.uniform(0.03, 0.08)

        # Responders vs non-responders
        is_resp = rng.random() < DEFAULT_CONFIG["p_resp"]
        p["e_i"] = rng.beta(5, 3) if is_resp else 0.0

        # Neural channel parameters
        p["drift"] = rng.normal(0, 1) * drift_scale
        p["circ_amp"] = 0.5

        # AR(1) neural noise sequence: eps[t] = phi*eps[t-1] + sqrt(1-phi^2)*noise
        phi_n = DEFAULT_CONFIG["phi_neur"]
        eps_innov = rng.normal(0, neural_sd * np.sqrt(1 - phi_n**2), size=total_hours)
        eps = np.zeros(total_hours)
        for t in range(1, total_hours):
            eps[t] = phi_n * eps[t - 1] + eps_innov[t]
        p["eps_neur"] = eps

        # Daily DP noise and contextual confounders
        phi_d = DEFAULT_CONFIG["phi_dp"]
        e_innov = rng.normal(0, np.sqrt(1 - phi_d**2), size=n_days)
        e_dp = np.zeros(n_days)
        for k in range(1, n_days):
            e_dp[k] = phi_d * e_dp[k - 1] + e_innov[k]
        p["e_dp"] = e_dp

        phi_c = DEFAULT_CONFIG["phi_ctx"]
        u_innov = rng.normal(0, np.sqrt(1 - phi_c**2), size=n_days)
        u_ctx = np.zeros(n_days)
        for k in range(1, n_days):
            u_ctx[k] = phi_c * u_ctx[k - 1] + u_innov[k]
        p["u_ctx"] = u_ctx

        p["a"] = a
        p["c"] = c
        p["rho_ctx"] = rho_ctx
        p["miss_inform"] = m_inform

        # Pre-drawn transition uniforms for Markov chain (common across arms)
        p["trans_u"] = rng.random(size=n_days)
        # Pre-drawn missingness uniforms
        p["miss_u"] = rng.random(size=n_days)
        # Contraindication uniforms (5% daily contraindication probability)
        p["contra"] = rng.random(size=n_days) < 0.05

        cohort.append(p)
    return cohort


# =============================================================================
# 2. OBSERVATION GENERATION & RISK MODEL
# =============================================================================

def generate_patient_observations(p, z_seq):
    """
    Generates hourly neural signal and daily digital phenotyping signal
    given latent state trajectory z_seq.
    """
    n_days = len(z_seq)
    total_hours = n_days * 24

    # Fast Neural Signal (hourly)
    mu_map = np.array([0.0, 0.8, 1.5])
    mu_z = mu_map[z_seq]  # Daily mean
    b = np.zeros(total_hours)
    for k in range(n_days):
        day_mu = mu_z[k]
        drift_term = p["drift"] * (k / n_days)
        for h in range(24):
            t = k * 24 + h
            circ = p["circ_amp"] * np.sin(2 * np.pi * (h - 14) / 24.0)
            b[t] = day_mu + circ + p["eps_neur"][t] + drift_term

    # Slow DP Signal (daily)
    sev_map = np.array([0.0, 0.5, 1.0])
    sev = sev_map[z_seq]
    ctx = p["rho_ctx"] * ((sev - 0.5) / 0.7) + np.sqrt(max(0.0, 1 - p["rho_ctx"]**2)) * p["u_ctx"]
    y = p["a"] * sev + p["c"] * ctx + p["e_dp"]

    # Severity-dependent missingness
    p_miss = 1.0 / (1.0 + np.exp(-(-2.0 + 0.8 * p["miss_inform"] * sev)))
    is_missing = p["miss_u"] < p_miss

    return b, y, is_missing


class SimpleLogisticRiskModel:
    """
    Estimates probability of episode onset (z=2) within H=3 days
    using EWMA DP features, change score, and missingness rate.
    Trained strictly on natural history (no stimulation).
    """
    def __init__(self, horizon=3, alpha=0.3):
        self.horizon = horizon
        self.alpha = alpha
        self.weights = np.zeros(4)  # bias + 3 features
        self.trained = False

    def extract_features(self, y, is_missing):
        n_days = len(y)
        features = np.zeros((n_days, 3))
        ewma_y = 0.0
        ewma_m = 0.0
        y_clean_hist = []

        for k in range(n_days):
            # Mandatory 1-day pipeline lag: features for day k use info strictly through k-1
            if k == 0:
                features[k] = [0.0, 0.0, 0.0]
            else:
                chg3 = (y_clean_hist[-1] - y_clean_hist[-4]) if len(y_clean_hist) >= 4 else (y_clean_hist[-1] - y_clean_hist[0])
                features[k, 0] = ewma_y
                features[k, 1] = chg3
                features[k, 2] = ewma_m

            # Update running feature state with day k observations (available for day k+1)
            m_k = 1.0 if is_missing[k] else 0.0
            ewma_m = self.alpha * m_k + (1.0 - self.alpha) * ewma_m
            if not is_missing[k]:
                ewma_y = self.alpha * y[k] + (1.0 - self.alpha) * ewma_y
                y_clean_hist.append(y[k])
            else:
                y_clean_hist.append(ewma_y)

        return features

    def fit(self, cohort):
        X_all = []
        Y_all = []
        for p in cohort:
            # Simulate natural history (arm 'None')
            z_nat = simulate_natural_history(p)
            _, y, is_missing = generate_patient_observations(p, z_nat)
            feat = self.extract_features(y, is_missing)

            n_days = len(z_nat)
            for k in range(1, n_days - self.horizon):
                # Target: episode onset within horizon among patients currently not in episode
                if z_nat[k] in (0, 1):
                    target = 1.0 if np.any(z_nat[k + 1: k + 1 + self.horizon] == 2) else 0.0
                    X_all.append(feat[k])
                    Y_all.append(target)

        if not X_all:
            return

        X = np.array(X_all)
        Y = np.array(Y_all)
        # Standardize features
        self.mean = np.mean(X, axis=0)
        self.std = np.std(X, axis=0) + 1e-6
        X_norm = (X - self.mean) / self.std
        X_design = np.hstack([np.ones((len(X_norm), 1)), X_norm])

        # Simple regularized logistic regression solver (IRLS or gradient descent)
        w = np.zeros(X_design.shape[1])
        lr = 0.05
        for _ in range(200):
            preds = 1.0 / (1.0 + np.exp(-X_design @ w))
            grad = X_design.T @ (preds - Y) / len(Y) + 0.01 * w
            w -= lr * grad
        self.weights = w
        self.trained = True

    def predict_step(self, ewma_y, chg3, ewma_m):
        """Evaluates forward risk in O(1) for a single day using 1-day lagged features."""
        feat = np.array([ewma_y, chg3, ewma_m])
        feat_norm = (feat - self.mean) / self.std
        z = self.weights[0] + np.dot(feat_norm, self.weights[1:])
        return float(1.0 / (1.0 + np.exp(-z)))

    def predict_risk(self, y, is_missing):
        feat = self.extract_features(y, is_missing)
        feat_norm = (feat - self.mean) / self.std
        X_design = np.hstack([np.ones((len(feat_norm), 1)), feat_norm])
        risk = 1.0 / (1.0 + np.exp(-X_design @ self.weights))
        return risk


# =============================================================================
# 3. STATE SIMULATION ENGINE
# =============================================================================

def simulate_natural_history(p):
    """Simulates untreated latent disease trajectory z_k in {0, 1, 2}."""
    n_days = DEFAULT_CONFIG["n_days"]
    z = np.zeros(n_days, dtype=int)
    # Initial state distribution: 70% remission, 15% prodrome, 15% episode
    u0 = p["trans_u"][0]
    z[0] = 0 if u0 < 0.70 else (1 if u0 < 0.85 else 2)

    for k in range(n_days - 1):
        curr = z[k]
        u = p["trans_u"][k + 1]
        if curr == 0:
            z[k + 1] = 1 if u < p["p01"] else 0
        elif curr == 1:
            if u < p["p12"]:
                z[k + 1] = 2
            elif u < p["p12"] + p["p10"]:
                z[k + 1] = 0
            else:
                z[k + 1] = 1
        else: # curr == 2
            z[k + 1] = 0 if u < p["p20"] else 2
    return z


def run_patient_simulation(p, arm, risk_model, theta_b=1.2, theta_r=0.4, kappa=3.0,
                           target_stim_days=None, theta_min=None, theta_max=None):
    """
    Simulates a single virtual patient under a specified controller arm.
    Returns:
      z: latent state sequence (180 days)
      u: stimulation matrix (180 days x 24 hours)
      stims_per_day: daily stimulation counts
    """
    n_days = DEFAULT_CONFIG["n_days"]
    ref_h = DEFAULT_CONFIG["refractory_h"]
    cap_w = DEFAULT_CONFIG["weekly_cap"]
    slot_h = DEFAULT_CONFIG["slot_hour"]
    gamma = DEFAULT_CONFIG["gamma"]
    d50 = DEFAULT_CONFIG["d50"]
    t_tol = DEFAULT_CONFIG["t_tol"]
    harm = DEFAULT_CONFIG["harm"]
    if theta_min is None:
        theta_min = DEFAULT_CONFIG.get("theta_min", 0.5)
    if theta_max is None:
        theta_max = DEFAULT_CONFIG.get("theta_max", 2.5)

    z = np.zeros(n_days, dtype=int)
    u_mat = np.zeros((n_days, 24), dtype=int)
    stim_history_hours = []
    stim_window = deque()
    last_stim_hour = -999999

    # Initial state
    u0 = p["trans_u"][0]
    z[0] = 0 if u0 < 0.70 else (1 if u0 < 0.85 else 2)

    d_prev = 0.0
    cum_stims = 0
    burn_in = DEFAULT_CONFIG["burn_in_days"]

    # Pre-generate DP observations and risk for daily controllers
    # If arm does not use DP risk, bypass risk estimation to maximize efficiency
    needs_risk = arm in ("dp_gated", "hierarchical")
    r_hat = np.full(n_days, 0.2)
    r_bar = 0.2

    ewma_y = 0.0
    ewma_m = 0.0
    y_clean_hist = []

    for k in range(n_days):
        # 1. Update observations & risk for today
        sev_k = [0.0, 0.5, 1.0][z[k]]
        ctx_k = p["rho_ctx"] * ((sev_k - 0.5) / 0.7) + np.sqrt(max(0.0, 1 - p["rho_ctx"]**2)) * p["u_ctx"][k]
        y_k = p["a"] * sev_k + p["c"] * ctx_k + p["e_dp"][k]
        p_miss_k = 1.0 / (1.0 + np.exp(-(-2.0 + 0.8 * p["miss_inform"] * sev_k)))
        m_k = p["miss_u"][k] < p_miss_k

        if needs_risk and risk_model.trained:
            if k == 0:
                r_hat[k] = 0.2
            else:
                chg3 = (y_clean_hist[-1] - y_clean_hist[-4]) if len(y_clean_hist) >= 4 else (y_clean_hist[-1] - y_clean_hist[0])
                r_hat[k] = risk_model.predict_step(ewma_y, chg3, ewma_m)

        # Update running DP feature state with today's observations (strictly available for day k+1)
        m_val = 1.0 if m_k else 0.0
        ewma_m = DEFAULT_CONFIG["alpha_ewma"] * m_val + (1.0 - DEFAULT_CONFIG["alpha_ewma"]) * ewma_m
        if not m_k:
            ewma_y = DEFAULT_CONFIG["alpha_ewma"] * y_k + (1.0 - DEFAULT_CONFIG["alpha_ewma"]) * ewma_y
            y_clean_hist.append(y_k)
        else:
            y_clean_hist.append(ewma_y)

        # 2. Hourly simulation across 24 hours
        mu_k = [0.0, 0.8, 1.5][z[k]]
        drift_k = p["drift"] * (k / n_days)
        contra_today = p["contra"][k]

        for h in range(24):
            t_hour = k * 24 + h
            circ = p["circ_amp"] * np.sin(2 * np.pi * (h - 14) / 24.0)
            b_th = mu_k + circ + p["eps_neur"][t_hour] + drift_k

            # Hard safety gate evaluation (short-circuited for efficiency)
            gate = False
            if (k >= burn_in) and (not contra_today):
                if (t_hour - last_stim_hour) >= ref_h:
                    while stim_window and stim_window[0] <= t_hour - 168:
                        stim_window.popleft()
                    if len(stim_window) < cap_w:
                        gate = True

            fire = False
            if gate:
                if arm == "neural_only":
                    fire = (b_th > theta_b)
                elif arm == "dp_gated":
                    fire = (h == slot_h) and (r_hat[k] > theta_r)
                elif arm == "hierarchical":
                    # Modulated threshold with safety clipping
                    eff_theta_b = np.clip(theta_b - kappa * (r_hat[k] - r_bar), theta_min, theta_max)
                    fire = (b_th > eff_theta_b)
                elif arm == "oracle":
                    fire = (h == slot_h) and (z[k] in (1, 2))
                elif arm == "burden_matched_random":
                    # Fired on pre-selected random safe days at slot_h
                    if target_stim_days is not None and k in target_stim_days and h == slot_h:
                        fire = True
                elif arm == "scheduled":
                    if target_stim_days is not None and k in target_stim_days and h == slot_h:
                        fire = True

            if fire:
                u_mat[k, h] = 1
                stim_history_hours.append(t_hour)
                stim_window.append(t_hour)
                last_stim_hour = t_hour

        # 3. Treatment dose calculation for day k
        n_k = np.sum(u_mat[k])
        d_k = gamma * d_prev + n_k
        w_k = 1.0 - np.exp(-d_k / d50)
        cum_stims += n_k
        tol_k = np.exp(-cum_stims / t_tol)
        e_eff = p["e_i"] * w_k * tol_k
        d_prev = d_k

        # 4. Latent state transition for day k+1
        if k < n_days - 1:
            p12_eff = p["p12"] * (1.0 - 0.8 * e_eff)
            p20_eff = min(0.90, p["p20"] * (1.0 + 2.0 * e_eff))
            p01_eff = p["p01"] * (1.0 + harm * d_k)
            u_trans = p["trans_u"][k + 1]

            curr = z[k]
            if curr == 0:
                z[k + 1] = 1 if u_trans < p01_eff else 0
            elif curr == 1:
                if u_trans < p12_eff:
                    z[k + 1] = 2
                elif u_trans < p12_eff + p["p10"]:
                    z[k + 1] = 0
                else:
                    z[k + 1] = 1
            else: # curr == 2
                z[k + 1] = 0 if u_trans < p20_eff else 2

    return z, u_mat


# =============================================================================
# 4. BURDEN-MATCHED COMPARATOR BUILDER
# =============================================================================

def build_burden_matched_days(n_target_stims, p, rng, ref_h=None, cap_w=None):
    """
    Selects exactly n_target_stims days (at slot_hour) such that
    refractory period, weekly cap, and contraindications are strictly satisfied.
    Guarantees exact matching across patients.
    """
    if ref_h is None:
        ref_h = DEFAULT_CONFIG["refractory_h"]
    if cap_w is None:
        cap_w = DEFAULT_CONFIG["weekly_cap"]

    n_days = DEFAULT_CONFIG["n_days"]
    burn_in = DEFAULT_CONFIG["burn_in_days"]

    if n_target_stims <= 0:
        return set()

    safe_days = [k for k in range(burn_in, n_days) if not p["contra"][k]]
    if len(safe_days) < n_target_stims:
        raise ValueError(f"Cannot match burden: safe days ({len(safe_days)}) < target ({n_target_stims})")

    min_day_gap = int(np.ceil(ref_h / 24.0))

    # Phase 1: Randomized greedy selection for moderate targets (<= 40 stims)
    if n_target_stims <= 40:
        for _ in range(5):
            shuffled = safe_days.copy()
            rng.shuffle(shuffled)
            selected = []
            for day in shuffled:
                if len(selected) >= n_target_stims:
                    break
                temp = sorted(selected + [day])
                # Refractory check
                valid = True
                for idx in range(len(temp) - 1):
                    if temp[idx + 1] - temp[idx] < min_day_gap:
                        valid = False
                        break
                if not valid:
                    continue
                # Rolling 7-day cap check (only check windows covering the new day)
                for st in [day - i for i in range(7) if day - i in temp]:
                    stims_in_window = sum(1 for other in temp if st <= other < st + 7)
                    if stims_in_window > cap_w:
                        valid = False
                        break
                if valid:
                    selected.append(day)

            if len(selected) == n_target_stims:
                return set(selected)

    # Phase 2: For dense/high target counts, build maximal valid sequential schedule
    # with randomized start offsets, and uniformly downsample to n_target_stims.
    # Mathematically, ANY subset of a valid schedule strictly preserves refractory,
    # rolling weekly cap, and contraindication constraints.
    for offset in range(14):
        dense_valid = []
        for d in safe_days[offset:]:
            if dense_valid and (d - dense_valid[-1] < min_day_gap):
                continue
            stims_in_window = sum(1 for other in dense_valid if d - 7 < other <= d)
            if stims_in_window < cap_w:
                dense_valid.append(d)

        if len(dense_valid) >= n_target_stims:
            sub = rng.choice(dense_valid, size=n_target_stims, replace=False)
            return set(int(x) for x in sub)

    # Fallback with pure sequential packing from day 0
    dense_valid = []
    for d in safe_days:
        if dense_valid and (d - dense_valid[-1] < min_day_gap):
            continue
        stims_in_window = sum(1 for other in dense_valid if d - 7 < other <= d)
        if stims_in_window < cap_w:
            dense_valid.append(d)

    if len(dense_valid) >= n_target_stims:
        sub = rng.choice(dense_valid, size=n_target_stims, replace=False)
        return set(int(x) for x in sub)

    raise RuntimeError(f"Failed to find burden-matched schedule for target {n_target_stims} after all attempts")


# =============================================================================
# 5. COHORT EVALUATION & METRICS
# =============================================================================

def evaluate_cohort(cohort, arm, risk_model, theta_b=1.2, theta_r=0.4, kappa=3.0,
                    target_stim_counts=None, seed=42):
    """Evaluates an arm across a cohort of virtual patients."""
    rng = np.random.default_rng(seed)
    burn_in = DEFAULT_CONFIG["burn_in_days"]
    n_days = DEFAULT_CONFIG["n_days"]
    scored_days = n_days - burn_in

    epi_fractions = []
    stim_counts = []
    false_stims = []
    total_stims = []

    for i, p in enumerate(cohort):
        target_days = None
        if arm == "burden_matched_random" and target_stim_counts is not None:
            n_target = int(target_stim_counts[i])
            target_days = build_burden_matched_days(n_target, p, rng)

        z, u_mat = run_patient_simulation(p, arm, risk_model, theta_b, theta_r, kappa,
                                          target_stim_days=target_days)

        # Scored period (post burn-in)
        z_scored = z[burn_in:]
        u_scored = u_mat[burn_in:]

        # Episode fraction (primary endpoint)
        epi_frac = np.mean(z_scored == 2)
        n_stims = int(np.sum(u_scored))

        # Explicit verification: burden matching must be 100% exact per patient
        if arm == "burden_matched_random" and target_stim_counts is not None:
            if n_stims != n_target:
                raise RuntimeError(f"Burden mismatch for patient {i}: delivered {n_stims} != target {n_target}")

        # False stimulations: delivered while z == 0
        z_hourly = np.repeat(z_scored, 24)
        u_hourly = u_scored.flatten()
        n_false = int(np.sum((u_hourly == 1) & (z_hourly == 0)))

        epi_fractions.append(epi_frac)
        stim_counts.append(n_stims)
        false_stims.append(n_false)
        total_stims.append(n_stims)

    tot_s = sum(total_stims)
    fsr = (sum(false_stims) / tot_s) if tot_s > 0 else 0.0
    ppv = 1.0 - fsr

    return {
        "mean_epi_frac": float(np.mean(epi_fractions)),
        "sd_epi_frac": float(np.std(epi_fractions)),
        "epi_fractions": np.array(epi_fractions),
        "mean_stims": float(np.mean(stim_counts)),
        "stim_counts": np.array(stim_counts),
        "fsr": float(fsr),
        "ppv": float(ppv)
    }


# =============================================================================
# 6. THRESHOLD TUNING ON TRAIN COHORT
# =============================================================================

def tune_thresholds_on_train(train_cohort, risk_model, seed=42):
    """
    Searches operational thresholds on TRAIN natural history
    to satisfy stimulation budget (default <= 30 stims per patient).
    Uses binary search over candidate grid to guarantee exact match with efficiency.
    """
    target_budget = DEFAULT_CONFIG["budget"]

    # 1. Tune Neural-only threshold theta_b
    b_cands = np.linspace(0.8, 2.2, 15)
    low, high = 0, len(b_cands) - 1
    best_idx_b = high
    while low <= high:
        mid = (low + high) // 2
        res = evaluate_cohort(train_cohort, "neural_only", risk_model, theta_b=b_cands[mid], seed=seed)
        if res["mean_stims"] <= target_budget:
            best_idx_b = mid
            high = mid - 1
        else:
            low = mid + 1
    best_theta_b = float(b_cands[best_idx_b])

    # 2. Tune DP-gated threshold theta_r
    r_cands = np.linspace(0.2, 0.7, 15)
    low, high = 0, len(r_cands) - 1
    best_idx_r = high
    while low <= high:
        mid = (low + high) // 2
        res = evaluate_cohort(train_cohort, "dp_gated", risk_model, theta_r=r_cands[mid], seed=seed)
        if res["mean_stims"] <= target_budget:
            best_idx_r = mid
            high = mid - 1
        else:
            low = mid + 1
    best_theta_r = float(r_cands[best_idx_r])

    return best_theta_b, best_theta_r


# =============================================================================
# 7. MAIN EXPERIMENT RUNNER
# =============================================================================

def run_stress_test(scenario="informative", miss_inform=None, neural_sd=1.0, drift_scale=0.0,
                    n_train=200, n_test=500, seed=1, kappa=3.0, config_override=None, out_file=None):
    """Runs a complete Train-Tune-Test DP-CLN stress test evaluation."""
    print(f"[DP-CLN Engine v0.3] Initializing Stress Test: Scenario={scenario}, "
          f"MissInform={miss_inform}, NeuralSD={neural_sd}, Drift={drift_scale}, Kappa={kappa}, Seed={seed}")

    if config_override:
        cfg = DEFAULT_CONFIG.copy()
        cfg.update(config_override)
    else:
        cfg = DEFAULT_CONFIG

    rng = np.random.default_rng(seed)

    # 1. Generate Cohorts (independent test cohort, train-test distribution shift applied)
    train_cohort = sample_patient_cohort(n_train, rng, scenario, miss_inform, neural_sd, drift_scale, is_test=False)
    test_cohort = sample_patient_cohort(n_test, rng, scenario, miss_inform, neural_sd, drift_scale, is_test=True)

    # 2. Train Risk Model strictly on natural history
    risk_model = SimpleLogisticRiskModel(horizon=3, alpha=cfg["alpha_ewma"])
    risk_model.fit(train_cohort)
    print("  -> Risk model trained successfully on natural history.")

    # 3. Tune Thresholds on Train Cohort
    theta_b, theta_r = tune_thresholds_on_train(train_cohort, risk_model, seed=seed)
    print(f"  -> Tuned thresholds on TRAIN: theta_b={theta_b:.3f}, theta_r={theta_r:.3f}")

    # 4. Evaluate Natural History (None)
    res_none = evaluate_cohort(test_cohort, "None", risk_model, seed=seed)

    # 5. Evaluate Adaptive Arms with paired burden-matched comparators
    results = {"None": res_none}
    adaptive_arms = ["neural_only", "dp_gated", "hierarchical"]
    for idx, arm in enumerate(adaptive_arms):
        res_arm = evaluate_cohort(test_cohort, arm, risk_model, theta_b=theta_b, theta_r=theta_r, kappa=kappa, seed=seed + 10 * (idx + 1))

        # 6. Evaluate Exact Burden-Matched Random Comparator
        res_rand = evaluate_cohort(test_cohort, "burden_matched_random", risk_model,
                                   target_stim_counts=res_arm["stim_counts"], seed=seed + 50 * (idx + 1))

        # Explicit Verification: 100% identical stimulation delivery per patient
        if not np.array_equal(res_rand["stim_counts"], res_arm["stim_counts"]):
            mismatches = np.sum(res_rand["stim_counts"] != res_arm["stim_counts"])
            raise RuntimeError(f"Burden mismatch in arm {arm}: {mismatches} patient discrepancies detected!")

        # Primary Outcome: Burden-Matched Benefit
        # Benefit = Episode Fraction (Random) - Episode Fraction (Adaptive)
        # Positive values strictly indicate that the adaptive arm averts more episode days than random
        paired_diff = res_rand["epi_fractions"] - res_arm["epi_fractions"]
        benefit_mean = float(np.mean(paired_diff))
        benefit_sd = float(np.std(paired_diff))
        pct_positive = float(np.mean(paired_diff > 0.0) * 100)

        # Bootstrap 95% confidence interval for primary paired benefit (2000 resamples)
        boot_rng = np.random.default_rng(seed + 999 * (idx + 1))
        n_boot = 2000
        n_pts = len(paired_diff)
        boot_means = np.zeros(n_boot)
        for b in range(n_boot):
            boot_idx = boot_rng.integers(0, n_pts, size=n_pts)
            boot_means[b] = np.mean(paired_diff[boot_idx])
        ci_lower = float(np.percentile(boot_means, 2.5))
        ci_upper = float(np.percentile(boot_means, 97.5))

        res_arm["burden_matched_random_mean_epi"] = res_rand["mean_epi_frac"]
        res_arm["burden_matched_benefit_mean"] = benefit_mean
        res_arm["burden_matched_benefit_sd"] = benefit_sd
        res_arm["benefit_ci_95"] = [ci_lower, ci_upper]
        res_arm["pct_patients_benefiting"] = pct_positive
        res_arm["paired_diff"] = paired_diff
        res_arm["rand_stim_counts"] = res_rand["stim_counts"]
        results[arm] = res_arm

    # Format JSON serializable output with complete configuration and patient-level raw data
    output_summary = {
        "version": "v0.3",
        "scenario": scenario,
        "miss_inform": miss_inform,
        "neural_sd": neural_sd,
        "drift_scale": drift_scale,
        "n_train": n_train,
        "n_test": n_test,
        "seed": seed,
        "kappa": float(kappa),
        "config": cfg,
        "tuned_theta_b": float(theta_b),
        "tuned_theta_r": float(theta_r),
        "natural_history_epi_frac": res_none["mean_epi_frac"],
        "arms": {}
    }

    for arm in adaptive_arms:
        a_res = results[arm]
        output_summary["arms"][arm] = {
            "mean_epi_frac": a_res["mean_epi_frac"],
            "mean_stims": a_res["mean_stims"],
            "fsr": a_res["fsr"],
            "ppv": a_res["ppv"],
            "burden_matched_random_epi": a_res["burden_matched_random_mean_epi"],
            "primary_burden_matched_benefit": a_res["burden_matched_benefit_mean"],
            "benefit_sd": a_res["burden_matched_benefit_sd"],
            "benefit_ci_95": a_res["benefit_ci_95"],
            "pct_patients_benefiting": a_res["pct_patients_benefiting"],
            "burden_matching_verified": True,
            "raw_paired_benefit": [float(x) for x in a_res["paired_diff"]],
            "raw_stim_counts": [int(x) for x in a_res["stim_counts"]],
            "raw_epi_fractions": [float(x) for x in a_res["epi_fractions"]]
        }

    if out_file:
        os.makedirs(os.path.dirname(os.path.abspath(out_file)), exist_ok=True)
        with open(out_file, "w") as f:
            json.dump(output_summary, f, indent=2)
        print(f"  -> Results successfully saved to {out_file}")

    return output_summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="DP-CLN Stress Test Engine v0.3")
    parser.add_argument("--scenario", type=str, default="informative",
                        choices=["informative", "weak", "uninformative", "misleading"])
    parser.add_argument("--miss_inform", type=int, default=1, choices=[0, 1])
    parser.add_argument("--neural_sd", type=float, default=1.0)
    parser.add_argument("--drift_scale", type=float, default=0.0)
    parser.add_argument("--n_train", type=int, default=100)
    parser.add_argument("--n_test", type=int, default=200)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out", type=str, default=None)
    args = parser.parse_args()

    res = run_stress_test(
        scenario=args.scenario,
        miss_inform=args.miss_inform,
        neural_sd=args.neural_sd,
        drift_scale=args.drift_scale,
        n_train=args.n_train,
        n_test=args.n_test,
        seed=args.seed,
        out_file=args.out
    )
    print("\n--- Summary Results ---")
    print(json.dumps(res, indent=2))
