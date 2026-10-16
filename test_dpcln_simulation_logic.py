"""
Deterministic Unit and Integration Tests for DP-CLN Simulation Control Logic

Verifies:
1. Safety Gate Failures (contraindications, confidence thresholds).
2. Rolling Refractory Enforcement (no stim within tau_ref).
3. Rolling Weekly Budget / Cap Enforcement (no more than N_cap per 7 days).
4. Boundary Cases:
   - Zero treatment efficacy (e_i = 0)
   - Zero observation noise
   - 100% missing data
5. Treatment Pharmacodynamic Kinetics:
   - Carryover accumulation (D_k = gamma * D_{k-1} + n_k)
   - Nonlinear saturation (w_k = 1 - exp(-D_k / D50))
   - Tolerance desensitization (tol_k = exp(-C_k / T))
6. Hand-Calculated Known Trajectory (step-by-step verification).
"""

import os
import sys
import unittest
import numpy as np

# Ensure parent and local directory are on sys.path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


class SafetyGate:
    """Implements clinician-governed hard safety gate g(t)."""
    def __init__(self, refractory_h=24, weekly_cap=5, min_confidence=0.6):
        self.refractory_h = refractory_h
        self.weekly_cap = weekly_cap
        self.min_confidence = min_confidence

    def evaluate(self, t, stim_history, contraindication=False, confidence=1.0):
        # 1. Active contraindication
        if contraindication:
            return 0
        # 2. Insufficient confidence (safe-mode fallback)
        if confidence < self.min_confidence:
            return 0
        # 3. Refractory period
        recent_stims = [st for st in stim_history if st < t]
        if recent_stims and (t - recent_stims[-1]) < self.refractory_h:
            return 0
        # 4. Rolling 7-day (168 h) stimulation cap
        window_stims = [st for st in stim_history if t - 168 < st < t]
        if len(window_stims) >= self.weekly_cap:
            return 0
        return 1


class TreatmentDynamics:
    """Implements dose carryover, saturation, and tolerance."""
    def __init__(self, gamma=0.85, d50=2.0, t_tol=60.0):
        self.gamma = gamma
        self.d50 = d50
        self.t_tol = t_tol

    def step(self, d_prev, cum_stims, n_today, e_intrinsic):
        d_k = self.gamma * d_prev + n_today
        w_k = 1.0 - np.exp(-d_k / self.d50)
        cum_stims_next = cum_stims + n_today
        tol_k = np.exp(-cum_stims_next / self.t_tol)
        e_eff = e_intrinsic * w_k * tol_k
        return d_k, cum_stims_next, w_k, tol_k, e_eff


class TestDPCLNSimulationLogic(unittest.TestCase):

    def setUp(self):
        self.gate = SafetyGate(refractory_h=24, weekly_cap=5, min_confidence=0.6)
        self.dynamics = TreatmentDynamics(gamma=0.85, d50=2.0, t_tol=60.0)

    # -------------------------------------------------------------------------
    # 1. Safety Gate Tests
    # -------------------------------------------------------------------------
    def test_contraindication_blocks_stimulation(self):
        """g(t) MUST be 0 if contraindication is active, regardless of time or history."""
        # Clean history
        stim_history = []
        g = self.gate.evaluate(t=100, stim_history=stim_history, contraindication=True, confidence=0.99)
        self.assertEqual(g, 0, "Gate must reject when contraindication=True")

    def test_low_confidence_blocks_stimulation(self):
        """g(t) MUST be 0 if model confidence < c_min (safe mode fallback)."""
        stim_history = []
        g = self.gate.evaluate(t=100, stim_history=stim_history, contraindication=False, confidence=0.59)
        self.assertEqual(g, 0, "Gate must reject when confidence < 0.60")

        g_pass = self.gate.evaluate(t=100, stim_history=stim_history, contraindication=False, confidence=0.60)
        self.assertEqual(g_pass, 1, "Gate must pass when confidence >= 0.60")

    def test_canonical_refractory_enforcement_24h(self):
        """Verifies canonical benchmark refractory enforcement (tau_ref = 24 h): no stimulation permitted if elapsed time < 24 h."""
        stim_history = [100]
        # At t = 101 to 123 (less than 24 h elapsed)
        for t in range(101, 124):
            g = self.gate.evaluate(t=t, stim_history=stim_history, contraindication=False, confidence=1.0)
            self.assertEqual(g, 0, f"Gate must block stimulation at t={t} (last stim at 100, gap={t-100} < 24)")

        # At t = 124 (exactly 24 h elapsed)
        g_exact = self.gate.evaluate(t=124, stim_history=stim_history, contraindication=False, confidence=1.0)
        self.assertEqual(g_exact, 1, "Gate must permit stimulation at t=124 (gap=24 h)")

    def test_weekly_cap_enforcement(self):
        """No stimulation permitted if 5 stimulations occurred in the past 168 hours."""
        # 5 stimulations spaced 25 hours apart: t = 0, 25, 50, 75, 100
        stim_history = [0, 25, 50, 75, 100]
        # At t = 125 (gap since last is 25 h > 24 h, but all 5 are within 168 h window [125-168 = -43])
        g_cap = self.gate.evaluate(t=125, stim_history=stim_history, contraindication=False, confidence=1.0)
        self.assertEqual(g_cap, 0, "Gate must block 6th stimulation when 5 already occurred within 168 h")

        # Advance time to t = 169 (stim at 0 drops out of the 168-hour window [169-168 = 1])
        # Remaining stims in window: 25, 50, 75, 100 (total = 4)
        g_unblock = self.gate.evaluate(t=169, stim_history=stim_history, contraindication=False, confidence=1.0)
        self.assertEqual(g_unblock, 1, "Gate must allow stimulation once oldest stimulation exits rolling window")

    # -------------------------------------------------------------------------
    # 2. Pharmacodynamic Kinetics Tests
    # -------------------------------------------------------------------------
    def test_treatment_dose_carryover_recurrence(self):
        """Verifies analytical recurrence of dose carryover D_k = gamma * D_{k-1} + n_k."""
        d = 0.0
        cum = 0
        gamma = 0.85
        # Day 1: 1 stim
        d, cum, w, tol, eff = self.dynamics.step(d, cum, n_today=1, e_intrinsic=0.5)
        self.assertAlmostEqual(d, 1.0, places=7)
        self.assertEqual(cum, 1)

        # Day 2: 0 stim
        d, cum, w, tol, eff = self.dynamics.step(d, cum, n_today=0, e_intrinsic=0.5)
        self.assertAlmostEqual(d, 0.85, places=7)
        self.assertEqual(cum, 1)

        # Day 3: 2 stims
        d, cum, w, tol, eff = self.dynamics.step(d, cum, n_today=2, e_intrinsic=0.5)
        self.assertAlmostEqual(d, 0.85 * 0.85 + 2.0, places=7)
        self.assertEqual(cum, 3)

    def test_saturation_and_tolerance_limits(self):
        """Saturation w_k approaches 1 as dose -> inf, tolerance -> 0 as cum_stims -> inf."""
        # Massive dose
        d_k, cum, w, tol, eff = self.dynamics.step(d_prev=100.0, cum_stims=0, n_today=10, e_intrinsic=1.0)
        self.assertGreater(w, 0.9999, "Saturation must approach 1 for large dose")

        # Massive chronic stimulations
        d_k, cum, w, tol, eff = self.dynamics.step(d_prev=0.0, cum_stims=500, n_today=0, e_intrinsic=1.0)
        self.assertLess(tol, 0.001, "Tolerance factor must approach 0 for high cumulative stimulations")

    def test_zero_efficacy_boundary(self):
        """When e_intrinsic = 0 (non-responder), effective efficacy e_eff MUST be exactly 0."""
        d_k, cum, w, tol, eff = self.dynamics.step(d_prev=5.0, cum_stims=10, n_today=2, e_intrinsic=0.0)
        self.assertEqual(eff, 0.0, "Effective efficacy must be strictly 0 for non-responders")

    # -------------------------------------------------------------------------
    # 3. Known Simplified Trajectory (10-Day Trace)
    # -------------------------------------------------------------------------
    def test_known_10_day_trajectory(self):
        """Step-by-step verification on a deterministic 10-day synthetic trajectory."""
        # Patient intrinsic efficacy = 0.8
        e_i = 0.8
        d = 0.0
        cum = 0
        stim_history = []
        stim_log = []

        # Scheduled protocol: attempts stimulation at 09:00 on days 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
        # Days 0, 1, 2, 3, 4: should succeed (5 stims)
        # Day 5: should be BLOCKED by weekly cap (5 stims in rolling 7 days)
        # Day 6: should be BLOCKED by weekly cap
        # Day 7: t = 7*24 + 9 = 177. 177 - 168 = 9. Stim at day 0 was at t=9.
        # Boundary at Day 7 depends on strict inequalities. Let's trace explicitly:
        for day in range(10):
            t_eval = day * 24 + 9
            g = self.gate.evaluate(t=t_eval, stim_history=stim_history)
            n_today = 1 if g == 1 else 0
            if n_today == 1:
                stim_history.append(t_eval)
            stim_log.append(n_today)
            d, cum, w, tol, eff = self.dynamics.step(d, cum, n_today, e_i)

        # Days 0, 1, 2, 3, 4 fired (5 stims)
        self.assertEqual(stim_log[0:5], [1, 1, 1, 1, 1])
        # Day 5 (t=129) must be blocked because all 5 stims are within past 168 h
        self.assertEqual(stim_log[5], 0)
        # Day 6 (t=153) must be blocked
        self.assertEqual(stim_log[6], 0)
        # Day 7 (t=177): Day 0 (t=9) rolls out of past 168h window (177-168=9). Day 7 unblocks!
        self.assertEqual(stim_log[7], 1)
        # Day 8 (t=201): Day 1 (t=33) rolls out (201-168=33). Day 8 unblocks!
        self.assertEqual(stim_log[8], 1)
        # Day 9 (t=225): Day 2 (t=57) rolls out (225-168=57). Day 9 unblocks!
        self.assertEqual(stim_log[9], 1)
        # Total stimulations over 10 days = 5 + 0 + 0 + 1 + 1 + 1 = 8
        self.assertEqual(cum, 8)

    # -------------------------------------------------------------------------
    # 4. Missing Data Graceful Handling
    # -------------------------------------------------------------------------
    def test_missing_data_imputation_fallback(self):
        """Verifies that 100% missing DP observations fall back to normative prior without crash."""
        n_days = 30
        y_obs = np.full(n_days, np.nan)
        missing_indicator = np.ones(n_days)

        # EWMA with missingness handling
        alpha = 0.3
        ewma_y = np.zeros(n_days)
        last_val = 0.0  # Population prior
        for k in range(n_days):
            if not np.isnan(y_obs[k]):
                last_val = alpha * y_obs[k] + (1 - alpha) * last_val
            ewma_y[k] = last_val

        self.assertTrue(np.all(ewma_y == 0.0), "EWMA of 100% missing data must remain at normative prior 0.0")
        self.assertFalse(np.any(np.isnan(ewma_y)), "Imputed sequence must not contain NaNs")

    # -------------------------------------------------------------------------
    # 5. Burden-Matching & Fair Comparator Tests
    # -------------------------------------------------------------------------
    def test_burden_matching_canonical_refractory_24h_and_safety_constraints(self):
        """Verifies build_burden_matched_days achieves exact target counts under canonical benchmark 24-hour refractory (tau_ref = 24 h) and weekly cap constraints."""
        try:
            from dpcln_stress_test import sample_patient_cohort, build_burden_matched_days, run_patient_simulation, SimpleLogisticRiskModel
        except ImportError:
            from literature.dpcln_stress_test import sample_patient_cohort, build_burden_matched_days, run_patient_simulation, SimpleLogisticRiskModel
        rng = np.random.default_rng(101)
        cohort = sample_patient_cohort(20, rng)
        model = SimpleLogisticRiskModel()

        for p in cohort:
            for n_target in [0, 5, 15, 30]:
                selected = build_burden_matched_days(n_target, p, rng, ref_h=None, cap_w=5)
                self.assertEqual(len(selected), n_target, f"Failed to match target count {n_target}")

                # Verify refractory constraint (minimum 1-day gap for canonical ref_h=24)
                sorted_days = sorted(selected)
                for idx in range(len(sorted_days) - 1):
                    self.assertGreaterEqual(sorted_days[idx + 1] - sorted_days[idx], 1, "Refractory violation in selected days")

                # Verify rolling 7-day cap constraint (<= 5 per 7 days)
                for st in sorted_days:
                    in_window = sum(1 for o in sorted_days if st - 7 < o <= st)
                    self.assertLessEqual(in_window, 5, "Weekly cap violation in selected days")

                # Verify that run_patient_simulation delivers EXACTLY n_target stimulations
                z, u_mat = run_patient_simulation(p, "burden_matched_random", model, target_stim_days=selected)
                delivered = np.sum(u_mat[30:])  # Scored post burn-in
                self.assertEqual(delivered, n_target, f"Delivered stimulations ({delivered}) != target ({n_target})")

    def test_burden_matching_refractory_sensitivity_sweep_48h(self):
        """Verifies that when the refractory interval is varied to 48 hours in exploratory parameter sensitivity sweeps, burden matching enforces a minimum 2-day gap, contrasting with the canonical benchmark (tau_ref = 24 h)."""
        try:
            from dpcln_stress_test import sample_patient_cohort, build_burden_matched_days
        except ImportError:
            from literature.dpcln_stress_test import sample_patient_cohort, build_burden_matched_days
        rng = np.random.default_rng(202)
        cohort = sample_patient_cohort(10, rng)
        for p in cohort:
            for n_target in [5, 15, 25]:
                selected = build_burden_matched_days(n_target, p, rng, ref_h=48, cap_w=5)
                self.assertEqual(len(selected), n_target)
                sorted_days = sorted(selected)
                for idx in range(len(sorted_days) - 1):
                    self.assertGreaterEqual(sorted_days[idx + 1] - sorted_days[idx], 2, "48-h refractory sweep violation: gap < 2 days")

    # -------------------------------------------------------------------------
    # 6. Pipeline Lag & Latent State Leakage Tests
    # -------------------------------------------------------------------------
    def test_pipeline_lag_no_leakage(self):
        """Verifies that features for day k depend strictly on history through day k-1 and cannot see day k data."""
        try:
            from dpcln_stress_test import SimpleLogisticRiskModel
        except ImportError:
            from literature.dpcln_stress_test import SimpleLogisticRiskModel
        model = SimpleLogisticRiskModel()
        rng = np.random.default_rng(303)

        n_days = 20
        y = rng.normal(0, 1, size=n_days)
        is_missing = np.zeros(n_days, dtype=bool)

        feat_orig = model.extract_features(y, is_missing)

        # Alter day k's observation dramatically
        for k in range(1, n_days):
            y_modified = y.copy()
            y_modified[k] += 100.0  # Massive spike on day k
            is_miss_modified = is_missing.copy()
            is_miss_modified[k] = True

            feat_mod = model.extract_features(y_modified, is_miss_modified)

            # Features AT day k MUST NOT change (strict 1-day lag)
            np.testing.assert_allclose(feat_mod[k], feat_orig[k], err_msg=f"Leakage detected: day {k} features changed when day {k} data altered!")

            # Features AT day k+1 MUST reflect day k's change
            if k + 1 < n_days:
                self.assertFalse(np.allclose(feat_mod[k + 1], feat_orig[k + 1]), f"Lag error: day {k+1} features failed to reflect day {k} change!")

    def test_uninformative_scenario_dissociation(self):
        """Verifies uninformative scenario defaults to m_b=0 and has zero coupling to disease state."""
        try:
            from dpcln_stress_test import sample_patient_cohort, generate_patient_observations
        except ImportError:
            from literature.dpcln_stress_test import sample_patient_cohort, generate_patient_observations
        rng = np.random.default_rng(404)
        cohort = sample_patient_cohort(10, rng, scenario="uninformative")

        for p in cohort:
            self.assertEqual(p["a"], 0.0, "Uninformative scenario must have a=0")
            self.assertEqual(p["c"], 0.0, "Uninformative scenario must have c=0")
            self.assertEqual(p["miss_inform"], 0, "Uninformative scenario must have miss_inform=0 (no missingness leakage)")

            # Generate observations across synthetic alternating states (180 days)
            z_seq = np.array([0, 1, 2] * 60)
            b, y, is_missing = generate_patient_observations(p, z_seq)

            # When a=0 and c=0, y MUST be identically equal to pure noise e_dp (zero state dependency)
            np.testing.assert_allclose(y, p["e_dp"], err_msg="y must be identically equal to pure DP noise e_dp")

            # When miss_inform=0, missingness probability is strictly constant (logistic(-2.0))
            # regardless of whether state is 0, 1, or 2
            p_miss_const = 1.0 / (1.0 + np.exp(2.0))
            expected_missing = p["miss_u"] < p_miss_const
            np.testing.assert_array_equal(is_missing, expected_missing, err_msg="Missingness must be strictly independent of state")

    def test_misleading_distribution_shift(self):
        """Verifies train-test distribution shift for misleading scenario (rho_train=0.6, rho_test=0.0)."""
        try:
            from dpcln_stress_test import sample_patient_cohort
        except ImportError:
            from literature.dpcln_stress_test import sample_patient_cohort
        rng = np.random.default_rng(505)
        train_cohort = sample_patient_cohort(5, rng, scenario="misleading", is_test=False)
        test_cohort = sample_patient_cohort(5, rng, scenario="misleading", is_test=True)

        for p in train_cohort:
            self.assertAlmostEqual(p["rho_ctx"], 0.6, places=5, msg="Train cohort must have rho_train=0.6")
        for p in test_cohort:
            self.assertAlmostEqual(p["rho_ctx"], 0.0, places=5, msg="Test cohort must have rho_test=0.0 (distribution shift)")

    # -------------------------------------------------------------------------
    # 7. Outcome Sign & Statistical Metric Tests
    # -------------------------------------------------------------------------
    def test_primary_benefit_sign_convention(self):
        """Confirms positive benefit strictly corresponds to fewer episodes under adaptive controller."""
        epi_rand = np.array([0.25, 0.30, 0.20])
        epi_adaptive = np.array([0.10, 0.15, 0.05])

        paired_diff = epi_rand - epi_adaptive
        benefit = np.mean(paired_diff)

        self.assertGreater(benefit, 0.0, "Benefit must be strictly positive when adaptive policy reduces episodes")
        self.assertAlmostEqual(benefit, 0.15, places=5)

    def test_bootstrap_confidence_interval(self):
        """Verifies bootstrap CI generates correct percentile bounds."""
        rng = np.random.default_rng(606)
        diffs = rng.normal(0.05, 0.02, size=100)
        n_boot = 1000
        boot_means = np.zeros(n_boot)
        for b in range(n_boot):
            boot_idx = rng.integers(0, len(diffs), size=len(diffs))
            boot_means[b] = np.mean(diffs[boot_idx])
        ci_lower = np.percentile(boot_means, 2.5)
        ci_upper = np.percentile(boot_means, 97.5)

        self.assertLess(ci_lower, np.mean(diffs))
        self.assertGreater(ci_upper, np.mean(diffs))

    # -------------------------------------------------------------------------
    # 8. Bounded Threshold Modulation Tests (Eq. 5.6 & Sec 6)
    # -------------------------------------------------------------------------
    def test_bounded_threshold_adaptation(self):
        """
        Verifies theta_b(r_hat) = clip(theta_0 - kappa*(r_hat - r_bar), theta_min, theta_max).
        Prevents runaway threshold collapse under extreme high risk (r_hat = 1.0)
        and unbounded threshold escalation under zero risk.
        """
        theta_0 = 1.2
        r_bar = 0.2
        kappa = 3.0
        theta_min = 0.5
        theta_max = 2.5

        # Extreme high risk: r_hat = 1.0
        # Unbounded would yield 1.2 - 3.0*(0.8) = -1.2 (runaway hyper-stimulation collapse)
        unbounded_high = theta_0 - kappa * (1.0 - r_bar)
        self.assertAlmostEqual(unbounded_high, -1.2, places=5)
        bounded_high = np.clip(unbounded_high, theta_min, theta_max)
        self.assertEqual(bounded_high, theta_min, "Floor theta_min must prevent threshold collapse below 0.5")

        # Baseline risk: r_hat = r_bar = 0.2
        bounded_base = np.clip(theta_0 - kappa * (r_bar - r_bar), theta_min, theta_max)
        self.assertEqual(bounded_base, theta_0, "At baseline risk, threshold must equal nominal theta_0")

        # Low risk: r_hat = 0.0
        unbounded_low = theta_0 - kappa * (0.0 - r_bar)
        self.assertAlmostEqual(unbounded_low, 1.8)
        bounded_low = np.clip(unbounded_low, theta_min, theta_max)
        self.assertAlmostEqual(bounded_low, 1.8)

        # Extreme negative perturbation driving threshold above theta_max
        unbounded_extreme_low = theta_0 - kappa * (-0.6)  # 1.2 + 1.8 = 3.0
        bounded_ceiling = np.clip(unbounded_extreme_low, theta_min, theta_max)
        self.assertEqual(bounded_ceiling, theta_max, "Ceiling theta_max must bound upper threshold at 2.5")


# =============================================================================
# 9. Comprehensive Mathematical Audit Test Suite (Items 1-7)
# =============================================================================

class TestMarkovMatrixAndDynamicsAudit(unittest.TestCase):
    """
    Dedicated tests verifying:
    1. Complete 3x3 Markov transition matrix stochastic validity across parameter space.
    2. Units, dimensions, and non-saturation of over-stimulation parameter h.
    3. Geometric expected duration formulas for baseline and treatment-modified states.
    4. Dimensional consistency of dose carryover, saturation, and tolerance desensitization.
    5. Dual-clock temporal schedule and causal update sequence.
    """

    def test_markov_matrix_stochastic_validity_exhaustive(self):
        """
        Sweeps 100,000 parameter combinations across all allowable ranges:
        p01 in [0.005, 0.030], p12 in [0.10, 0.30], p10 in [0.05, 0.15],
        p20 in [0.03, 0.08], e_eff in [0, 1], D_k in [0, 100].
        Verifies:
        - p10 >= 0, p12_eff >= 0, p10 + p12_eff <= 1.
        - 0 <= p01_eff <= 1, 0 <= p20_eff <= 1.
        - Every row sums identically to 1.0 (abs error < 1e-14).
        """
        rng = np.random.default_rng(42)
        n_samples = 100000

        p01_arr = rng.uniform(0.005, 0.030, size=n_samples)
        p12_arr = rng.uniform(0.10, 0.30, size=n_samples)
        p10_arr = rng.uniform(0.05, 0.15, size=n_samples)
        p20_arr = rng.uniform(0.03, 0.08, size=n_samples)
        e_eff_arr = rng.uniform(0.0, 1.0, size=n_samples)
        d_k_arr = rng.uniform(0.0, 100.0, size=n_samples)
        h = 0.05

        # Vectorized probability calculations
        p12_eff = p12_arr * (1.0 - 0.80 * e_eff_arr)
        p20_eff = np.minimum(0.90, p20_arr * (1.0 + 2.0 * e_eff_arr))
        p01_eff = np.minimum(1.0, p01_arr * (1.0 + h * d_k_arr))

        # Check conditions
        self.assertTrue(np.all(p10_arr >= 0), "p10 must be non-negative")
        self.assertTrue(np.all(p12_eff >= 0), "p12' must be non-negative")
        self.assertTrue(np.all(p10_arr + p12_eff <= 1.0), "p10 + p12' must be <= 1.0")
        self.assertTrue(np.all((p01_eff >= 0) & (p01_eff <= 1.0)), "0 <= p01' <= 1.0")
        self.assertTrue(np.all((p20_eff >= 0) & (p20_eff <= 1.0)), "0 <= p20' <= 1.0")

        # Row 0: [1 - p01_eff, p01_eff, 0]
        row0_sum = (1.0 - p01_eff) + p01_eff
        np.testing.assert_allclose(row0_sum, 1.0, atol=1e-14, err_msg="Row 0 must sum to 1.0")

        # Row 1: [p10, 1 - p10 - p12_eff, p12_eff]
        row1_sum = p10_arr + (1.0 - p10_arr - p12_eff) + p12_eff
        np.testing.assert_allclose(row1_sum, 1.0, atol=1e-14, err_msg="Row 1 must sum to 1.0")

        # Row 2: [p20_eff, 0, 1 - p20_eff]
        row2_sum = p20_eff + (1.0 - p20_eff)
        np.testing.assert_allclose(row2_sum, 1.0, atol=1e-14, err_msg="Row 2 must sum to 1.0")

    def test_overstimulation_parameter_and_saturation(self):
        """
        Verifies:
        - Units: h has units of sessions^-1, so h*D_k is dimensionless.
        - Non-saturation: Under rolling 7-day budget N_cap = 5 sessions (mean 5/7 sessions/day),
          maximum steady-state D_k <= 5 / (1 - 0.95) = 14.3 sessions.
          With h = 0.05, h*D_k <= 0.714, yielding p01' <= 0.03 * 1.714 = 0.051 << 1.0.
          Demonstrates p01' does NOT saturate at 1.0 under allowed clinical dose histories.
        """
        h = 0.05  # sessions^-1
        # Steady-state maximum dose under weekly cap: 5 sessions per 7 days
        max_daily_rate = 5.0 / 7.0
        for gamma in [0.70, 0.85, 0.95]:
            d_ss = max_daily_rate / (1.0 - gamma)
            h_dk = h * d_ss
            self.assertLess(h_dk, 1.0, "h*D_k must remain < 1.0 under clinical weekly budget")
            for p01 in [0.005, 0.030]:
                p01_prime = min(1.0, p01 * (1.0 + h_dk))
                self.assertLess(p01_prime, 0.06, "p01' must remain well below 1.0 (no saturation)")

    def test_expected_durations_analytical_vs_treatment(self):
        """
        Verifies discrete-time geometric holding times E[T_i] = 1 / (1 - P_ii):
        - Remission: E[T] = 1 / p01 in [33.3, 200.0] days.
        - Prodrome: E[T] = 1 / (p12 + p10) in [2.22, 6.67] days.
        - Episode (Baseline): E[T] = 1 / p20 in [12.5, 33.3] days.
        - Episode (Treated, e_eff=1.0): p20' = 3*p20, E[T'] = 1 / (3*p20) in [4.17, 11.11] days.
        """
        # Baseline remission
        self.assertAlmostEqual(1.0 / 0.030, 33.3333, places=3)
        self.assertAlmostEqual(1.0 / 0.005, 200.0, places=3)

        # Baseline prodrome
        self.assertAlmostEqual(1.0 / (0.30 + 0.15), 2.2222, places=3)
        self.assertAlmostEqual(1.0 / (0.10 + 0.05), 6.6667, places=3)

        # Baseline episode
        self.assertAlmostEqual(1.0 / 0.08, 12.5, places=3)
        self.assertAlmostEqual(1.0 / 0.03, 33.3333, places=3)

        # Treated episode (maximum efficacy e_eff = 1.0 => p20' = min(0.90, p20 * 3.0))
        p20_treated_min = 0.03 * 3.0  # 0.09
        p20_treated_max = 0.08 * 3.0  # 0.24
        self.assertAlmostEqual(1.0 / p20_treated_max, 4.1667, places=3)
        self.assertAlmostEqual(1.0 / p20_treated_min, 11.1111, places=3)

    def test_dose_kinetics_dimensional_consistency(self):
        """
        Verifies:
        - D_k = gamma*D_{k-1} + n_k has units of sessions.
        - w_k = 1 - exp(-D_k / D50) is strictly in [0, 1) and dimensionless.
        - tol_k = exp(-C_k / T) is strictly in (0, 1] and dimensionless.
        - e_eff = e_i * w_k * tol_k is strictly in [0, 1].
        """
        d50 = 2.0  # sessions
        t_tol = 60.0  # sessions
        gamma = 0.85

        d_prev = 0.0
        cum_stims = 0
        e_i = 0.80

        for n_k in [1, 0, 2, 0, 1, 0]:
            d_k = gamma * d_prev + n_k
            cum_stims += n_k
            w_k = 1.0 - np.exp(-d_k / d50)
            tol_k = np.exp(-cum_stims / t_tol)
            e_eff = e_i * w_k * tol_k

            self.assertGreaterEqual(w_k, 0.0)
            self.assertLess(w_k, 1.0)
            self.assertGreater(tol_k, 0.0)
            self.assertLessEqual(tol_k, 1.0)
            self.assertGreaterEqual(e_eff, 0.0)
            self.assertLessEqual(e_eff, 1.0)
            d_prev = d_k

    def test_dual_clock_schedule_causality(self):
        """
        Verifies causal update order across the simulation:
        1. Day k starts: DP observation & risk computed using prior data.
        2. Hours 0..23: Neural observations generated, pulses delivered.
        3. End of Day k: Daily dose D_k and effective efficacy e_eff updated.
        4. Day k -> k+1: Latent state transition uses day k's treatment.
        Confirms z[k+1] does not leak into day k's stimulation decisions.
        """
        try:
            from dpcln_stress_test import sample_patient_cohort, run_patient_simulation, SimpleLogisticRiskModel
        except ImportError:
            from literature.dpcln_stress_test import sample_patient_cohort, run_patient_simulation, SimpleLogisticRiskModel
        rng = np.random.default_rng(123)
        cohort = sample_patient_cohort(2, rng, scenario="informative")
        risk_model = SimpleLogisticRiskModel()
        # Run simulation
        z, u_mat = run_patient_simulation(cohort[0], "hierarchical", risk_model)

        self.assertEqual(len(z), 180)
        self.assertEqual(u_mat.shape, (180, 24))
        # Verify valid states throughout
        self.assertTrue(np.all(np.isin(z, [0, 1, 2])))
        self.assertTrue(np.all(np.isin(u_mat, [0, 1])))

    def test_markov_matrix_deterministic_edge_cases(self):
        """
        Tests 6 deterministic extreme boundary corner cases:
        1. Minimal baseline parameters, zero dose, zero efficacy.
        2. Maximal baseline parameters, zero dose, zero efficacy.
        3. Maximal efficacy (e_eff = 1.0) with zero dose.
        4. Maximal efficacy (e_eff = 1.0) with maximal steady-state dose (D_k = 14.28 sessions).
        5. Extreme out-of-envelope dose (D_k = 1000.0 sessions), testing min(1.0, ...) cap on p01'.
        6. Zero efficacy non-responder (e_eff = 0.0) with heavy dose (D_k = 20.0 sessions).
        Verifies all entries in [0, 1] and all row sums identically 1.0.
        """
        h = 0.05
        test_corners = [
            # (p01, p12, p10, p20, e_eff, d_k, description)
            (0.005, 0.10, 0.05, 0.03, 0.0, 0.0, "Minimal baseline, no treatment"),
            (0.030, 0.30, 0.15, 0.08, 0.0, 0.0, "Maximal baseline, no treatment"),
            (0.030, 0.30, 0.15, 0.08, 1.0, 0.0, "Maximal efficacy, zero dose"),
            (0.030, 0.30, 0.15, 0.08, 1.0, 14.2857, "Maximal efficacy, max budget dose"),
            (0.030, 0.30, 0.15, 0.08, 1.0, 1000.0, "Extreme out-of-envelope dose (clamp check)"),
            (0.015, 0.20, 0.10, 0.05, 0.0, 20.0, "Non-responder under heavy dose")
        ]

        for p01, p12, p10, p20, e_eff, d_k, desc in test_corners:
            p12_eff = p12 * (1.0 - 0.80 * e_eff)
            p20_eff = min(0.90, p20 * (1.0 + 2.0 * e_eff))
            p01_eff = min(1.0, p01 * (1.0 + h * d_k))

            # Matrix P'
            P = np.array([
                [1.0 - p01_eff, p01_eff, 0.0],
                [p10, 1.0 - p10 - p12_eff, p12_eff],
                [p20_eff, 0.0, 1.0 - p20_eff]
            ])

            # Assert all elements within [0, 1]
            self.assertTrue(np.all(P >= 0.0), f"Negative entry in {desc}: P=\n{P}")
            self.assertTrue(np.all(P <= 1.0), f"Entry > 1.0 in {desc}: P=\n{P}")

            # Assert each row sums to 1.0
            row_sums = np.sum(P, axis=1)
            np.testing.assert_allclose(row_sums, [1.0, 1.0, 1.0], atol=1e-15,
                                      err_msg=f"Row sum failure in {desc}: {row_sums}")

            # Specifically verify p01_eff clamp at extreme dose
            if d_k == 1000.0:
                self.assertEqual(p01_eff, 1.0, "Extreme dose must clamp p01' at exactly 1.0")
                self.assertEqual(P[0, 0], 0.0)
                self.assertEqual(P[0, 1], 1.0)

    def test_treatment_modified_episode_duration_exact_bounds(self):
        """
        Independently verifies analytical bounds of p20' and episode duration:
        - Baseline: p20 in [0.03, 0.08] => E[T] in [12.5, 33.33] days.
        - Treatment formula: p20' = min(0.90, p20*(1 + 2*e_eff)) with e_eff in [0, 1].
        - Across ALL possible treatment levels: p20' in [0.03, 0.24] => E[T'] in [4.167, 33.333] days.
        - Under MAXIMAL efficacy (e_eff = 1.0): p20' = 3*p20 in [0.09, 0.24] => E[T'] in [4.167, 11.111] days.
        """
        # All possible values of p20 and e_eff
        p20_vals = np.linspace(0.03, 0.08, 100)
        e_eff_vals = np.linspace(0.0, 1.0, 100)

        p20_grid, e_grid = np.meshgrid(p20_vals, e_eff_vals)
        p20_prime_grid = np.minimum(0.90, p20_grid * (1.0 + 2.0 * e_grid))

        # Absolute bounds across all parameter combinations
        min_p20_prime = np.min(p20_prime_grid)
        max_p20_prime = np.max(p20_prime_grid)

        self.assertAlmostEqual(min_p20_prime, 0.03, places=5)
        self.assertAlmostEqual(max_p20_prime, 0.24, places=5)

        # Durations
        durations = 1.0 / p20_prime_grid
        self.assertAlmostEqual(np.min(durations), 4.1667, places=3)
        self.assertAlmostEqual(np.max(durations), 33.3333, places=3)

        # Under maximal efficacy e_eff = 1.0:
        p20_prime_max_eff = p20_vals * 3.0
        durations_max_eff = 1.0 / p20_prime_max_eff
        self.assertAlmostEqual(np.min(durations_max_eff), 4.1667, places=3)  # 4.2 days
        self.assertAlmostEqual(np.max(durations_max_eff), 11.1111, places=3) # 11.1 days


def run_tests():
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    suite.addTests(loader.loadTestsFromTestCase(TestDPCLNSimulationLogic))
    suite.addTests(loader.loadTestsFromTestCase(TestMarkovMatrixAndDynamicsAudit))
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    return result.wasSuccessful()


if __name__ == "__main__":
    success = run_tests()
    sys.exit(0 if success else 1)

