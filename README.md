# DP-CLN: Integrating Digital Phenotyping with Neuromodulation for Personalized Psychiatry

This repository provides open-science supplementary materials, reproducible literature search protocols, and the mathematical in silico benchmark specification for the manuscript:

> **"Integrating Digital Phenotyping with Neuromodulation for Personalized Psychiatry: A Structured Narrative Review and Conceptual Hierarchical Framework"**  
> *Target Journal:* Discover Mental Health (Springer Nature)  
> *Authors:* Vibha Tiwari, Prasanna Saxena, Akshra Tiwari, Rebakah Geddam, Muhammad Awais, Hemant Ghayvat.

---

## Repository Contents

| File | Description |
| :--- | :--- |
| `dpcln_stress_test.py` | Complete Python implementation of the in silico hierarchical stress-test benchmark specification (Section 6 of the manuscript), featuring multi-state Markov disease dynamics, circadian variation, baseline signal drift, severity-dependent informative missingness, treatment carryover, tolerance decay, and exact burden-matched random comparators. |
| `test_dpcln_simulation_logic.py` | Comprehensive deterministic unit test suite (24 tests across 2 test classes) verifying control policies, clinician safety veto $g(t)$, refractory enforcement (canonical $\tau_{\text{ref}} = 24$ h, alongside exploratory 48 h sensitivity sweeps), 7-day budget caps ($N_{\text{cap}} = 5$), treatment kinetics, bounded threshold modulation $\theta_b(\hat{r})$, complete $3 \times 3$ Markov stochastic validity, and geometric duration bounds. |
| `protocol_v2_search_strategy.md` | Complete PRISMA-S aligned documentation of the four-stream PubMed/MEDLINE search architecture (Protocol v2) executed via NCBI Entrez E-utilities, detailing full Boolean query strings, retrieval dates, and deduplication criteria. |
| `screening_audit_trail.md` | Screening and adjudication logs across the 8,137 unique deduplicated corpus, detailing Stage 1 rule-based automated triage (990 exclusions) and Stage 2 manual adjudication of 73 edge cases (+5 net adjustment $\to$ 7,152 retained candidates). |
| `reference_claim_verification_table.csv` | Evidence extraction table capturing study design, patient sample size, sensing modalities, neuromodulation parameters, validation strategies, and quantitative findings across core benchmark studies. |
| `requirements.txt` | Python package dependencies (`numpy`, `scipy`). |

---

## Important Scientific Scope & Declarations

1. **Benchmark Specification vs. Empirical Runs:** In accordance with reporting integrity standards, **no simulation runs, numerical performance claims, or clinical efficacy findings are reported or claimed in this repository or the accompanying manuscript**. The computational model is provided strictly as a testable benchmark specification to guide future prospective research.
2. **Review Methodology:** This work constitutes a **structured narrative review and conceptual framework**, adhering to SANRA quality principles and PRISMA-S search transparency standards. It is not a systematic review or meta-analysis. Independent dual screening was not performed; single-investigator automated triage and edge-case adjudication are transparently documented in `screening_audit_trail.md`.

---

## Installation & Running the Unit Tests

The test suite validates all control-theoretic gates, mathematical kinetics, and stochastic matrix properties deterministically without requiring external data:

### Prerequisites
- Python 3.9+
- `numpy`
- `scipy`

### Setup and Execution
```bash
# Clone the repository
git clone https://github.com/PRASANNA-THE-PRASANN1/dp-cln-benchmark.git
cd dp-cln-benchmark

# Install dependencies
pip install -r requirements.txt

# Run all 24 deterministic unit tests
python -m unittest -v test_dpcln_simulation_logic.py
```

The 24 tests verify:
- **Control & Safety Governance (17 tests):** Clinician safety veto enforcement (contraindications, confidence thresholds), rolling refractory intervals (canonical $\tau_{\text{ref}} = 24$ h, alongside exploratory 48 h sensitivity sweeps), 7-day weekly dose budgets ($N_{\text{cap}} = 5$), pharmacodynamic saturation ($w_k$), chronic tolerance desensitization ($\text{tol}_k$), bounded risk-adaptive threshold clipping $\theta_b(\hat{r})$, and absence of state leakage.
- **Mathematical Matrix & Kinetic Audits (7 tests):** Exhaustive stochastic validity across 100,000 parameter sweeps (row sums identically 1.0 $\pm 10^{-14}$), 6 deterministic extreme corner cases, dimensional consistency of $h D_k$, non-saturation of over-stimulation risk ($p_{01}' \le 0.051 \ll 1.0$), analytical episode holding times ($p_{20}' \in [0.03, 0.24]$, duration $4.2$--$11.1$ days under maximal efficacy), and dual-clock causal update order.

---

## License

This project is licensed under the [MIT License](LICENSE).
