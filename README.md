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
| `test_dpcln_simulation_logic.py` | Deterministic unit test suite (17 tests) verifying control policies, clinician safety veto $g(t)$, refractory enforcement ($\tau_{\text{ref}} = 24$ h), 7-day budget caps ($N_{\text{cap}} = 5$), treatment kinetics, and bounded threshold modulation $\theta_b(\hat{r})$. |
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

The test suite validates all control-theoretic gates and mathematical kinetics deterministically without requiring external data:

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

# Run the 17 deterministic unit tests
python -m unittest test_dpcln_simulation_logic.py
```

All 17 tests verify:
- Clinician safety veto enforcement (contraindications, confidence thresholds).
- Rolling refractory intervals ($\tau_{\text{ref}}$) and 7-day weekly dose budgets ($N_{\text{cap}}$).
- Pharmacodynamic saturation ($w_k$) and chronic tolerance desensitization ($\text{tol}_k$).
- Bounded risk-adaptive threshold clipping: $\theta_b(\hat{r}) = \operatorname{clip}\big(\theta_0 - \kappa(\hat{r} - \bar{r}),\, \theta_{\min},\, \theta_{\max}\big)$.
- Absence of state leakage during uninformative and missingness scenarios.

---

## License

This project is licensed under the [MIT License](LICENSE).
