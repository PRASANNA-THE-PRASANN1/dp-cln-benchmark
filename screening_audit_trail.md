# Screening Reconciliation Report & Literature Selection Audit

**Manuscript:** *Integrating Digital Phenotyping with Neuromodulation for Personalized Psychiatry: A Structured Narrative Review and Conceptual Hierarchical Framework*  
**Target Journal:** *Discover Mental Health* (Springer Nature)  
**Date:** October 2026  
**Document Status:** Final Verified Audit Trail  

---

## 1. Executive Summary & Purpose

This document provides a complete, forensic reconciliation of the literature search, deduplication, automated triage, manual adjudication, and study selection counts for the revised manuscript submitted to *Discover Mental Health*. 

Previous audit drafts reported two sets of screening figures:
1. An initial triage breakdown of **7,147 retained records** and **990 excluded records**;
2. A post-adjudication breakdown of **7,152 retained records** and **985 excluded records**.

This audit confirms that **both sets of numbers are authentic, empirically documented milestones** within a two-stage screening protocol. The difference represents exactly 73 audited reclassifications recorded in `literature/full_exports/screening_audit_trail.csv` (34 false inclusions reclassified to exclude; 39 false exclusions reclassified to include; net $+5$ retained). Both stages sum perfectly to the identical deduplicated corpus denominator of **8,137 records**.

---

## 2. Quantitative Search & Screening Reconciliation

```
========================================================================================
RAW DATABASE RETRIEVAL (PubMed/MEDLINE via NCBI E-utilities, Oct 8-15, 2026)
----------------------------------------------------------------------------------------
Stream S1:  Digital Phenotyping Validity & Relapse Prediction          = 2,522 records
Stream S3a: Open-Loop Psychiatric Neuromodulation Clinical Trials       = 3,401 records
Stream S3b: Closed-Loop & Personalized Psychiatric Neuromodulation      = 1,197 records
Stream S5:  Neurological Technical Precedents (Parkinson's & Epilepsy)  = 1,229 records
----------------------------------------------------------------------------------------
TOTAL RAW EXPORTS                                                      = 8,349 records
========================================================================================
                                       │
                                       ▼
                     [Hierarchical Deduplication Algorithm]
                 (Canonical PMID -> Normalized DOI -> Title Prefix)
                       - Cross-stream duplicates removed: 212
                                       │
                                       ▼
========================================================================================
UNIFIED DEDUPLICATED CORPUS DENOMINATOR                                = 8,137 records
========================================================================================
                                       │
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────────────┐
│ STAGE 1: INITIAL AUTOMATED KEYWORD TRIAGE (Protocol v2 Algorithmic Filter)           │
├──────────────────────────────────────────────────────────────────────────────────────┤
│ Evaluated against formal inclusion/exclusion keyword dictionaries:                   │
│                                                                                      │
│ - Excluded Records (Initial):                                           990 (12.17%) │
│     * EXC_OFF_TARGET_INTERVENTION (No DP sensing / no neuromodulation): 719          │
│     * EXC_ANIMAL (Preclinical animal models lacking human data):        151          │
│     * EXC_OFF_TARGET_CONDITION (Non-psychiatric somatic conditions):    112          │
│     * EXC_PUB_TYPE (Errata, corrections, meeting abstracts):              8          │
│                                                                                      │
│ - Retained Records (Initial In-Scope Pool):                            7,147 (87.83%) │
│                                                                                      │
│ Stage 1 Sum: 7,147 + 990 = 8,137 records                                             │
└──────────────────────────────────────────────────────────────────────────────────────┘
                                       │
                                       ▼
┌──────────────────────────────────────────────────────────────────────────────────────┐
│ STAGE 2: AUDITED MANUAL ADJUDICATION (Screening Audit Trail of Edge Cases)           │
├──────────────────────────────────────────────────────────────────────────────────────┤
│ Forensically audited across 73 flagged borderline cases and keyword anomalies        │
│ (Empirical source: `literature/full_exports/screening_audit_trail.csv`):             │
│                                                                                      │
│ 1. False Inclusions Corrected to Exclude:                                 34 records │
│    - Nutritional/pharmacological trials with incidental actigraphy (e.g.,            │
│      PeptiSleep trial [PMID 42591960], melatonin supplements, sedative trials)       │
│                                                                                      │
│ 2. False Exclusions Rescued to Include:                                   39 records │
│    - Rescued hyphenated compound terms and novel acronyms missed by initial          │
│      exact keyword matching (e.g., "theta-burst", "iTBS", "tACS", "HD-tDCS")         │
│                                                                                      │
│ Net Adjudication Adjustment: -34 + 39 = +5 records                                   │
└──────────────────────────────────────────────────────────────────────────────────────┘
                                       │
                                       ▼
========================================================================================
FINAL AUDITED TRIAGE BASELINE (Post-Adjudication Baseline)
----------------------------------------------------------------------------------------
Post-Adjudication Retained Records:  7,147 + 5                          = 7,152 records
Post-Adjudication Excluded Records:  990 - 5                            =   985 records
----------------------------------------------------------------------------------------
FINAL RECONCILED DENOMINATOR:        7,152 + 985                        = 8,137 records
========================================================================================
```

---

## 3. Disaggregated Stream Distribution

### Post-Triage Stream Breakdown
The 7,152 post-adjudication retained records are distributed across review streams as follows:
- **Stream S1 (Digital Phenotyping & Prediction):** 2,121 records
- **Stream S3a (Open-Loop Psychiatric Trials):** 3,153 records
- **Stream S3b (Closed-Loop & Adaptive Psychiatric Studies):** 878 records (875 initial + 3 rescued)
- **Stream S5 (Neurological Precedents in PD & Epilepsy):** 1,000 records (998 initial + 2 rescued)
- **Total In-Scope Corpus Pool:** **7,152 records**

---

## 4. Evidence Extraction & Benchmark Study Accounting

### 4.1 Narrative Scope vs. Full-Text Screening Declaration
In strict accordance with reporting ethics and international guidelines:
- **No Claim of Comprehensive Full-Text Review:** The 7,152 retained records were **not** all individually retrieved or evaluated at full text. Claiming full-text screening of 7,152 papers by a single investigator would be physically implausible and scientifically fraudulent.
- **Role of the Triage Pool:** The 7,152 records constituted an indexed, keyword-triaged reference corpus from which eligible primary studies and technical benchmarks were systematically identified.

### 4.2 Benchmark Studies Accounting (23 Total Studies)
Detailed full-text data extraction was conducted on **23 core benchmark studies**:
1. **Table 2 (Digital Phenotyping Empirical Evidence; $n=7$ studies):**
   - Saeb et al. (2015) — Community depression sensing ($n=28$; PMC4526997)
   - Chen et al. (2025) — Adolescent MDD nocturnal ECG-HRV ($n=63$)
   - Valenza et al. (2016) — Bipolar mood state transition ($n=14$)
   - Faurholt-Jepsen et al. (2015) — Bipolar smartphone calls/SMS ($n=28$; PMID 26395972)
   - Barnett et al. (2018) — Schizophrenia relapse anomaly detection ($n=17$; PMC6006240)
   - Boukhechba et al. (2018) — Social anxiety mobility entropy ($n=20$)
   - Obuchi et al. (2020) — Smartphone sensing to fMRI functional connectivity ($n=105$)
2. **Table 3 (Neuromodulation Clinical Trials & Engineering Precedents; $n=14$ studies):**
   - *Open-Loop Psychiatric DBS:*
     * Mayberg et al. (2005) — SCC DBS open-label pilot ($n=6$)
     * Holtzheimer et al. (2017) — SCC DBS BROADEN multicenter sham RCT ($n=90$)
     * Dougherty et al. (2015) — VC/VS DBS RECLAIM multicenter sham RCT ($n=30$)
     * Bergfeld et al. (2016) — vALIC DBS crossover trial ($n=25$)
     * Alonso et al. (2015) / Dougherty (2016) — OCD DBS meta-analysis & cohort ($n=268 / 26$)
   - *Non-Invasive Neuromodulation:*
     * O'Reardon et al. (2007) — 10 Hz rTMS sham RCT ($n=301$)
     * Blumberger et al. (2018) — THREE-D iTBS non-inferiority RCT ($n=414$)
     * Cole et al. (2022) — Stanford Neuromodulation Therapy (SNT) RCT ($n=29$)
   - *Closed-Loop & Invasive Psychiatric Biomarkers:*
     * Scangos et al. (2021) — Amygdala gamma-triggered VC/VS DBS ($n=1$; PMC11219029)
     * Sendi et al. (2026) — Intracranial LFP aperiodic exponent in depression ($n=20$)
     * Alagapan et al. (2023) — Chronic SCC LFP sensing of recovery/relapse ($n=10$; PMC10751270)
   - *Neurological Control-Theoretic Precedents:*
     * Little et al. (2013) — STN beta-triggered aDBS in Parkinson's ($n=8$; PMC3886292)
     * Morrell et al. (2011) — ECoG-responsive neurostimulation (RNS) in epilepsy ($n=191$)
     * Oehrn et al. (2024) — Chronic at-home gamma aDBS in Parkinson's ($n=4$; PMC11826929)
3. **Section 10 (Computational Microcircuit & Control Theory Precedents; $n=2$ studies):**
   - Moosavi et al. (2023) — Localized state feedback in nonlinear neural networks
   - Yin et al. (2024 / 2026) — Prefrontal monoaminergic attractor models & closed-loop disruption in OCD

**Total Benchmark Synthesis:** $7 + 14 = 21$ empirical studies in Tables 2 and 3, plus $2$ computational modeling studies in Section 10 = **23 benchmark studies**.

---

## 5. Methodological Justification: PRISMA-S vs. PRISMA 2020 Flow Diagram

- **PRISMA-S Adapted:** The literature search reporting strictly follows PRISMA-S (Rethlefsen et al., 2021) standards: full query strings, platform interfaces (NCBI E-utilities), execution dates (October 8–15, 2026), deduplication algorithms, and single-screener constraints are fully documented.
- **PRISMA 2020 Flow Diagram Rejected:** A traditional PRISMA 2020 flow diagram requires reporting universal full-text retrieval counts and itemized full-text exclusion reasons across all retained records. Because the 7,152 records served as an indexed pool for structured narrative synthesis rather than exhaustive systematic review, inserting a synthetic PRISMA 2020 flow diagram would misrepresent the study design. Narrative reviews are formally guided by SANRA criteria (Baethge et al., 2019), which this manuscript satisfies in full.

---

## 6. Audit Verification Log Tracing

The counts documented above are verified against the following primary files in the repository:
1. `literature/full_exports/title_abstract_screening_summary.json` — Stage 1 counts ($8,137$ total, $7,147$ included, $990$ excluded).
2. `literature/full_exports/screening_audit_trail.csv` — Stage 2 reclassifications ($73$ rows: $34$ to exclude, $39$ to include).
3. `literature/full_exports/title_abstract_screening_decisions_audited.csv` — Stage 2 audited corpus ($8,137$ rows: $7,152$ included, $985$ excluded).
4. `final_submission_package/revised_manuscript_final.tex` — Updated Sections 2.3 & 2.4.
5. `final_submission_package/response_to_reviewers_final.docx` — Updated Comment 2.1 author response.
