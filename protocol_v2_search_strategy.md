# Literature Search Methods Supplement (PRISMA-S Aligned)

**Manuscript Title:** Integrating Digital Phenotyping with Neuromodulation for Personalized Psychiatry: A Structured Narrative Review and Conceptual Hierarchical Framework  
**Target Journal:** *Discover Mental Health* (Springer Nature)  
**Search Governance:** Protocol v2 (Frozen October 2026)  
**Execution Timing:** October 2026 (Initiated on 8 October 2026; complete final export executed on 15 October 2026, 04:37:32 to 04:37:52 IST)  
**Execution Interface:** National Center for Biotechnology Information (NCBI) Entrez E-utilities API (`esearch` and `efetch`)  
**Primary Database:** PubMed / MEDLINE  

---

## 1. Methodological Overview & Search History

### 1.1 Scope and Purpose
This supplement documents the multi-stream literature search architecture developed to inform the structured narrative review and conceptual framework. In accordance with the reporting principles of the **PRISMA-S extension** (Rethlefsen et al., *Syst Rev* 2021;10:39), this document provides full, unedited search strings, platform specifications, execution dates, deduplication logic, and source limitations to ensure complete transparency and algorithmic reproducibility.

### 1.2 Literature Search History and Date Reconciliation
- **Original Submission:** The original submission stated an informal search cutoff of 2005--2024. However, that search was conducted without a frozen protocol or exported decision logs, and the resulting bibliography inadvertently included citations published in 2025 and 2026.
- **Current Revision:** During the present revision, the literature retrieval was systematically reconstructed from the ground up under a frozen protocol (Protocol v2) and executed via automated API scripts in **October 2026** (initiated on 8 October 2026, with complete final export on 15 October 2026).
- **Coverage Window:** The search strategy captured literature from 2010 to October 2026 for digital phenotyping and machine learning (Streams S1 and S3b), and from database inception through October 2026 for neuromodulation trials and neurological technical precedents (Streams S3a and S5).
- **Transparency Statement:** Adopting PRISMA-S search reporting standards does **not** convert this manuscript into a systematic review. The review remains a structured narrative review and conceptual framework.

---

## 2. Information Sources and Platform Details

| Parameter | Specification |
|---|---|
| **Database** | PubMed / MEDLINE |
| **Platform / Provider** | National Center for Biotechnology Information (NCBI), U.S. National Library of Medicine |
| **Interface** | Entrez Programming Utilities (E-utilities) API (`esearch.fcgi`, `efetch.fcgi`) |
| **API Parameters** | `db=pubmed`, `retmode=xml`, `usehistory=y`, `tool=DP_CLN_Review_Protocol_v2` |
| **Search Timing** | October 2026 (Initiated on 8 October 2026; complete final export on 15 October 2026, 04:37 IST) |
| **Language Limit** | English only |
| **Publication Types Excluded at Query** | Errata, author corrections, retractions, book reviews (`EXC_PUB_TYPE`) |

---

## 3. Exact Finalized Boolean Search Strategies

### Stream S1: Digital Phenotyping Validity & Prediction
- **Target Concept:** Passive smartphone sensing, wearable physiological monitoring, symptom tracking, and machine learning prediction in psychiatric cohorts.
- **Exact Executed Query:**
```text
("digital phenotyping" OR "passive sensing" OR "smartphone sensing" OR "mobile sensing" OR "phone sensor*" OR "smartphone sensor*" OR "mobile phone sensor*" OR wearable*) AND (psychiatr* OR depression OR bipolar OR schizophrenia OR anxiety OR OCD OR PTSD) AND (validity OR reliability OR "symptom severity" OR relapse OR prediction OR diagnosis)
```
- **Raw Records Exported:** **2,522**
- **Date Range Applied:** 2010/01/01 to 2026/10/15

### Stream S3a: Open-Loop Psychiatric Neuromodulation Trials
- **Target Concept:** Randomized sham-controlled, comparative, and open-label clinical trials of invasive and non-invasive neuromodulation in psychiatric indications.
- **Exact Executed Query:**
```text
("deep brain stimulation" OR DBS OR "transcranial magnetic stimulation" OR TMS OR rTMS OR "theta burst stimulation" OR TBS OR "transcranial direct current stimulation" OR tDCS OR "transcranial electrical stimulation" OR tES OR "vagus nerve stimulation" OR VNS) AND ("treatment-resistant depression" OR TRD OR "major depressive disorder" OR MDD OR "bipolar disorder" OR schizophrenia OR "obsessive-compulsive disorder" OR OCD) AND ("clinical trial" OR RCT OR randomized OR efficacy OR sham-controlled OR pivotal)
```
- **Raw Records Exported:** **3,401**
- **Date Range Applied:** Inception to 2026/10/15

### Stream S3b: Closed-Loop & Personalized Psychiatric Neuromodulation
- **Target Concept:** State-dependent, adaptive, biomarker-guided, and personalized neuromodulation in psychiatric disorders ($n \ge 1$).
- **Exact Executed Query:**
```text
("closed-loop" OR "adaptive deep brain stimulation" OR "adaptive DBS" OR aDBS OR "responsive neurostimulation" OR "biomarker-guided" OR "personalized stimulation" OR "personalized neuromodulation" OR "personalized DBS" OR "personalized TMS" OR "state-dependent stimulation" OR "closed-loop neuromodulation") AND (depression OR psychiatric OR OCD OR "obsessive-compulsive disorder" OR bipolar OR schizophrenia)
```
- **Raw Records Exported:** **1,197**
- **Date Range Applied:** 2010/01/01 to 2026/10/15

### Stream S5: Technical Engineering Precedents in Neurology
- **Target Concept:** Adaptive DBS in Parkinson's disease and responsive neurostimulation (RNS) in refractory epilepsy, evaluated strictly as technical and control-theoretic precedents.
- **Exact Executed Query:**
```text
("adaptive deep brain stimulation" OR "responsive neurostimulation" OR RNS OR "closed-loop stimulation") AND (Parkinson* OR epilepsy OR seizure)
```
- **Raw Records Exported:** **1,229**
- **Date Range Applied:** Inception to 2026/10/15

---

## 4. Multi-Stream Deduplication & Record Triage Protocol

1. **Raw Records Pool:** Cumulative sum of records exported across Streams S1, S3a, S3b, and S5 = **8,349 records**.
2. **Deduplication Priority:**
   - **Tier 1:** Exact canonical PubMed Identifier (`pmid:<id>`).
   - **Tier 2:** Normalized Digital Object Identifier (`doi:<doi>`), stripping protocol prefixes and lowercasing.
   - **Tier 3:** Normalized alphanumeric title prefix (first 60 alphanumeric characters, whitespace and punctuation stripped).
3. **Cross-Stream Duplicates Identified:** **212 instances**.
4. **Unified Deduplicated Corpus:** **8,137 unique records**.
5. **Two-Stage Record Triage and Manual Adjudication:**
   - **Stage 1 (Initial Automated Keyword Triage):** Evaluated all 8,137 deduplicated records using rule-based keyword matching against Protocol v2 criteria.
     * **Excluded (Initial):** **990 records** (12.17%) across 4 explicit criteria: off-target interventions ($n=719$), preclinical animal models ($n=151$), somatic conditions outside review scope ($n=112$), and non-empirical publication types ($n=8$).
     * **Retained after Initial Triage:** **7,147 candidate records** (87.83%).
   - **Stage 2 (Audited Manual Adjudication):** Forensic manual audit was conducted across 73 flagged borderline cases and keyword anomalies (fully documented in `literature/full_exports/screening_audit_trail.csv`):
     * *False Inclusions Excluded:* **34 records** (e.g., nutritional or pharmacological trials with incidental actigraphy).
     * *False Exclusions Rescued:* **39 records** (e.g., rescuing hyphenated neuromodulation terms such as "theta-burst" and tACS).
     * *Net Adjudication Adjustment:* -34 + 39 = +5 records.
   - **Final Audited Triage Baseline:** Established **7,152 retained candidate records** and **985 excluded records** ($7,152 + 985 = 8,137$).
6. **Study Selection for Evidence Extraction:**
   - The 7,152 retained candidate records were **not** individually reviewed at full text.
   - In accordance with structured narrative review methodology, detailed evidence extraction was prioritized for **23 core clinical trials, benchmark studies, and engineering precedents** (with 21 studies directly synthesized in Tables 2 and 3, and 2 computational control and microcircuit modeling studies informing Section 10).

---

## 5. Benchmark Recall Verification (100% Recall Across All 6 Known Studies)

Prior to triage, 6 landmark benchmark papers specified in the protocol were verified against the deduplicated corpus:
- **Onnela & Rauch (2016)**, *Neuropsychopharmacology*, PMID: 26818126 (Stream S1) $\to$ **RECALLED**
- **Saeb et al. (2015)**, *J Med Internet Res*, PMID: 26180009 (Stream S1) $\to$ **RECALLED**
- **Scangos et al. (2021)**, *Nature Medicine*, PMID: 34608328 (Stream S3b) $\to$ **RECALLED**
- **Little et al. (2013)**, *Annals of Neurology*, PMID: 23852650 (Stream S5) $\to$ **RECALLED**
- **Morrell et al. (2011)**, *Neurology*, PMID: 21917777 (Stream S5) $\to$ **RECALLED**
- **Oehrn et al. (2024)**, *Nature Medicine*, PMID: 39160351 (Stream S5) $\to$ **RECALLED**

**Benchmark Recall Rate:** **6 / 6 (100.0%)**.

---

## 6. Access Limitations & Excluded Databases

In strict compliance with PRISMA-S Item 2 and reporting integrity standards, the following database access boundaries are documented:
1. **Subscription-Gated Databases:** Scopus, Web of Science Core Collection, PsycINFO, and IEEE Xplore were not accessible without institutional subscription credentials and were not searched.
2. **Clinical Trial Registries:** A separate, dedicated database export of ClinicalTrials.gov was not conducted; clinical trial registrations (NCT identifiers) were evaluated only where reported within primary publications.
3. **Preprint / Gray Literature Servers:** Pilot queries on Europe PMC and OpenAlex were assessed during protocol validation but were not subjected to full systematic extraction.
