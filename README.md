# Environmental Influences on Stress Resilience Through Brain Gray Matter

This repository contains the analysis scripts and code used for the manuscript:

**"Environmental Influences on Stress Resilience Through Brain Gray Matter in Early Adolescence: Longitudinal Insights from the ABCD Study"**

**Status:** *Under Review*

---

## 📌 Overview

This study investigates the enviornmental and neurodevelopmental mechanisms enabling resilience in early adolescence. Using longitudinal data from the Adolescent Brain Cognitive Development (ABCD) Study (N = 8,814), we modeled resilience as a dynamic "stressor-reactivity" phenotype.

**Key Analyses in this Repo:**
1.  **Data Preprocessing:** Cleaning, quality control, and merging of ABCD instruments (demographics, CBCL, SLE, and MRI metrics).
2.  **Environmental Factor Analysis:** Reduction of 344 environmental indicators into 69 latent factors.
3.  **Resilience Scoring:** Calculation of residual-based resilience scores (normative modeling of Stress vs. Psychopathology).
4.  **Environmental Predictors:** Linear Mixed Models (LMM) to identify robust environmental predictors.
5.  **Network Analysis:** Gaussian Graphical Models (GGM) to map the ecosystem of risk.
6.  **Neurobiological Correlates:** Investigation of associations between resilience scores and gray matter volume (GMV).
7.  **Longitudinal Mediation:** Causal mediation analysis to test GMV as a mediator between environmental exposures and future resilience.

---

## ⚠️ Data Access & Privacy

**The raw data cannot be shared in this repository.**

The data used in this study were obtained from the [Adolescent Brain Cognitive Development (ABCD) Study](https://abcdstudy.org) and are stored in the NIMH Data Archive (NDA).

* **Data Release:** 5.1 (DOI: 10.15154/z563-zd24)
* **Access:** Researchers must obtain a Data Use Certification (DUC) from the [NIMH Data Archive](https://nda.nih.gov/).

**Instructions for users:**
1.  Download the raw `.csv` files from NDA.
2.  Place them in a local folder named `data/raw/` (this folder is ignored by git).
3.  Run the cleaning scripts to generate the processed dataframes.

---

## 📂 Repository Structure

```text
├── data/                                # (Not included in repo)
│   ├── raw/                             # Raw ABCD files
│   └── processed/                       # Cleaned Rds/CSV files
├── scripts/
│   ├── 01_preprocess_env.R/             # Merging ABCD instruments & QC
│   ├── 02_reduce_env_dim.R/             # Normative modeling (Residuals)
│   ├── 03_generate_srs.R/               # EFA on environmental vars
│   ├── 04_run_univariate_regression.R/  # Mixed-effects models (Discovery/Holdout)
│   ├── 05_run_networkanalysis.R/        # Network analysis (qgraph/bootnet)
│   ├── 06_brain_development.R/          # Brain development in relation to resilience (Discovery/Holdout)
│   └── 07_run_cme.R/                    # Neuroimaging mediation (causal mediation)
├── outputs/                             # (Not included in repo)
│   ├── caches/
│   ├── figures/
│   └── tables/
├── src/                                
│   └── R/                               # Helper functions
├── 2025_resilience_exposome.Rproj       # R project file
└── README.md
```

---

## 💻 System Requirements & Dependencies

The analysis was performed using R (version 4.4.1).

---

## 📄 Citation

If you use this code or the derivation of the Stressor Reactivity Score in your work, please cite:

    Wong, T. Y., et al. (Under Review). Environmental Influences on Stress Resilience Through Brain Gray Matter in Early Adolescence: Longitudinal Insights from the ABCD Study.

---

## ✉️ Contact

For questions regarding the code or methodology, please open an issue in this repository or contact:

Ting-Yat Wong

Email: yatyat0321@gmail.com