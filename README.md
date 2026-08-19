# Recreational Demand and Access Costs in National Forest of Silvânia

Replication package for the Individual Travel Cost Method (ITCM) analysis of recreational demand in the National Forest of Silvânia (FLONA de Silvânia), Goiás, Brazil.

**Manuscript:** *Recreational Demand and Access Costs in National Forest of Silvânia*
**Dossiê:** Geoturismo e o Cerrado Brasileiro: paisagem, patrimônio natural e desenvolvimento territorial
**Authors:** Laura Andreina Matos Márquez, Karine Borges Machado, Philip Teles Soares, Claudiano Carneiro da Cruz Neto

This repository contains the survey data, the full R pipeline, and the diagnostic outputs used to estimate recreational demand and its sensitivity to access costs at the FLONA de Silvânia. The pipeline is designed to run end-to-end from the raw field data to the final model results, so every number reported in the manuscript can be traced back to code.

---

## Repository structure

```
.
├── database.csv                       # Raw field data (157 responses, as collected)
├── Paper_MCV_REVISED.R                # Full analysis pipeline (single entry point)
├── database_final_with_dummies.csv    # Cleaned analytic dataset (output of the pipeline, n = 154)
├── table_models_final.txt             # Final model + diagnostics table (text output)
├── fig_correlation_revised.png        # Pearson correlations: access cost, distance, travel time
├── fig_residuals_vs_fitted.png        # Deviance residuals vs. fitted values (final NB model)
├── fig_qq_residuals.png               # Normal Q-Q plot of deviance residuals
├── fig_cooks_distance.png             # Cook's distance / influence diagnostics
└── README.md
```

## Requirements

- R ≥ 4.3
- Packages: `tidyverse`, `readr`, `MASS`, `AER`, `lmtest`, `sandwich`, `performance`, `car`, `GGally`

```r
install.packages(c("tidyverse", "readr", "MASS", "AER", "lmtest",
                    "sandwich", "performance", "car", "GGally"))
```

## Reproducing the analysis

The entire pipeline — data cleaning, variable construction, model estimation, and diagnostics — runs from a single script, starting from the raw, uncleaned field data:

```r
setwd("path/to/this/repository")
source("Paper_MCV_REVISED.R")
```

Running the script reproduces, in order:

1. **Data cleaning** (documented and auditable — see *Data notes* below).
2. **Variable construction**, including the generalized access-cost formula.
3. **Descriptive statistics** matching Section 3 of the manuscript.
4. **Poisson and Negative Binomial count-data models**, baseline and full specifications.
5. **The full diagnostic battery**: dispersion test, likelihood-ratio test, AIC/BIC, VIF, robust (sandwich) standard errors, Cook's distance, zero-inflation check, and a residual-clustering test by municipality of origin as a spatial-dependence proxy.
6. **Output files**: the cleaned dataset (`database_final_with_dummies.csv`), the diagnostic figures, and the final results table (`table_models_final.txt`).

No manual steps are required between the raw data and the reported results.

## Data notes

`database.csv` is the raw field dataset exactly as exported from data collection (157 responses). It is **not** analysis-ready: `Paper_MCV_REVISED.R` documents, in code, every correction applied before estimation.

- **Corrupted records (2 removed).** Two respondents carried a Unix-timestamp value that had leaked into the `fuel_consumption` field during data export, producing an implausible access-cost outlier (~US$6.2 billion). These records are dropped in Section 2.1 of the script.
- **Zero-distance records (20 corrected).** Twenty respondents (12.8% of the cleaned sample) had `distance` recorded as `0`, traced to a failed name-matching step between accented and non-accented spellings of the same municipality (e.g., *Silvânia* vs. *Silvania*, *Goiânia* vs. *Goiania*) in the original data-construction pipeline. Distances are recovered from validated in-sample values for the same municipality where available, and from official road-distance figures otherwise. See Section 2.2 of the script for the full mapping.
- **Invalid income code (1 removed).** One respondent carried an unclassifiable income code and is dropped.
- **Final analytic sample:** n = 154.

## Access-cost construction

The generalized access-cost variable (`travel_cost_new`) is built from a transparent, auditable formula rather than taken directly from the raw survey fields:

```
TC_i = 2 · d_i · c_km + k · (w_i / h_m) · (2 · d_i / v)
```

| Parameter | Meaning | Value |
|---|---|---|
| `d_i` | One-way distance, respondent *i* (km) | self-reported / corrected (see above) |
| `c_km` | Monetary cost per km (fuel + vehicle wear) | US$0.26/km |
| `k` | Opportunity-cost-of-time coefficient | 0.30 |
| `w_i` | Monthly income, income-class midpoint (US$) | 121 / 485 / 970 / 1,818 |
| `h_m` | Hours per month | 220 |
| `v` | Assumed average travel speed | 60 km/h |

This replaces an earlier construction of the cost variable that used fuel expenditure alone (with no time-opportunity-cost term) and that was structurally collinear with distance and travel time. Full derivation and rationale are documented in Section 2.3 of the manuscript and in the code comments of `Paper_MCV_REVISED.R`.

## Model specification

Recreational demand is estimated with a Negative Binomial count-data model (log link), with a Poisson baseline for comparison:

```
E(Trips_i | X_i) = exp( β₀ + β₁·Gender_i + β₂·Education_i + β₃·Age_i + β₄·Income_i + β₅·ln(TC_i) )
```

Socioeconomic dummies (age, education, income, gender) are rebuilt directly from the raw ordinal survey variables using the reference categories stated in Section 2.3 of the manuscript, rather than using the pre-existing dummy columns in the raw data, which did not consistently match those definitions.

## Key diagnostics (final model, n = 154)

| Diagnostic | Value |
|---|---|
| AIC (Poisson / Negative Binomial) | 539.86 / 537.24 |
| BIC (Poisson / Negative Binomial) | 558.08 / 558.50 |
| LR test, Poisson vs. NB | χ² = 4.61, df = 1, p = 0.032 |
| Dispersion test (Poisson) | z = 1.15, p = 0.125 |
| Dispersion parameter θ (NB) | 10.70 (SE = 5.98) |
| Variance Inflation Factors | 1.03–1.54 |
| McFadden pseudo-R² | 0.064 |
| Cook's distance (max) | 0.25 |
| Residual clustering by municipality | F(10,143) = 2.16, p = 0.024 |

Full derivations and discussion are in Section 2.4 and 3.1 of the manuscript, and in the response letter to reviewers included with the submission.

## Known limitations

- The sample is an on-site (intercept) survey; visitation counts are truncated at 1 by design (Shaw, 1988; Englin & Shonkwiler, 1995). A truncated-NB re-estimation is a natural extension not implemented here.
- Only municipality-level (not point) origin data are available, so the spatial-dependence test reported here is a residual-clustering proxy rather than a formal geostatistical test (e.g., Moran's I).
- Value-of-time parameters (`c_km`, `k`, `v`) are fixed across respondents rather than individually elicited.

## Citation

If you use this dataset or code, please cite the manuscript: [UNDER REVIEW]

> Márquez, L. A. M., Machado, K. B., Soares, P. T., & Cruz Neto, C. C. (2026). Recreational Demand and Access Costs in National Forest of Silvânia. *[Journal / Dossiê Geoturismo e o Cerrado Brasileiro]*.

## License

Data and code are shared for the purpose of scientific replication. Please contact the corresponding author before any commercial use.

## Contact

Claudiano Carneiro da Cruz Neto — cneto@ufrb.edu.br
Federal University of Recôncavo da Bahia (UFRB)
