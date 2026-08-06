# BNP to NT-proBNP Conversion Study

An R analysis for validating and recalibrating the Kasahara BNP-to-NT-proBNP conversion formula in a perioperative cohort, with secondary analyses examining associations with clinical outcomes.

## Background

B-type natriuretic peptide (BNP) and N-terminal pro-BNP (NT-proBNP) are cardiac biomarkers used for perioperative vascular risk stratification. While prognostic thresholds have been established for NT-proBNP, many laboratories measure BNP instead. This study validates an existing conversion formula (Kasahara et al., 2019) and derives corresponding BNP thresholds for established NT-proBNP risk categories.

**Clinical Trial Registration:** [NCT05352698](https://clinicaltrials.gov/study/NCT05352698)

## Study Objectives

**Primary:** Validate and recalibrate the Kasahara BNP-to-NT-proBNP conversion formula in a perioperative population.

**Secondary:** Examine associations between biomarker categories and the composite outcome of myocardial injury after noncardiac surgery (MINS) or vascular death within 30 days.

## Project Structure

```
BNP conversion/
├── R/
│   ├── 00_config.R           # Configuration parameters
│   ├── 01_libraries.R        # Package dependencies
│   ├── 02_utils.R            # Utility functions and model definitions
│   ├── 03_data_prep.R        # Data loading and cleaning
│   ├── 04_baseline_table.R   # Table 1 generation
│   ├── 05_kasahara_eval.R    # External validation of Kasahara model
│   ├── 06_new_model.R        # Nested model comparison
│   ├── 07_outcomes.R         # Clinical outcome associations
│   ├── 08_report.R           # Report generation
│   └── main.R                # Orchestrator script
├── assets/                   # Output directory (plots, tables, report)
├── checks/
│   └── verify_release.R      # Cohort, event, and release-output checks
├── renv.lock                 # Package version lockfile
└── BNP NTproBNP protocol.md  # Study protocol
```

## Installation

### Prerequisites

- R 4.5.2 (the version recorded in `renv.lock`)
- RStudio (recommended)

### Setup

1. Clone the repository:
   ```bash
   git clone <repository-url>
   cd "BNP conversion"
   ```

2. Restore the exact package versions used in this analysis:
   ```r
   # renv will auto-activate when you open the project
   # If prompted, run:
   renv::restore()
   ```

   This installs all dependencies with the exact versions recorded in `renv.lock`, ensuring reproducibility.

### Private input data

Patient-level CSV, Excel, and Word files are intentionally excluded from version control. Place the current `BNP Data Analysis Jan 9 2025.csv` in the project root before running the analysis.

## Running the Analysis

From the project root directory in R:

```r
source("R/main.R")
```

This executes all analysis scripts in order and produces outputs in the `assets/` directory.

For the current locked analysis dataset, the expected model cohort is 448 participants. The secondary outcome cohort contains 423 participants, including 20 events. After running the analysis, verify the release outputs with:

```r
source("checks/verify_release.R")
```

### Configuration

Key parameters can be modified in `R/00_config.R`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `seed` | 20250409 | Random seed for reproducibility |
| `bootstrap_reps` | 1000 | Number of bootstrap iterations |
| `bootstrap_cores` | 10 | Parallel cores for bootstrapping |
| `binary_threshold` | 300 | NT-proBNP threshold for binary metrics (pg/mL) |
| `class_breaks` | c(100, 200, 1500) | NT-proBNP category thresholds |
| `rcs_knots` | 3 | Knots for restricted cubic splines |

## Outputs

### Report
- `assets/analysis_report.Rmd` — Reproducible methods and results report source
- `assets/analysis_report.docx` — Generated internal Word report (not tracked)

### Figures
| File | Description |
|------|-------------|
| `calibration_plot_kasahara_log.png` | Kasahara model calibration (log scale) |
| `calibration_plot_kasahara_recal_ols_loglog.png` | Recalibrated model calibration |
| `BA_kasahara_*.png` | Bland-Altman plots (relative and absolute error) |
| `outcome_assoc_combined.png` | Biomarker-outcome curves with matched NT-proBNP x-axis ranges and complete confidence bands |

### Data Files (.rds)
| File | Description |
|------|-------------|
| `kasahara_validation_metrics_with_ci.rds` | External validation metrics with bootstrap CIs |
| `kasahara_recal_optimism_corrected.rds` | Optimism-corrected performance metrics |
| `nested_model_comparison.rds` | Likelihood ratio test results |
| `bnp_threshold_sensitivity.rds` | Derived BNP thresholds by patient profile |
| `ntprobnp_cat_summary.rds` | Measured NT-proBNP outcome counts by category |
| `pred_cat_summary.rds` | Predicted NT-proBNP outcome counts by category |

## Statistical Methods

### Primary Analysis

1. **External validation** of the Kasahara conversion formula
2. **Log-log OLS recalibration** to adjust for systematic miscalibration
3. **Bootstrap optimism correction** (Harrell's method) for internal validation
4. **Nested model comparison** via likelihood ratio test to evaluate whether interaction terms improve prediction

### Validation Metrics

- **Continuous:** Pearson r, R², RMSE (all on log₁₀ scale)
- **Calibration:** Intercept (calibration-in-the-large), slope
- **Classification:** Sensitivity, specificity, PPV, NPV at ≥300 pg/mL
- **Agreement:** Weighted kappa (quadratic) across 4 NT-proBNP categories

### Secondary Analysis

Logistic regression models (with Firth correction for sparse data) examining associations between:
- Measured NT-proBNP categories and composite outcome
- Predicted NT-proBNP categories and composite outcome
- Derived BNP categories and composite outcome

## Key Packages

| Package | Purpose |
|---------|---------|
| `rms` | Regression modeling with restricted cubic splines |
| `mcr` | Passing-Bablok regression (sensitivity analysis) |
| `rsample` | Bootstrap resampling |
| `gtsummary` | Publication-ready summary tables |
| `logistf` | Firth logistic regression |
| `nephro` | Creatinine clearance calculation (Cockcroft-Gault) |

See `renv.lock` for complete package versions.

## Data Requirements

Input data should contain the following variables:

| Variable | Description |
|----------|-------------|
| `Age` | Age in years |
| `Sex` | "M" or "F" |
| `Wt` | Weight in kg |
| `Ht` | Height in cm |
| `Cr` | Serum creatinine in µmol/L |
| `Hb` | Hemoglobin in g/L |
| `BNP` | BNP concentration in pg/mL |
| `NTproBNP` | NT-proBNP concentration in pg/mL |
| `AF` | Atrial fibrillation (0/1) |
| `MINS30` | MINS within 30 days (0/1) — for secondary analysis |
| `vascular_death` | Vascular death within 30 days (0/1) — for secondary analysis |

## References

- Kasahara S, et al. Conversion formula from B-type natriuretic peptide to N-terminal proBNP values in patients with cardiovascular diseases. *Int J Cardiol*. 2019;280:184-189.
- Duceppe E, et al. Preoperative N-Terminal Pro-B-Type Natriuretic Peptide and Cardiovascular Events After Noncardiac Surgery. *Ann Intern Med*. 2020;172(12):843.
