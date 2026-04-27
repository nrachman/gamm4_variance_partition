# Plan: `gamm4` Integration for Variance Partitioning and Hypothesis Testing

## 1. Overview
The goal is to create a "drop-in" suite of R functions that extend the `variancePartition` workflow to support `gamm4`. This will allow for non-linear modeling of covariates (like Age) while still performing rigorous variance partitioning and hypothesis testing across thousands of genes. The workflow will mirror the `limma` approach: **Weight Estimation -> Model Fitting -> Result Extraction**.

## 2. Updated Function Suite

### A. `voomWithGamm4Weights(counts, formula, random_effects, data, ...)`
A `gamm4`-aware version of `voomWithDreamWeights`.
- **Action:** 
  - Converts counts to log2CPM.
  - Fits `gamm4` models to estimate the mean-variance relationship, accounting for smooth terms and random effects.
  - Calculates precision weights based on the residual variance and a loess-fitted trend.
- **Output:** An `EList` object with precision weights, compatible with downstream fitting.

### B. `gamm4_dream(exprObj, formula, random_effects, data, BPPARAM = SerialParam(), ...)`
The core fitting engine, mimicking `limma::lmFit` / `variancePartition::dream`.
- **Features:** 
  - **Parallelization:** Uses `BiocParallel` for efficient processing across genes.
  - **Flexibility:** Supports multiple smooth terms (e.g., `~ s(Age) + s(Time)`) and linear fixed effects.
  - **Weights:** Automatically utilizes precision weights if `exprObj` is an `EList`.
- **Output:** A `gamm4_ModelList` object.

### C. `gamm4_topTable(modelList, coefs = NULL, smooths = NULL, number = Inf)`
Extracts hypothesis testing results, mimicking `limma::topTable`.
- **Action:** 
  - Extracts statistics for specified terms of interest.
  - For **Fixed Effects** (`coefs`): Reports coefficients, t-statistics, and p-values from the `gam` summary.
  - For **Smooth Terms** (`smooths`): Reports Effective Degrees of Freedom (EDF) and p-values.
- **Output:** A single data frame with results for all genes.

### D. `gamm4_extractVarPart(modelList, spline_vars = "Age")`
Calculates variance partitioning, mimicking `variancePartition::extractVarPart`.
- **Logic:** Uses "Realized Predictions" to compare fixed and random effects on the same scale.
- **Spline Handling:** Automatically combines the linear component (fixed) and non-linear component (random) for each spline variable into a single "Total Variance" contribution.
- **Output:** A data frame of variance percentages.

---

## 3. Implementation Logic

### Variance Decomposition
We calculate the variance of the fitted values for each component to ensure comparability:
1. **Fixed Effects:** $Var(X \hat{\beta})$
2. **Random Effects:** $Var(Z \hat{u})$
3. **Splines:** $Var(X_{spline} \hat{\beta}_{spline} + Z_{spline} \hat{u}_{spline})$
4. **Residuals:** $\sigma^2$ (residual variance)

### Hypothesis Testing
We leverage `mgcv::summary.gam()` for robust inference:
- **P-values:** Based on the Bayesian posterior covariance matrix, which accounts for the smoothing penalty.
- **Directionality:** Extracted from the linear component of the spline for global trend assessment.

---

## 4. Integration & Testing Plan

1. **Development:** Create `gamm4_varpart_utils.R`.
2. **Validation:** 
   - Verify `voomWithGamm4Weights` correctly captures the mean-variance trend in simulated count data.
   - Use the 3-gene simulation (Linear, Non-linear, Subject-heavy) to validate `gamm4_extractVarPart`.
   - Ensure `gamm4_topTable` correctly identifies significant smooth effects.
3. **Deployment:** Source the utility script in the user's project directory.
4. **Git:** Commit and push the updated plan and implementation.
