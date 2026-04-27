# Plan: `gamm4` Integration for Variance Partitioning and Hypothesis Testing

## Clarifying Questions
1. **Parallelization:** Should the `gamm4_dream` function include built-in support for parallel processing (e.g., using `BiocParallel` or `future`), or will you handle the iteration outside the function?

Response: BiocParallel is good.

2. **Multiple Smooths:** Do you anticipate models with multiple smooth terms (e.g., `s(Age) + s(Time)`), or should we optimize for the single-spline case?

Response: If possible, please allow for multiple smooths, but if overly complicated, you can default to a single smooth.

3. **Contrast Handling:** In `dream`, contrasts are often used to compare groups (e.g., `GroupA - GroupB`). For GAMs, hypothesis testing usually focuses on the significance of the smooth term or specific linear coefficients. Do you need a formal contrast matrix interface, or is extracting standard coefficients/p-values sufficient?

Response: Let's leave out the contrast matrix for now and focus on the smooth terms and fixed effects.

4. **Data Scale:** Is the input expression data expected to be log-transformed (e.g., log2CPM) or raw counts? (Standard `dream` logic assumes log-transformed or `voom` data).

Response: The data is expected to be log-transformed (e.g., log2CPM), with VoomWeights, so we can use the `voom` logic from `dream`. You may find it helpful to clone https://github.com/GabrielHoffman/variancePartition/tree/devel/R and review the code. If possible, I would something like VoomWithDreamWeights that allows for smooth terms in the model fitting.

5. **Output Structure:** For the "single table for hypothesis testing," should it include results for all terms in the model for every gene, or just the primary effect of interest?

Response: It should include results for terms of interest. These will be specified ahead of time. I would slightly like to diverge from the logic in variance Partition and focus on developing a first function that fits the models and then a second function that extracts the results (either variance partition or hypothesis tests) in a way that is more akin to the lmFit, eBayes approach of Limma - https://github.com/cran/limma

---

## 1. Overview
The goal is to create a "drop-in" suite of R functions that extend the `variancePartition` workflow to support `gamm4`. This will allow for non-linear modeling of covariates (like Age) while still performing rigorous variance partitioning and hypothesis testing across thousands of genes.

## 2. Proposed Function Suite

### A. `gamm4_dream(exprObj, formula, random_effects, data, ...)`
This is the core engine, mimicking `variancePartition::dream`.
- **Input:** 
  - `exprObj`: Gene expression matrix (genes as rows, samples as columns).
  - `formula`: The fixed effects and smooth terms (e.g., `~ s(Age)`).
  - `random_effects`: The random effects formula (e.g., `~ (1|Subject)`).
  - `data`: Metadata data frame.
- **Action:** 
  - Fits a `gamm4` model for every gene.
  - Returns a `gamm4_ModelList` object (a list of `gamm4` results with custom class attributes).

### B. `gamm4_topTable(modelList, coef = NULL, number = Inf)`
Mimics `limma::topTable` for the GAMM context.
- **Action:** 
  - Iterates through the model list.
  - Extracts p-values, effective degrees of freedom (edf), and coefficients from the `gam` portion of each model.
  - For splines, it will report the `s.table` statistics. For linear terms, it will report `p.table` statistics.
- **Output:** A single data frame containing testing results for all genes.

### C. `gamm4_varPart(modelList, spline_var = "Age", combine_age = TRUE)`
Mimics `variancePartition::extractVarPart`.
- **Action:** 
  - Implements the "Realized Predictions" logic: `var(Xb)` for fixed effects and `var(Zu)` for random effects.
  - Specifically handles `gamm4` splines by combining the linear component (from `mer` fixed effects) and the non-linear component (from `mer` random effects) into a single "Total Variance" for the variable.
- **Output:** A data frame where rows are genes and columns are variance components (percentages).

---

## 3. Implementation Logic (Based on Prototyping)

### Variance Decomposition
To ensure fixed and random effects are comparable, we will calculate the variance of the fitted values for each component:
1. **Fixed Effects:** $Var(X \hat{\beta})$
2. **Random Effects:** $Var(Z \hat{u})$
3. **Splines:** $Var(X_{spline} \hat{\beta}_{spline} + Z_{spline} \hat{u}_{spline})$
4. **Residuals:** $\sigma^2$

### Hypothesis Testing
We will leverage the `summary.gam()` output from `mgcv`. This provides:
- **Approximate p-values** for smooth terms based on the Bayesian posterior covariance matrix.
- **Linear coefficients** for global trend determination (e.g., is the gene overall increasing or decreasing with Age?).

---

## 4. Integration & Testing Plan

1. **Development:** Create a single R script (e.g., `gamm4_varpart_utils.R`) containing the three functions.
2. **Validation:** 
   - Use the simulated data from `prototyping/varpart_gamm4_exploration.Rmd` to ensure the functions recover the expected proportions (e.g., Gene 1 should show high linear Age variance, Gene 2 high non-linear Age variance, Gene 3 high Subject variance).
   - Compare `gamm4_topTable` results against a standard linear `dream` fit for linear-only scenarios.
3. **Deployment:** Provide instructions to source this utility script in:
   `/Users/nicholasrachmaninoff/Library/CloudStorage/GoogleDrive-nrachman.bio@gmail.com/My Drive/nicaragua_paper_include_single_cell/2024_04_02_bulk/Nicaragua_immune_intrinsicness_bulk/scripts/variancePartition/run_varpart`
