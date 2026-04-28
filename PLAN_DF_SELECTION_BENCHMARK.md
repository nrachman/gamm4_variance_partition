# Plan: Rigorous Spline Selection and Benchmarking

## 1. Objective
The goal of this plan is to establish a rigorous, empirical framework for modeling non-linear age trajectories (ages 0–14) in longitudinal RNA-seq data. 

Specifically, we aim to:
1.  **Benchmark Spline Types:** Compare Natural Cubic Splines (`ns`) against Piecewise Linear Splines (`lspline`) to definitively determine which basis better captures biological reality without overfitting.
2.  **Determine Optimal Complexity:** Rigorously select the optimal degrees of freedom (df) or number of knots ($n \in [1, 10]$) using the **Global Empirical Prior** approach.
3.  **Assess Robustness:** Evaluate how the choice of spline type and complexity impacts the stability of other critical model parameters, specifically the variance explained by Subject (individuality) and the effect sizes of other covariates.
4.  **Ensure Compatibility:** Guarantee that the winning approach integrates seamlessly with the `variancePartition` and `dream` ecosystem for both variance partitioning and hypothesis testing.

---

## 2. Methodology: The Global Empirical Prior Benchmark

Evaluating thousands of models across 15,000 genes is computationally wasteful and statistically noisy. Instead, we will establish a Global Empirical Prior using a representative subset.

### Step 1: Gene Selection
Select a subset of **500 highly variable genes (HVGs)**. These genes are the most likely to exhibit meaningful temporal dynamics and individual variation, providing a strong signal-to-noise ratio for benchmarking.

### Step 2: The Modeling Tournament
For each of the 500 genes, we will fit a suite of linear mixed-effects models using the `variancePartition::dream` framework. The baseline formula will be:
`~ (1|Subject.ID) + [Age Spline] + [Other Covariates]`

We will test two families of splines across 10 levels of complexity:
*   **Family A: Natural Cubic Splines (`splines::ns`)**
    *   Models: `ns(Age, df=1)` through `ns(Age, df=10)`.
    *   *Note:* `df=1` is equivalent to a standard linear model.
*   **Family B: Piecewise Linear Splines (`lspline::elspline`)**
    *   Models: `elspline(Age, n=1)` through `elspline(Age, n=10)`.
    *   *Note:* `n=1` is equivalent to a standard linear model. `n=8` represents the 8 equally-spaced bins used in previous analyses.

*Total Models per Gene: 20*

---

## 3. Evaluation Metrics

To select the winning model, we will evaluate the tournament results across three distinct axes: Statistical Fit, Variance Robustness, and Inference Robustness.

### Metric A: Statistical Fit and Parsimony (AIC / BIC)
We must balance the model's ability to fit the data against the penalty for adding unnecessary complexity.
*   **Akaike Information Criterion (AIC):** Good for predictive models, but tends to favor slightly more complex models.
*   **Bayesian Information Criterion (BIC):** Penalizes extra parameters more heavily than AIC. In high-dimensional omics data with limited sample sizes (6 per subject), BIC is strongly preferred to protect against overfitting.
*   **Benchmark Output:** A bar chart showing the percentage of the 500 genes that achieved their minimum BIC at each `df` (1 through 10) for both `ns` and `lspline`. The peak of this distribution defines our Global Empirical Prior.

### Metric B: Robustness of Subject Identity (Variance and Effect Sizes)
A critical danger of overfitting the Age term is "variance stealing." If a spline becomes too wiggly (e.g., `lspline(n=10)`), it may begin fitting random noise that actually belongs to the subject's unique baseline.
*   **Analysis:** For each gene, we will extract the Percentage of Variance Explained by `Subject.ID` using `variancePartition::extractVarPart`, as well as the actual effect size estimates (BLUPs - Best Linear Unbiased Predictors) for the `Subject.ID` random intercepts.
*   **Benchmark Output:** A line plot tracking the median Subject Variance (across the 500 genes) as a function of spline complexity (df 1-10), and a measure of BLUP stability (e.g., correlation of Subject effects between df=2 and df=10).
*   **Interpretation:** A robust spline method will show stabilization. If Subject Variance drops precipitously as `df` increases, or if the subject-specific BLUPs become highly volatile, the spline is structurally overfitting and absorbing individual-level biological identity.

### Metric C: Robustness of Inference (Age P-values and Covariate Stability)
We need to understand how the choice of spline impacts our statistical power to detect Age effects, and ensure it doesn't destabilize the rest of the model.
*   **Analysis (Age P-values):** We will perform joint F-tests (ANOVA) on the spline terms across all models to extract the overarching p-value for the Age effect. We will track how the number of "significant" age-associated genes changes as complexity increases.
*   **Analysis (Covariate Stability):** We will track the estimated coefficient (log-fold change) and p-value of a stable covariate (e.g., `sex`) across the 20 models.
*   **Benchmark Output:** 
    1. A plot showing the distribution of Age p-values (or the count of significant genes at FDR < 0.05) across the `df` spectrum.
    2. Variance/Standard Deviation of the `sex` log-fold change estimates across the models.
*   **Interpretation:** As `df` increases, power to detect true Age effects may dilute because the test is penalized for the extra degrees of freedom. We want to find the "sweet spot" of maximum power. Additionally, a highly volatile `sex` effect indicates the spline is introducing collinearity or structural instability into the design matrix.

---

## 4. Alternative Methods Considered

While `ns` and `lspline` are the primary candidates due to their compatibility with standard LMMs, other methods exist:

1.  **Fractional Polynomials:** An automated way to find the best fitting polynomial (e.g., $Age^{-1} + Age^{0.5}$). *Cons:* Highly computationally intensive to tune per gene, and the resulting coefficients are biologically uninterpretable. Not natively supported by `variancePartition`.
2.  **Orthogonal Polynomials (`poly(Age, degree=3)`):** Fits smooth global curves. *Cons:* Unlike splines, polynomials are global—a change in expression at age 14 mathematically forces a change in the curve at age 1. Splines are local (knots isolate effects), making them vastly superior for segmented developmental stages.
3.  **Penalized Splines (`gamm4`):** The "gold standard" for preventing overfitting by penalizing wiggliness. *Cons:* As estimated previously, this increases compute time from ~20 minutes to ~10 hours and complicates variance partitioning math.

**Conclusion on Alternatives:** `ns` remains the strongest theoretical candidate because it offers local, smooth biological realism (unlike polynomials) and computational efficiency (unlike `gamm4`), while mathematical bounds (the BIC tournament) prevent the overfitting seen in `lspline`.

---

## 5. Implementation & Deployment

Once the benchmark is complete, the workflow for the main analysis will be identical to the standard `variancePartition` pipeline:

1.  **Run Benchmark:** Execute the 500-gene tournament script.
2.  **Declare Winner:** e.g., "Natural splines with df=3 won the BIC tournament and maintained robust Subject variance."
3.  **Final Fit:**
    ```R
    # The winner drops perfectly into the existing dream ecosystem
    form <- ~ (1|Subject.ID) + ns(Age.months, df=3) + sex + (1|Batch)
    vobj <- voomWithDreamWeights(counts, form, metadata)
    fit  <- dream(vobj, form, metadata)
    ```
4.  **Extract Results:** 
    *   `extractVarPart(fit)` handles `ns` and `lspline` natively, aggregating the terms automatically.
    *   A custom wrapper `extract_spline_pvalues(fit)` will be provided to run the joint F-test (ANOVA) across the multi-column spline term to return a single, clean p-value for the Age effect.
