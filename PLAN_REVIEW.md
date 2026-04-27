# Thorough Review: `gamm4` Variance Partitioning Plan

This document provides a critical review of the logic, mathematics, and architecture of the proposed plan to integrate `gamm4` with the `variancePartition` workflow.

## 1. Mathematical Critique: Variance Partitioning

### The Core Problem: Splines as Random Effects
In standard linear mixed models (LMMs), a random effect (like Subject) is modeled as drawn from a normal distribution: $u \sim N(0, \sigma^2_{Subject})$. The variance parameter $\sigma^2_{Subject}$ represents biological variation between individuals. `variancePartition` directly uses these $\sigma^2$ values to calculate the percentage of variance explained.

`gamm4` uses a clever mathematical trick to fit splines (like `s(Age)`): it translates the smoothing penalty into a random effect variance component. Thus, the model contains a variance term for the spline, let's call it $\sigma^2_{Xr}$. 
**Critical Flaw:** $\sigma^2_{Xr}$ does **not** represent biological variance explained by Age. It is inversely related to the smoothing penalty ($\lambda = \sigma^2_{residual} / \sigma^2_{Xr}$). A huge $\sigma^2_{Xr}$ just means the curve is very wiggly (low penalty); it does not mean Age explains a lot of variance in the gene expression. 

### The Solution: "Realized Predictions" (Validated)
The plan proposes using "Realized Predictions" instead of variance components. This is mathematically correct and absolutely mandatory for GAMMs.
- Instead of using $\sigma^2$, we calculate the variance of the actual predicted values in our sample: $Var(\hat{y}_{component})$.
- **Fixed Effects:** $Var(X \hat{\beta})$
- **Random Effects:** $Var(Z \hat{u})$

### Correction Needed: Summing Spline Components
The prototype code in the plan combines the linear and non-linear parts of the spline by adding their variances: `linear_var + nonlinear_var` ($Var(X_{age}\hat{\beta}_{age}) + Var(Z_{age}\hat{u}_{age})$).
**Mathematical Correction:** The variance of a sum is only equal to the sum of the variances if the two components are perfectly orthogonal (covariance = 0). While spline bases are designed to be orthogonal theoretically, in a finite sample, $Cov(X_{age}\hat{\beta}_{age}, Z_{age}\hat{u}_{age})$ is rarely exactly zero. 
- **Refined Logic:** To get the true total variance explained by `Age`, the implementation must calculate the combined prediction first, then take the variance: $Var(X_{age}\hat{\beta}_{age} + Z_{age}\hat{u}_{age})$.

## 2. Mathematical Critique: Hypothesis Testing

### Degrees of Freedom and P-values
The `dream` function in `variancePartition` relies heavily on Satterthwaite or Kenward-Roger approximations to calculate the denominator degrees of freedom (DDF) for accurate p-values in small samples.

**The `gamm4` Reality:** 
- The p-values for smooth terms derived from `summary(fit$gam)` do not use Satterthwaite. They use a Bayesian posterior covariance matrix approach (developed by Simon Wood). This is statistically rigorous and the gold standard for GAMs, but it is fundamentally different from `dream`'s underlying math.
- For purely linear fixed effects (e.g., `Sex`), `summary(fit$gam)` uses standard frequentist approximations, which may be less conservative than `dream`'s Satterthwaite approach in small sample sizes (N < 20).

**Conclusion:** The plan's reliance on `mgcv`'s `summary.gam()` is correct because standard mixed-model DDF approximations cannot be applied to penalized splines. However, users must be aware that p-values for linear covariates might differ slightly from standard `dream` outputs due to the underlying engine change.

## 3. Mathematical Critique: Voom Precision Weights

The plan introduces `voomWithGamm4Weights`. 
`voom` estimates the mean-variance relationship of RNA-seq data to generate precision weights. To do this for mixed models, `voomWithDreamWeights` fits an LMM for *every* gene just to extract the residual variance ($\sigma^2_{res}$) and fitted values.

**Computational Bottleneck:** `gamm4` is computationally demanding. Fitting `gamm4` on 20,000 genes to get weights, and then fitting `gamm4` *again* to get the final statistics (`gamm4_dream`), will be extremely slow.

**Refinement Strategy:** The mean-variance trend in RNA-seq is a global property driven by count sizes, not by subtle non-linear age effects. 
- The implementation of `voomWithGamm4Weights` should optionally allow fitting a faster approximation for the weight estimation step (e.g., using natural cubic splines `ns(Age)` in a standard `lmer` model) to calculate weights rapidly, before passing those weights to the rigorous `gamm4` engine for the final fit.

## 4. Architectural Review

### "Fit then Extract" Workflow (limma style)
The plan proposes moving away from `variancePartition`'s all-in-one approach to a `limma`-style workflow:
1. `gamm4_dream()` -> returns a model list.
2. `gamm4_topTable()` -> extracts p-values.
3. `gamm4_extractVarPart()` -> extracts variance fractions.

**Assessment:** Highly recommended. 
- GAMMs are complex objects. Storing the fitted models allows the user to plot the predicted curves (`plot.gam()`), check diagnostics (`gam.check()`), and extract various statistics without having to refit the computationally expensive models. 
- This architecture provides maximum flexibility and aligns better with standard bioconductor practices for complex modeling.

## Summary of Changes Required for Implementation
1. **Variance Calculation:** Update the math in `gamm4_extractVarPart` to calculate $Var(Prediction_{linear} + Prediction_{nonlinear})$ directly to account for covariance.
2. **Voom Optimization:** Build a "fast-path" approximation into `voomWithGamm4Weights` to prevent doubling the computational runtime.
3. **Documentation:** Explicitly document the divergence from Satterthwaite approximations so users understand the p-value derivation.

---

## 5. Alternative Strategy: Unpenalized Natural Splines (`lme4` + `ns()`)

Given the computational bottleneck of `gamm4` (fitting a penalized REML model per gene), a pragmatic and highly effective alternative is to use unpenalized natural cubic splines (`ns()`) within the standard `lme4` / `variancePartition` framework.

### The Trade-off: Penalized vs. Unpenalized
*   **`gamm4` (Penalized):** Automatically tunes the "wiggliness" (degrees of freedom) for every gene. If a gene is perfectly linear, the penalty shrinks the curve to a straight line (1 df). If it's complex, it allows more df. The cost is massive computational overhead.
*   **`lme4` + `ns()` (Unpenalized):** You pre-specify the degrees of freedom (e.g., `y ~ ns(Age, df=3) + (1|Subject)`). It is extremely fast and integrates seamlessly into `voomWithDreamWeights`.

### Mitigating Overfitting
The primary concern with unpenalized splines is overfitting—forcing a complex curve onto noisy linear data. This is mitigated by:
1.  **Using Natural Splines (`ns`) instead of B-splines (`bs`):** Natural splines are constrained to be strictly linear beyond the boundary knots (the youngest and oldest timepoints). This prevents the wild "tail-wagging" characteristic of standard polynomial regression.
2.  **Restricting Degrees of Freedom:** For longitudinal RNA-seq, we rarely have the temporal resolution to support highly complex waves. Setting `df = 2` or `df = 3` is usually optimal. A `df=2` allows for exactly one "bend" (e.g., a U-shape or an inverted U-shape), capturing most meaningful non-linear biology without providing enough flexibility to overfit to noise.

### Pipeline Adaptation for `ns()`
If adopting the `ns()` approach, the pipeline simplifies drastically:
*   **Variance Partitioning (`extractVarPart`):** Works out-of-the-box. `variancePartition` intelligently groups the columns generated by `ns(Age, df=3)` (e.g., `ns1`, `ns2`, `ns3`) and calculates the variance of their combined linear predictor. You get a single "% variance explained" for the overall Age effect.
*   **Hypothesis Testing (`dream`):** Works out-of-the-box, but requires a joint test. Standard `topTable` will output separate p-values for `ns1`, `ns2`, and `ns3`. To ask "Is Age significant?", we perform an F-test (or Likelihood Ratio Test) comparing the full model (`~ ns(Age) + (1|Subject)`) against the null model (`~ (1|Subject)`), jointly testing the null hypothesis that all spline coefficients are zero.

---

## 6. Computational Estimation for Complex Study Designs

To illustrate the stark computational differences between these approaches, consider the following study design:
*   **Dataset:** 115 Subjects × 6 longitudinal samples = **690 Total Samples**
*   **Transcriptome:** ~15,000 Genes
*   **Model Formula:** `~ (1|Subject.ID) + sex.numeric + (1|RNA.isolation.Batch) + (1|Lib_prep_batches) + (1|Year.Drawn) + [Age Term] + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent + cellfreq1 + cellfreq2 + cellfreq3 + cellfreq4 + cellfreq5`

This is a highly complex model with **4 distinct random intercepts** and **10 continuous/categorical fixed effects** (including the 5 cell frequencies).

### A. Linear Model (`Age.months`)
*   **Mechanics:** Fits standard `lmer` models via `variancePartition::dream` (which uses parallelization). 
*   **Speed:** Very fast. The 10 fixed effects slightly increase the size of the design matrix $X$, but fixed-effect estimation is computationally cheap in `lme4`. ~0.3 – 0.6 seconds per gene per core.
*   **Estimated Total Runtime (8 Cores):** **~10 to 20 minutes**.

### B. Natural Splines (`ns(Age.months, df=3)`)
*   **Mechanics:** Fits standard `lmer` models, adding 2 additional columns to the fixed effects design matrix (total 12 fixed effects). The underlying REML optimization for the 4 random intercepts is practically identical to the linear model.
*   **Speed:** ~0.4 – 0.8 seconds per gene per core.
*   **Estimated Total Runtime (8 Cores):** **~15 to 25 minutes**.

### C. Penalized GAMM (`gamm4` with `s(Age.months, bs='cr')`)
*   **Mechanics:** `gamm4` translates the smoothing penalty into a 5th variance component and iteratively calls `lmer` / `optim` to find the optimal smoothing parameter. Juggling 4 random intercepts *plus* the smoothing penalty across 690 samples creates a highly complex, multidimensional optimization surface.
*   **Speed:** Conservatively 10 to 30 seconds per gene per core (could be much worse for genes that struggle to converge due to the high number of covariates).
*   **Estimated Total Runtime (8 Cores):** **~8 to 15 hours**.

**Verdict:** The `ns(df=3)` approach provides the ability to model non-linear age trajectories with virtually zero computational penalty compared to a standard linear model, whereas the `gamm4` approach turns a quick coffee-break analysis into an overnight compute job.

---

## 7. Rigorous Selection of Degrees of Freedom and Knot Placement

When moving away from penalized splines (`gamm4`), which auto-tune complexity, the burden falls on the researcher to rigorously define the flexibility of the spline model.

### A. The Flaw of Arbitrary Linear Splines
In previous scripts, **Linear Splines** (piecewise linear models) were used via `lspline::elspline(Age, n = 8)`. This essentially cuts the age range into 8 equal bins and fits straight lines connecting them.

**Pros of Linear Splines:**
*   Coefficients have direct interpretation ("The rate of change per month between age 20 and 30 is X").

**Cons (Why Reviewers Dislike Them):**
*   **Arbitrary Knots:** The choice of 8 knots is arbitrary. Why not 7? Why not 9? If the knots do not align with true biological transitions (like puberty at age 12), the model will perform poorly.
*   **Non-Differentiable (Sharp Edges):** Linear splines create harsh "corners" at every knot. Biology rarely changes in sharp geometric angles.
*   **Massive Overfitting Risk:** Using `n = 8` consumes 8 degrees of freedom. In a dataset with only 6 timepoints per subject, 8 degrees of freedom will wildly overfit to random noise or batch effects, destroying statistical power.

### B. Natural Cubic Splines (`ns`)
**Natural Cubic Splines** solve the "sharp edges" problem by fitting smooth, cubic polynomials between knots, and forcing the curve to be strictly linear at the outer edges to prevent tail-wagging.

**Pros over Linear Splines:**
*   Biologically realistic (smooth curves).
*   Safer at the boundaries (doesn't wag wildly for the youngest/oldest subjects).
*   Instead of picking specific "knot locations" based on domain knowledge, you just pick the **degrees of freedom (`df`)**, and the `ns()` function mathematically optimally places the internal knots at quantiles of your data.

### C. Rigorous `df` Selection: The Global Empirical Prior
To rigorously justify the choice of `df` to reviewers, you must prove that your choice is optimal for the dataset as a whole, rather than arbitrary. 

Because fitting and tuning 15,000 independent models via Cross-Validation or AIC/BIC is computationally ruinous, the standard practice in genomics is to establish a **Global Empirical Prior**.

**The Workflow:**
1.  **Subset:** Randomly sample 500 highly variable genes from the expression matrix.
2.  **The AIC/BIC Tournament:** For each of those 500 genes, fit a mixed model using `df=1` (linear), `df=2`, `df=3`, and `df=4`.
3.  **Score:** Extract the Akaike Information Criterion (AIC) or Bayesian Information Criterion (BIC) for each model. BIC is generally preferred as it penalizes complex models more harshly, protecting against overfitting.
4.  **Aggregate & Justify:** Calculate the percentage of genes that "voted" for each `df` (e.g., "75% of genes minimized BIC at df=2").
5.  **Apply Globally:** Use the winning `df` for the entire 15,000-gene differential expression run.

**Why this satisfies reviewers:**
Instead of saying "We picked df=3 because it looked nice," you can state: *"To prevent overfitting while capturing non-linear biological trajectories, we empirically determined the optimal degrees of freedom by evaluating BIC across a random subset of 500 highly variable genes. We found that a Natural Cubic Spline with df=2 minimized BIC for the majority of genes, and this global prior was applied to the full differential expression analysis."*

**The Biological Ceiling Rule & Overlapping Cohorts:** In a simple longitudinal study, a hard mathematical ceiling exists: the maximum `df` should not exceed `(Number of unique timepoints per subject) - 2`. However, your study utilizes an **overlapping cohort design** (e.g., cohorts aged 1-6, 2-7, 4-9, 9-14). 

Because the mixed model estimates the spline at the *population level* while accounting for individual baselines via random intercepts (`1|Subject.ID`), the model pools information across all subjects. Therefore, the effective range of the spline is the full 14 years, densely sampled across the different cohorts. 

This means you are **not** strictly limited to `df=2` just because an individual only has 6 timepoints. The population-level curve might easily support `df=3` or `df=4` to capture waves of development across the 14-year lifespan. The `tune_spline_df()` (BIC tournament) approach is especially critical here, as it will mathematically determine exactly how much complexity the pooled cohort data can support without overfitting the individual-level sparsity.