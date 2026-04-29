# Variance Partitioning Methods: Standard Parameters vs. Realized Predictions

When fitting mixed-effects models to RNA-seq data, the goal of variance partitioning is to determine what percentage of the total variation in gene expression is driven by each covariate (e.g., Age, Subject, Batch). 

There are two primary mathematical approaches to calculate these percentages. While they often yield similar results in simple models, understanding the conceptual difference between them is critical—especially when modeling complex, non-linear trajectories like splines.

---

## 1. The Standard Approach: Variance Parameters ($\sigma^2$)
*(This is the default method used internally by `variancePartition::calcVarPart` for random effects)*

In a linear mixed model (LMM), random effects are assumed to be drawn from a normal distribution. For example, the effect of `Subject.ID` is modeled as:
$u_{subject} \sim N(0, \sigma^2_{subject})$

### How it works:
1.  **Random Effects:** The algorithm directly extracts the estimated variance parameter ($\sigma^2_{subject}$, $\sigma^2_{batch}$) from the fitted model. 
2.  **Residuals:** It extracts the residual variance ($\sigma^2_{residual}$).
3.  **Fixed Effects:** Because fixed effects (like `Sex` or linear `Age`) do not have a "variance parameter" in the same way, `variancePartition` calculates the variance of the linear predictor for the specific sample: $Var(X \hat{\beta})$.
4.  **Total Variance:** It sums these components: $V_{total} = \sigma^2_{subject} + \sigma^2_{batch} + Var(X \hat{\beta}) + \sigma^2_{residual}$
5.  **Fractions:** It divides each component by $V_{total}$.

### Pros & Cons:
*   **Pro:** Extremely fast computationally because it mostly just reads pre-calculated parameters from the model summary.
*   **Con (The Mismatch):** It mixes two different philosophical scales. It compares the *theoretical population variance* of the random effects ($\sigma^2$) against the *empirical sample variance* of the fixed effects ($Var(X \hat{\beta})$). 
*   **Con (The Sample Size Problem):** $\sigma^2_{subject}$ represents the variance of the infinite population of all possible subjects. If your specific sample happens to randomly draw subjects that are very similar to each other, the model will still estimate a large $\sigma^2_{subject}$ for the population, even if `Subject.ID` doesn't explain much variance in your *actual data*.

---

## 2. The Empirical Approach: Realized Predictions ($Var(\hat{y})$)
*(This is the method we must use to properly evaluate complex splines and GAMMs)*

Instead of looking at the abstract population parameters, this approach looks at the actual predictions the model makes for the specific data points in your experiment.

### How it works:
It calculates the specific fitted value ($\hat{y}$) contributed by each component for every single observation, and then takes the sample variance of those predictions.
1.  **Fixed Effects (e.g., Age):** Generate the prediction based only on Age: $\hat{y}_{age} = X_{age} \hat{\beta}_{age}$. Then calculate the sample variance: $Var(\hat{y}_{age})$.
2.  **Random Effects (e.g., Subject):** Extract the actual Best Linear Unbiased Predictors (BLUPs) for each subject in your dataset, and multiply them by the design matrix: $\hat{y}_{subj} = Z_{subj} \hat{u}_{subj}$. Then calculate the sample variance: $Var(\hat{y}_{subj})$.
3.  **Residuals:** Calculate the variance of the actual model residuals: $Var(\epsilon)$.
4.  **Fractions:** Sum these sample variances and calculate the percentages.

### Pros & Cons:
*   **Pro (Scale Consistency):** It puts fixed effects, random effects, and residuals on the exact same scale: the empirical variance of the actual sample data. 
*   **Pro (Spline Compatibility):** This is the **only** mathematically valid way to partition variance for Generalized Additive Mixed Models (GAMMs). In `gamm4`, the $\sigma^2$ parameter for a spline represents a *smoothing penalty*, not biological variance. You *must* use the realized predictions to see how much variance the curve actually explains.
*   **Con:** Computationally heavier because you have to perform matrix multiplications for every observation across every model component.

---

## 3. Why This Matters for Your Study

In standard models without splines, `variancePartition`'s hybrid approach works wonderfully. 

However, you are dealing with a complex overlapping cohort design (ages 0-14) and considering multiple spline bases (`ns` vs `lspline`) to model non-linear development. 

1. **Spline Aggregation:** When you use `ns(Age, df=3)`, it creates three separate columns in the fixed effects matrix. To find the "Total Variance Explained by Age", you cannot just add their individual variances together (because the spline components covary). You must use the Realized Predictions approach: calculate the combined prediction for the whole spline ($\hat{y}_{age} = X_{ns1}\beta_1 + X_{ns2}\beta_2 + X_{ns3}\beta_3$) and take the variance of that total prediction. *(Note: modern versions of `variancePartition::calcVarPart` actually do this aggregation for fixed-effect splines correctly under the hood!)*
2. **The Subject vs. Age Debate:** Your reviewers are skeptical about the high variance attributed to `Subject.ID`. The standard $\sigma^2$ approach might theoretically inflate the apparent importance of Subject if the population variance is large, even if the specific subjects in your cohort don't vary that wildly. 
3. **The Solution:** By evaluating the **Realized Predictions**, you calculate the variance that `Subject.ID` *actually explained in your specific 690 samples*. If the Realized Prediction variance for Subject remains massively higher than the Realized Prediction variance for Age (even when using AIC to give Age the most flexible `df` possible), you have an irrefutable, empirical mathematical argument for the reviewers.