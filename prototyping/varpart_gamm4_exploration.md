---
title: "Variance Partitioning with gamm4 Exploration"
author: "Nicholas Rachmaninoff"
date: "2026-04-27"
output: html_document
---

# 1. Setup and Package Loading


``` r
library(lme4)
library(gamm4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

set.seed(123)
```

# 2. Simulate the Three Gene Scenarios

We simulate 30 subjects, each sampled 8 times over a period of years. We will generate three distinct gene expression profiles to test our variance partitioning logic:

*   **Gene 1 (Linear Age):** A clear, straight-line relationship with Age. Small subject effects.
*   **Gene 2 (Non-linear Age):** A strong sinusoidal wave over Age. Small subject effects.
*   **Gene 3 (High Subject Effect):** Massive variation between subjects, with almost zero effect from Age.


``` r
n_subj <- 30
n_timepoints <- 8

# Create base subjects with random starting ages
subjects <- data.frame(
  Subject = factor(paste0("S", sprintf("%02d", 1:n_subj))),
  Start_Age = runif(n_subj, 20, 65),
  Subj_Eff_Small = rnorm(n_subj, 0, 2),
  Subj_Eff_Large = rnorm(n_subj, 0, 15)
)

# Expand to multiple timepoints
df <- expand.grid(Subject = subjects$Subject, Timepoint = 1:n_timepoints)
df <- left_join(df, subjects, by = "Subject")

# Calculate precise age at each timepoint (assuming ~1.5 years between points)
df$Age <- df$Start_Age + (df$Timepoint * 1.5) + rnorm(nrow(df), 0, 0.2)

# Generate Gene Expression
# Gene 1: Strong linear age effect (slope = 0.8), small subject effect
df$Gene1_Linear <- 0.8 * df$Age + df$Subj_Eff_Small + rnorm(nrow(df), 0, 4)

# Gene 2: Strong sinusoidal age effect, small subject effect
df$Gene2_NonLinear <- 15 * sin(df$Age / 6) + df$Subj_Eff_Small + rnorm(nrow(df), 0, 4)

# Gene 3: Flat age effect (slope = 0.05), massive subject effect
df$Gene3_SubjectHeavy <- 0.05 * df$Age + df$Subj_Eff_Large + rnorm(nrow(df), 0, 4)

# Pivot data for easy plotting
df_long <- df %>%
  pivot_longer(
    cols = starts_with("Gene"),
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
  mutate(Gene = factor(Gene, levels = c("Gene1_Linear", "Gene2_NonLinear", "Gene3_SubjectHeavy")))
```

# 3. Visualize the Simulated Data

First, we look at the Age trajectories for the three genes. Then we look at the Subject boxplots.


``` r
# Plot 1: Age Trajectories
p_age <- ggplot(df_long, aes(x = Age, y = Expression)) +
  geom_point(aes(color = Subject), alpha = 0.4) +
  geom_line(aes(group = Subject, color = Subject), alpha = 0.2) +
  geom_smooth(method = "loess", color = "black", se = FALSE, linewidth = 1.2) +
  facet_wrap(~Gene, scales = "free_y", ncol = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "Age Trajectories by Gene", subtitle = "Black line is Loess smoother")

# Plot 2: Subject Boxplots (Ordered by median expression)
p_subj <- ggplot(df_long, aes(x = reorder(Subject, Expression, median), y = Expression, fill = Subject)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  facet_wrap(~Gene, scales = "free_y", ncol = 1) +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(title = "Variation by Subject", x = "Subjects (Ordered)")

# Combine plots using patchwork
p_age | p_subj
```

```
## `geom_smooth()` using formula = 'y ~ x'
```

![plot of chunk visualize_data](figure/visualize_data-1.png)

# 4. Upgraded Variance Partition Function

In `gamm4`, the smooth term `s(Age)` is decomposed into:
1.  **Linear part:** Included in the fixed effects design matrix `X`. Columns are typically named like `Xs(Age)Fx1`.
2.  **Non-linear part:** Included as a random effect, typically named `Xr`.

This function correctly identifies these components and combines them if requested.


``` r
calc_vp_gamm <- function(model, spline_var = "Age", combine_age = TRUE) {
  
  # Ensure we have the mer object
  if (inherits(model, "gamm4")) {
    mer_obj <- model$mer
  } else {
    mer_obj <- model
  }
  
  # 1. Random Effects Variance
  vc <- as.data.frame(VarCorr(mer_obj))
  re_vars <- numeric()
  for (grp in unique(vc$grp)) {
    if (grp != "Residual") re_vars[grp] <- sum(vc$vcov[vc$grp == grp])
  }
  
  # 2. Fixed Effects Variance
  X <- getME(mer_obj, "X")
  beta <- fixef(mer_obj)
  
  fe_vars <- numeric()
  # We iterate over the columns of X to calculate variance contribution of each term
  # Note: Intercept variance is 0
  for (colname in colnames(X)) {
    if (colname != "X(Intercept)") {
      term_pred <- as.numeric(X[, colname, drop = FALSE] %*% beta[colname])
      fe_vars[colname] <- var(term_pred)
    }
  }
  
  # 3. Residuals
  var_resid <- sigma(mer_obj)^2
  
  # 4. Handle Spline Combination Logic
  # gamm4 names linear fixed effect as Xs(Age)Fx1
  spline_fixed_name <- paste0("Xs(", spline_var, ")Fx1")
  
  linear_var <- fe_vars[spline_fixed_name]
  if (is.na(linear_var)) linear_var <- 0
  
  nonlinear_var <- re_vars["Xr"]
  if (is.na(nonlinear_var)) nonlinear_var <- 0
  
  # Remove the specific age parts from the general lists
  fe_vars <- fe_vars[names(fe_vars) != spline_fixed_name]
  re_vars <- re_vars[names(re_vars) != "Xr"]
  
  # Any other fixed effects (should be none in this simple model)
  other_fe_var <- sum(fe_vars)
  
  # Any other random effects (like Subject)
  other_re_vars <- re_vars
  
  if (combine_age) {
    # Combine linear and non-linear age
    combined_age_var <- linear_var + nonlinear_var
    all_vars <- c(other_re_vars, "Age (Total)" = combined_age_var, "Residuals" = var_resid)
    if (other_fe_var > 0) all_vars["Other (Fixed)"] <- other_fe_var
  } else {
    # Keep separated
    all_vars <- c(other_re_vars, 
                  "Age (Linear)" = linear_var, 
                  "Age (Non-linear)" = nonlinear_var, 
                  "Residuals" = var_resid)
    if (other_fe_var > 0) all_vars["Other (Fixed)"] <- other_fe_var
  }
  
  # 5. Calculate Percentages
  total_var <- sum(all_vars, na.rm = TRUE)
  
  res <- data.frame(
    Component = names(all_vars),
    Variance = unname(all_vars)
  ) %>%
    mutate(Pct_Explained = (Variance / total_var) * 100) %>%
    arrange(desc(Pct_Explained))
  
  return(res)
}
```

# 5. Fit Models and Partition Variance


``` r
genes <- c("Gene1_Linear", "Gene2_NonLinear", "Gene3_SubjectHeavy")
results_list <- list()

for (g in genes) {
  # Fit Model
  fmla <- as.formula(paste(g, "~ s(Age, bs = 'cr')"))
  fit <- gamm4(fmla, random = ~(1|Subject), data = df)
  
  # Partition (Separated)
  vp_sep <- calc_vp_gamm(fit, spline_var = "Age", combine_age = FALSE)
  vp_sep$Gene <- g
  vp_sep$Type <- "Separated"
  
  # Partition (Combined)
  vp_comb <- calc_vp_gamm(fit, spline_var = "Age", combine_age = TRUE)
  vp_comb$Gene <- g
  vp_comb$Type <- "Combined"
  
  results_list[[paste0(g, "_sep")]] <- vp_sep
  results_list[[paste0(g, "_comb")]] <- vp_comb
}

# Bind all results together
df_results <- bind_rows(results_list)
```

# 6. Compare the Results


``` r
# Plot 1: Separated Components
p_sep <- df_results %>%
  filter(Type == "Separated") %>%
  ggplot(aes(x = Gene, y = Pct_Explained, fill = Component)) +
  geom_bar(stat = "identity", position = "stack", color = "black", alpha = 0.9) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  labs(title = "Variance Partitioning: Separated Age Effects", y = "% Variance")

# Plot 2: Combined Components
p_comb <- df_results %>%
  filter(Type == "Combined") %>%
  ggplot(aes(x = Gene, y = Pct_Explained, fill = Component)) +
  geom_bar(stat = "identity", position = "stack", color = "black", alpha = 0.9) +
  scale_fill_manual(values = c("Age (Total)" = "#fc8d62", "Subject" = "#8da0cb", "Residuals" = "#e78ac3")) +
  theme_minimal(base_size = 14) +
  labs(title = "Variance Partitioning: Combined Age Effects", y = "% Variance")

p_sep / p_comb
```

![plot of chunk compare_results](figure/compare_results-1.png)
