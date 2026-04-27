# --- 0. Setup and Data ---
library(lme4)
library(Matrix) # For block matrix operations and sparse matrices
library(lmerTest)

# Toy Data (N=12, 6 donors, 2 conditions)
toy_data <- data.frame(
  DonorID = factor(c("D1", "D2", "D3","D4","D5","D6",
                     "D1", "D2", "D3","D4","D5","D6")),
  Condition = factor(c("PBS", "PBS", "PBS", "PBS", "PBS", "PBS",
                       "IFNG", "IFNG", "IFNG","IFNG", "IFNG", "IFNG"), 
                     levels = c("PBS", "IFNG")),
  Genotype = c(0,0,1,1,2,2,0,0,1,1,2,2),
  Expression = c(10.0, 12.0, 14.0, 11.0, 14.5, 18.0,
                 20, 22, 35, 40, 90,100)
)
Y <- as.matrix(toy_data$Expression) # Outcome Vector (6x1)

#### Step 1: Est. Variance Components ####
# this requires iterative numerical optimization.
# --- 1.1. LMM Fitting (REML) to get estimated parameters ---
fit_lmm <- lmer(Expression ~ Genotype * Condition + (1 | DonorID), data = toy_data)

# --- 1.2. Estimate Donor Variance and Residual Variance
# Estimated Variance Components (NEEDED for V)
var_comp <- as.data.frame(VarCorr(fit_lmm, comp=c("Variance","Std.Dev.")))
var_comp
sigma2_donor <- var_comp[var_comp$grp == "DonorID", "vcov"]
sigma2_epsilon <- var_comp[var_comp$grp == "Residual", "vcov"]

cat("Step 1 Results:\n")
cat("Estimated Donor Variance (sigma2_donor):", sigma2_donor, "\n")
cat("Estimated Residual Variance (sigma2_epsilon):", sigma2_epsilon, "\n\n")

# --- 1.3 Extract X matrix (Fixed Effects Design Matrix)
X <- getME(fit_lmm, "X")
X <- as.matrix(X)
cat("Fixed Effects Design Matrix (X, 12x4):\n")
print(X)
cat("\n")

# --- 1.4 Extract Z matrix (Random Effects Design Matrix)
Z <- getME(fit_lmm, "Z")
Z <- as.matrix(Z) # Convert to dense matrix for V calculation
cat("Random Effects Design Matrix (Z, 12x6):\n")
print(Z)
cat("\n")


#### --- 2. Construct the Variance Matrix (V) --- ####
# V = sigma2_epsilon * I + sigma2_donor * Z %*% t(Z)
V <- sigma2_epsilon * diag(nrow(X)) + sigma2_donor * (Z %*% t(Z))
cat("Step 2 Result: Estimated Variance Matrix (V, 12x12):\n")
print(V)
cat("\n")

##### --- 3. Compute the Inverse of V (V^-1) --- ####
V_inv <- solve(V)
cat("Step 3 Result: Inverse Variance Matrix (V_inv, 12x12):\n")
print(V_inv)
cat("\n")

##### --- 4. Calculate the Inner Term (X' * V^-1 * X) --- ####
XT_Vinv_X <- t(X) %*% V_inv %*% X
cat("Step 4 Result: Inner Precision Matrix (X'V^-1X, 4x4):\n")
print(XT_Vinv_X)
cat("\n")

#### --- 5. Invert the Result to get the Variance-Covariance Matrix (Sigma-hat) --- ####
# Sigma_hat = (X'V^-1X)^-1
Sigma_hat <- solve(XT_Vinv_X)
cat("Step 5 Result: Variance-Covariance Matrix (Sigma_hat, 4x4):\n")
print(Sigma_hat)
cat("\n")
# Note: Sigma_hat should be numerically identical to vcov(fit_lmm)

#### --- 6. Calculate the Right-Hand Side (X' * V^-1 * Y) --- ####
XT_Vinv_Y <- t(X) %*% V_inv %*% Y
cat("Step 6 Result: Weighted Covariance Vector (X'V^-1Y, 4x1):\n")
print(XT_Vinv_Y)
cat("\n")

#### --- 7. Final Solution for Fixed Effects (Beta-hat) --- ####
# Beta_hat = Sigma_hat %*% XT_Vinv_Y
Beta_hat <- Sigma_hat %*% XT_Vinv_Y
cat("Step 7 Result: Estimated Fixed Effects Vector (Beta_hat, 4x1):\n")
print(Beta_hat)
cat("\n")

#### Verify against lmer output ####
cat("--- Verification --- \n")
cat("lmer Coefficients (should match Beta_hat):\n")
print(fixef(fit_lmm))
cat("\n")
cat("lmer V-Cov Matrix (should match Sigma_hat):\n")
print(vcov(fit_lmm))


#### Calculate t-statistics for beta_hat of interaction term ####
# --- 1. Numerator is the beta: last fixed effect (Genotype:ConditionIFNG)
reQTL_estimate <- Beta_hat[4, 1] 
cat("1. Numerator (reQTL Estimate):", reQTL_estimate, "\n")

# --- 2. The Denominator: The Standard Error (SE)
# Variance is the element in the 4th row, 4th column of Sigma_hat
reQTL_variance <- Sigma_hat[4, 4] 
# SE is the square root of the variance
reQTL_SE <- sqrt(reQTL_variance)
cat("2. Denominator (reQTL SE):", reQTL_SE, "\n")

# --- 3. Calculate the t-statistic ---
t_statistic <- reQTL_estimate / reQTL_SE
cat("3. Final t-statistic (Estimate / SE):", t_statistic, "\n")

# --- 4. Calculate p-value of the t-statistic ---
# we need both the t-statistic and the Degrees of Freedom (df), but df is complicated
# here, with mixed effects, df is corrected by the Kenward-Roger method.
# but if we already have fit the lmm, we can get the Kenward-Roger corrected df from there
fit_lmm_KR <- lmerTest::as_lmerModLmerTest(fit_lmm)
anova(fit_lmm_KR, ddf = "Kenward-Roger")
df <- 4  # Degrees of Freedom

# 2. Use the Cumulative Distribution Function (pt)
# pt(q, df) gives the area (probability) to the left of 'q'.
# Since t-tests are typically two-sided, we look at the area in the tails.

# Get the area to the RIGHT of the positive t-statistic (one-sided P-value)
# Area_Right = 1 - Area_Left 
p_value_one_tail <- 1 - pt(t_statistic, df)

# 3. For a TWO-SIDED P-value, multiply by 2 (to cover both the positive and negative tails)
p_value_two_sided <- 2 * p_value_one_tail

cat("Calculated t-statistic:", t_statistic, "\n")
cat("Degrees of Freedom:", df, "\n")
cat("Final Two-Sided P-value:", p_value_two_sided, "\n")

# --- Verification (Compare to lme4's output) ---
lme4_summary <- summary(fit_lmm)
lme4_t_value <- lme4_summary$coefficients["Genotype:ConditionIFNG", "t value"]

cat("\n--- Verification Against lme4 Summary ---\n")
cat("   lme4 Summary t-value:", format(lme4_t_value, digits=4), "\n")
cat("   Calculated t-statistic:", format(t_statistic, digits=4), "\n")
# The values should be numerically identical (up to precision issues)!








