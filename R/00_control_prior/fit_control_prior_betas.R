#!/usr/bin/env Rscript
#
# Beta-Beta approximation to the pooled control-arm Weibull posterior.
#
# This is the previously-missing reproducible step between:
#   (1) the pooling + Weibull MCMC script (already found/exists), which
#       produces posterior samples lambda2sample, gamma2sample, and
#   (2) the control_model Beta-Beta hyperparameters actually used
#       throughout DTEAssurance (t1_Beta_a/b, diff_Beta_a/b).
#
# Run this AFTER sourcing/running the pooling+MCMC script, so that
# lambda2sample and gamma2sample already exist in the environment (or
# load them from wherever they were saved -- see the loading block below
# if you saved them to disk instead of keeping them in-session).
#
# Produces:
#   - the fitted Beta hyperparameters (compare against control_model's
#     current t1_Beta_a/b = 1499.487/1059.113, diff_Beta_a/b =
#     1639.044/8098.961 -- if this script reproduces those numbers
#     closely, that CONFIRMS this is indeed how they were derived)
#   - the validity check: % of simulated (S(t1), Delta) draws that are
#     probabilistically invalid (Delta > S(t1), i.e. S(t2) < 0)
#   - marginal fit comparison (fitted Beta vs original MCMC-derived
#     S(t1)/Delta samples)
#   - the independence check (correlation between S(t1) and Delta in the
#     true MCMC output, which the Beta-Beta construction assumes is ~0)

# --- Landmark times (must match control_model$t1, control_model$t2) ---
t1 <- 8
t2 <- 12

# --- If lambda2sample/gamma2sample aren't already in your session, load
#     them here instead (uncomment and adjust as needed): ---
# saved <- readRDS("path/to/saved_MCMC_output.rds")
# lambda2sample <- saved$lambda2sample
# gamma2sample  <- saved$gamma2sample

if (!exists("lambda2sample") || !exists("gamma2sample")) {
  stop("lambda2sample and/or gamma2sample not found in the environment. ",
       "Run the pooling + Weibull MCMC script first (or load saved samples).")
}

cat(sprintf("Posterior sample size: %d\n", length(lambda2sample)))
cat(sprintf("Landmark times: t1 = %d, t2 = %d\n\n", t1, t2))

# --- Derive S(t1), S(t2), Delta from each posterior draw ---
# Weibull survival: S(t) = exp(-(lambda * t)^gamma)
# (Note: confirm this matches the parameterisation used in the MCMC model
#  string -- the JAGS model here uses (lambda2 * datTimes[i])^gamma2 in
#  the hazard/log-likelihood, consistent with S(t) = exp(-(lambda*t)^gamma).)

S_t1 <- exp(-(lambda2sample * t1)^gamma2sample)
S_t2 <- exp(-(lambda2sample * t2)^gamma2sample)
Delta <- S_t1 - S_t2

cat("=== Posterior summary of derived quantities ===\n")
cat(sprintf("S(t1): mean = %.4f, sd = %.4f\n", mean(S_t1), sd(S_t1)))
cat(sprintf("S(t2): mean = %.4f, sd = %.4f\n", mean(S_t2), sd(S_t2)))
cat(sprintf("Delta: mean = %.4f, sd = %.4f\n\n", mean(Delta), sd(Delta)))

# --- Fit independent Beta distributions via method of moments ---
# For a Beta(alpha, beta): mean m = alpha/(alpha+beta),
# variance v = alpha*beta / ((alpha+beta)^2 * (alpha+beta+1)).
# Solving: n = m(1-m)/v - 1 ; alpha = m*n ; beta = (1-m)*n

fit_beta_mom <- function(x) {
  m <- mean(x)
  v <- var(x)
  n <- m * (1 - m) / v - 1
  alpha <- m * n
  beta <- (1 - m) * n
  c(alpha = alpha, beta = beta)
}

fit_S1 <- fit_beta_mom(S_t1)
fit_Delta <- fit_beta_mom(Delta)

cat("=== Fitted Beta hyperparameters (method of moments) ===\n")
cat(sprintf("S(t1) ~ Beta(%.3f, %.3f)\n", fit_S1["alpha"], fit_S1["beta"]))
cat(sprintf("Delta ~ Beta(%.3f, %.3f)\n\n", fit_Delta["alpha"], fit_Delta["beta"]))

cat("Compare against control_model's current values:\n")
cat("  t1_Beta_a/b   = 1499.487 / 1059.113\n")
cat("  diff_Beta_a/b = 1639.044 / 8098.961\n")
cat("(Close agreement confirms this is how those numbers were derived.)\n\n")

# --- Validity check: simulate from the fitted (independent) Betas and
#     check the probability of an invalid draw (Delta* > S(t1)*, i.e.
#     S(t2)* < 0) ---

n_check <- 1e6
S1_star <- rbeta(n_check, fit_S1["alpha"], fit_S1["beta"])
Delta_star <- rbeta(n_check, fit_Delta["alpha"], fit_Delta["beta"])
S2_star <- S1_star - Delta_star

invalid_rate <- mean(S2_star < 0)
gt_rate <- mean(Delta_star > S1_star)  # should be identical to invalid_rate

cat("=== Validity check (1,000,000 simulated draws) ===\n")
cat(sprintf("P(S(t2)* < 0)        = %.5f%%\n", 100 * invalid_rate))
cat(sprintf("P(Delta* > S(t1)*)   = %.5f%%  (should match the above)\n\n", 100 * gt_rate))

if (invalid_rate < 0.001) {
  cat("Invalid-draw rate is negligible (<0.001%) -- the independent-Beta\n")
  cat("approximation respects the survival-probability constraints in\n")
  cat("practice, consistent with the tight-posterior argument.\n\n")
} else {
  cat("WARNING: invalid-draw rate is non-negligible. Consider the\n")
  cat("constrained S(t1)*U reparameterisation discussed previously, or\n")
  cat("truncating/rejecting invalid simulated draws downstream.\n\n")
}

# --- Marginal fit check: fitted Beta vs original MCMC-derived samples ---

cat("=== Marginal fit check (Kolmogorov-Smirnov) ===\n")
ks_S1 <- ks.test(S_t1, "pbeta", fit_S1["alpha"], fit_S1["beta"])
ks_Delta <- ks.test(Delta, "pbeta", fit_Delta["alpha"], fit_Delta["beta"])
cat(sprintf("S(t1): D = %.4f, p = %.4f\n", ks_S1$statistic, ks_S1$p.value))
cat(sprintf("Delta: D = %.4f, p = %.4f\n\n", ks_Delta$statistic, ks_Delta$p.value))

# --- Independence check: correlation between S(t1) and Delta in the
#     TRUE MCMC output (the Beta-Beta construction assumes independence) ---

true_cor <- cor(S_t1, Delta)
cat("=== Independence check ===\n")
cat(sprintf("Correlation(S(t1), Delta) in the true MCMC posterior = %.4f\n", true_cor))
cat("(Close to 0 supports the independence assumption used in the\n")
cat(" Beta-Beta construction; the assumption matters less the tighter\n")
cat(" the posterior is, per the tight-posterior argument.)\n\n")

# --- Save everything for the manuscript/response letter ---

validity_summary <- list(
  t1 = t1, t2 = t2,
  n_posterior = length(lambda2sample),
  fit_S1 = fit_S1,
  fit_Delta = fit_Delta,
  invalid_draw_rate = invalid_rate,
  ks_S1 = ks_S1,
  ks_Delta = ks_Delta,
  correlation_S1_Delta = true_cor
)

saveRDS(validity_summary, "control_prior_beta_validity_check.rds")
cat("Saved full results to control_prior_beta_validity_check.rds\n")