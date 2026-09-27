#!/usr/bin/env Rscript
#
# Posterior state probability illustration (Section 4.2 / Section 5 / Appendix).
#
# Two parts, sharing the same setup and simulate-and-fit helper:
#
#   Part A: ONE representative interim dataset per scenario (seed=42),
#           reproducing the main-text Table (\label{tab:posterior_states}).
#   Part B: FIVE independent interim datasets per scenario (seeds 1-5, a
#           DIFFERENT seeding scheme from Part A), reproducing the Appendix
#           robustness-check table (\label{tab:posterior_states_multiseed}).
#
# IMPORTANT: Part A and Part B intentionally use different seeds/seeding
# schemes, matching how each table was originally generated. Do not
# "simplify" this to reuse one replicate from Part B as the Part A example
# -- doing so would produce different numbers from those already published
# in the manuscript text.
#
# Cheap (a handful of MCMC fits), no HPC needed.
# Run from this folder:  Rscript 01_representative_and_multiseed_check.R

library(DTEAssurance)

n_c <- 600
n_t <- 600
total_events <- 840
futility_IF <- 0.5
n_events_interim <- ceiling(futility_IF * total_events)

control_model <- list(dist = "Weibull", parameter_mode = "Distribution",
                      t1 = 8, t2 = 12,
                      t1_Beta_a = 1499.487, t1_Beta_b = 1059.113,
                      diff_Beta_a = 1639.044, diff_Beta_b = 8098.961)

effect_model <- list(delay_SHELF = SHELF::fitdist(c(2, 3, 4), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 12),
                     delay_dist = "gamma",
                     HR_SHELF = SHELF::fitdist(c(0.7, 0.8, 0.85), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1),
                     HR_dist = "gamma", P_S = 0.9, P_DTE = 0.7)

recruitment_model <- list(method = "power", period = 24, power = 1)

# Same three data-generating truths as S1/S2/S3 throughout the case study.
scenarios <- list(
  "Null (S1 truth)"      = list(delay_time = 0, post_delay_HR = 1),
  "Immediate (S3 truth)" = list(delay_time = 0, post_delay_HR = 0.8),
  "Delayed (S2 truth)"   = list(delay_time = 3, post_delay_HR = 0.8)
)

# --- Shared helper: simulate one interim dataset under a given truth,
#     fit via update_priors(), return the posterior state probabilities ---

simulate_and_fit <- function(sc) {
  trial_data <- sim_dte(n_c, n_t, lambda_c = 0.07452199, delay_time = sc$delay_time,
                        post_delay_HR = sc$post_delay_HR, dist = "Weibull", gamma_c = 1.210833)
  trial_data <- add_recruitment_time(trial_data, rec_method = recruitment_model$method,
                                     rec_period = recruitment_model$period,
                                     rec_power = recruitment_model$power)
  
  trial_data <- trial_data[order(trial_data$pseudo_time), ]
  t_interim <- trial_data$pseudo_time[n_events_interim]
  
  eligible_df <- trial_data[trial_data$rec_time <= t_interim, ]
  eligible_df$status <- eligible_df$pseudo_time < t_interim
  eligible_df$survival_time <- ifelse(eligible_df$status, eligible_df$time,
                                      t_interim - eligible_df$rec_time)
  
  posterior_samples <- update_priors(eligible_df,
                                     control_model = control_model,
                                     effect_model  = effect_model,
                                     n_samples     = 1000)
  
  attr(posterior_samples, "Z_probs")
}

# =============================================================================
# PART A: representative example (seed=42) -- reproduces the main-text table
# =============================================================================

cat("=== Part A: Representative example (seed=42) ===\n\n")

set.seed(42)
representative_results <- list()
for (label in names(scenarios)) {
  Z_probs <- simulate_and_fit(scenarios[[label]])
  representative_results[[label]] <- Z_probs
  cat(sprintf("%s: P(Z=1)=%.4f  P(Z=2)=%.4f  P(Z=3)=%.4f\n",
              label, Z_probs["P_Z1"], Z_probs["P_Z2"], Z_probs["P_Z3"]))
}

representative_table <- do.call(rbind, lapply(names(representative_results), function(label) {
  r <- representative_results[[label]]
  data.frame(Scenario = label, P_Z1 = r["P_Z1"], P_Z2 = r["P_Z2"], P_Z3 = r["P_Z3"])
}))

cat("\nCompare against manuscript Table (tab:posterior_states):\n")
cat("  S1 (null):      0.1685 / 0.1705 / 0.6610\n")
cat("  S3 (immediate): 0.0110 / 0.0765 / 0.9125\n")
cat("  S2 (delayed):   0.0010 / 0.2495 / 0.7495\n\n")

# =============================================================================
# PART B: 5-seed robustness check -- reproduces the Appendix table
# =============================================================================

cat("=== Part B: 5-seed robustness check ===\n\n")

n_seeds <- 5
multiseed_results <- data.frame()

for (label in names(scenarios)) {
  cat(sprintf("--- %s ---\n", label))
  for (seed in 1:n_seeds) {
    # Different seeding scheme from Part A -- intentional, matches how the
    # Appendix table was originally generated.
    set.seed(seed * 100 + which(names(scenarios) == label))
    
    Z_probs <- simulate_and_fit(scenarios[[label]])
    
    row <- data.frame(Scenario = label, Seed = seed,
                      P_Z1 = Z_probs["P_Z1"], P_Z2 = Z_probs["P_Z2"], P_Z3 = Z_probs["P_Z3"])
    multiseed_results <- rbind(multiseed_results, row)
    
    cat(sprintf("  seed %d: P(Z1)=%.3f  P(Z2)=%.3f  P(Z3)=%.3f\n",
                seed, Z_probs["P_Z1"], Z_probs["P_Z2"], Z_probs["P_Z3"]))
  }
  cat("\n")
}

cat("=== Mean posterior state probabilities across 5 seeds, per scenario ===\n\n")
means <- aggregate(cbind(P_Z1, P_Z2, P_Z3) ~ Scenario, data = multiseed_results, FUN = mean)
print(means, row.names = FALSE, digits = 3)

cat("\nCompare against manuscript Appendix Table (tab:posterior_states_multiseed) means:\n")
cat("  S1 (null):      0.312 / 0.113 / 0.575\n")
cat("  S3 (immediate): 0.028 / 0.575 / 0.397\n")
cat("  S2 (delayed):   0.069 / 0.270 / 0.662\n")

# --- Save both parts together ---

saveRDS(list(representative_table = representative_table,
             multiseed_results = multiseed_results,
             multiseed_means = means),
        "posterior_state_illustration.rds")
cat("\nSaved to posterior_state_illustration.rds\n")
