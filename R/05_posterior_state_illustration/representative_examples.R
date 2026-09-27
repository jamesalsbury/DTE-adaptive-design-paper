#!/usr/bin/env Rscript
#
# Representative posterior state probabilities (Section 4.2).
#
# For one representative interim dataset under each of the three
# data-generating truths used throughout the case study (null, immediate,
# delayed -- same parameters as S1/S2/S3), fits the posterior via
# update_priors() and reports P(Z=k | interim data), directly demonstrating
# that the posterior correctly favours the true underlying state even
# though only interim (IF=0.5) data are available.
#
# Cheap -- 3 MCMC fits, no HPC needed.

library(DTEAssurance)

set.seed(42)

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

# --- Same true data-generating parameters as S1/S2/S3 ---

scenarios <- list(
  "Null (S1 truth)"      = list(delay_time = 0, post_delay_HR = 1),
  "Immediate (S3 truth)" = list(delay_time = 0, post_delay_HR = 0.8),
  "Delayed (S2 truth)"   = list(delay_time = 3, post_delay_HR = 0.8)
)

results <- list()

for (label in names(scenarios)) {
  sc <- scenarios[[label]]
  cat(sprintf("=== %s ===\n", label))
  
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
  
  cat(sprintf("  Interim dataset: %d patients enrolled, %d events observed, t=%.2f months\n",
              nrow(eligible_df), sum(eligible_df$status), t_interim))
  
  posterior_samples <- update_priors(eligible_df,
                                     control_model = control_model,
                                     effect_model  = effect_model,
                                     n_samples     = 1000)
  
  Z_probs <- attr(posterior_samples, "Z_probs")
  converged <- attr(posterior_samples, "converged")
  
  cat(sprintf("  Converged: %s\n", converged))
  cat(sprintf("  P(Z=1, null)      = %.4f\n", Z_probs["P_Z1"]))
  cat(sprintf("  P(Z=2, immediate) = %.4f\n", Z_probs["P_Z2"]))
  cat(sprintf("  P(Z=3, delayed)   = %.4f\n\n", Z_probs["P_Z3"]))
  
  results[[label]] <- list(
    n_enrolled = nrow(eligible_df),
    n_events = sum(eligible_df$status),
    t_interim = t_interim,
    Z_probs = Z_probs,
    converged = converged
  )
}

cat("=== Summary table (for Section 4.2) ===\n\n")
summary_df <- do.call(rbind, lapply(names(results), function(label) {
  r <- results[[label]]
  data.frame(Scenario = label,
             P_Z1_null = r$Z_probs["P_Z1"],
             P_Z2_immediate = r$Z_probs["P_Z2"],
             P_Z3_delayed = r$Z_probs["P_Z3"])
}))
print(summary_df, row.names = FALSE, digits = 4)

saveRDS(results, "representative_posterior_states.rds")
cat("\nSaved to representative_posterior_states.rds\n")