#!/usr/bin/env Rscript
#
# D4 and D5 futility boundary calibration, using the REAL D3 null futility
# rate from the completed Eq11_Calibration.R run (kappa* = 0.20, D3's null
# P_early_fut = 0.8244667) as the matching target.
#
# No MCMC involved (calibrate_matched_futility_boundary/single_matched_futility_rep
# just simulate a trial and compute one Z at the futility look) -- this can
# run locally, no HPC needed. Cheap enough to use much higher precision than
# the original manuscript's 2,000-dataset calibration.

library(DTEAssurance)

set.seed(1)

n_c <- 600
n_t <- 600
total_events <- 840
futility_IF <- 0.5

recruitment_model <- list(method = "power", period = 24, power = 1)

# --- Target: D3's REAL null futility rate at kappa* = 0.20 (from the
#     completed Eq11_Calibration.R real HPC run) ---
target_null_futility_rate <- 0.8244667

cat(sprintf("Target null futility rate (matching D3 at kappa*=0.20): %.4f\n\n",
            target_null_futility_rate))

# --- Scenarios (same as used throughout: null, S2 delayed, S3 immediate) ---

scenarios <- list(
  null = list(lambda_c = 0.07452199, gamma_c = 1.210833,
              delay_time = 0, post_delay_HR = 1),
  S2   = list(lambda_c = 0.07452199, gamma_c = 1.210833,
              delay_time = 3, post_delay_HR = 0.8),
  S3   = list(lambda_c = 0.07452199, gamma_c = 1.210833,
              delay_time = 0, post_delay_HR = 0.8)
)

n_sims <- 20000   # cheap (no MCMC), so generous precision relative to the
# manuscript's original 2,000-dataset D4 calibration
n_cores <- max(1, parallel::detectCores() - 1)

cat(sprintf("n_sims = %d, n_cores = %d\n\n", n_sims, n_cores))

# --- D4: matched log-rank (LRT) futility design ---

analysis_model_D4 <- list(method = "LRT", alpha = 0.025,
                          alternative_hypothesis = "one.sided")

cat("=== Calibrating D4 (matched log-rank futility design) ===\n")
t0 <- Sys.time()
D4_result <- DTEAssurance:::calibrate_matched_futility_boundary(
  n_c = n_c, n_t = n_t,
  recruitment_model = recruitment_model,
  futility_IF = futility_IF, total_events = total_events,
  analysis_model = analysis_model_D4,
  target_null_futility_rate = target_null_futility_rate,
  scenarios = scenarios,
  n_sims = n_sims, n_cores = n_cores, seed = 1
)
cat(sprintf("Took %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat(sprintf("D4 boundary: stop for futility if Z < %.4f\n", D4_result$boundary))
cat("Futility-triggering rate under each scenario:\n")
print(D4_result$scenario_futility_rates)
cat("\n")

# --- D5: matched modestly-weighted log-rank (MW) futility design ---
# NOTE: t_star = 3 used here provisionally (prior median elicited delay).
# CONFIRM this is the committed value before treating this result as final.

t_star_D5 <- 3

analysis_model_D5 <- list(method = "MW", alpha = 0.025,
                          alternative_hypothesis = "one.sided",
                          t_star = t_star_D5, s_star = NULL)

cat(sprintf("=== Calibrating D5 (matched MW futility design, t_star=%.1f) ===\n", t_star_D5))
t0 <- Sys.time()
D5_result <- DTEAssurance:::calibrate_matched_futility_boundary(
  n_c = n_c, n_t = n_t,
  recruitment_model = recruitment_model,
  futility_IF = futility_IF, total_events = total_events,
  analysis_model = analysis_model_D5,
  target_null_futility_rate = target_null_futility_rate,
  scenarios = scenarios,
  n_sims = n_sims, n_cores = n_cores, seed = 2
)
cat(sprintf("Took %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat(sprintf("D5 boundary: stop for futility if Z < %.4f\n", D5_result$boundary))
cat("Futility-triggering rate under each scenario:\n")
print(D5_result$scenario_futility_rates)
cat("\n")

# --- Comparison ---

cat("=== Comparison ===\n")
cat(sprintf("%-10s %10s %10s\n", "", "D4 (LRT)", "D5 (MW)"))
cat(sprintf("%-10s %10.4f %10.4f\n", "boundary", D4_result$boundary, D5_result$boundary))
for (scen in names(scenarios)) {
  cat(sprintf("%-10s %10.4f %10.4f\n", paste0("fut_", scen),
              D4_result$scenario_futility_rates[scen],
              D5_result$scenario_futility_rates[scen]))
}
cat("\nBoth should show ~matching null futility rates (by construction) but\n")
cat("differ under S2 (delayed): D5's MW-based boundary should trigger futility\n")
cat("LESS often under a genuine delayed effect if MW is behaving as intended\n")
cat("(protecting power under the scenario it's specifically designed for).\n")

saveRDS(list(D4 = D4_result, D5 = D5_result, target = target_null_futility_rate,
             t_star_D5 = t_star_D5),
        "D4_D5_calibration_results.rds")
cat("\nSaved to D4_D5_calibration_results.rds\n")