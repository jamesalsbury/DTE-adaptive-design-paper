#!/usr/bin/env Rscript
#
# D4 and D5 futility boundary calibration, using the real kappa*=0.20
# result from the threshold-calibration stage (D3's null futility rate
# = 0.8245) as the matching target.
#
# No MCMC involved (calibrate_matched_futility_boundary()/
# single_matched_futility_rep() just simulate a trial and compute one Z
# at the futility look) -- runs locally, no HPC needed.

library(DTEAssurance)

set.seed(1)

n_c <- 600
n_t <- 600
total_events <- 840
futility_IF <- 0.5

recruitment_model <- list(method = "power", period = 24, power = 1)

# --- Target: D3's real null futility rate at kappa* = 0.20 (from
#     01_run_threshold_calibration.R) ---
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

n_sims <- 20000   # cheap (no MCMC), generous precision
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
#
# t_star = 2 months (CONFIRMED, not provisional). Rationale: accounting
# for the elicited probability of an immediate (T=0) rather than delayed
# effect (1-P_DTE=0.3), 42% of the prior mass over all effect-bearing
# scenarios corresponds to a true delay of 2 months or less -- close to
# the mixture median (~2.1 months) -- rather than the conditional median
# delay alone (3 months). Following Magirr and Burman's guidance to set
# t* closer to zero under genuine uncertainty about whether a delay
# exists at all. See manuscript Section 5.3.3 ("Choice of t* for D5").
#
# This value (t_star=2) is what produces D5_boundary_Z = 0.9424, the
# figure used throughout Table 4/5 (01-04_run_scenario_*.R). An earlier
# provisional run at t_star=3 gave a DIFFERENT boundary (0.9393) -- do
# not reuse that value; t_star=2 is the committed, final choice.

t_star_D5 <- 2

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

cat("\nCompare against the values actually used in 01-04_run_scenario_*.R:\n")
cat("  D4_boundary_Z = 0.9386,  D5_boundary_Z = 0.9424 (t_star=2)\n")
cat("If this run's D4/D5 boundaries differ from the above, investigate\n")
cat("before trusting downstream results -- these are the numbers every\n")
cat("Table 4/5 script assumes.\n")

saveRDS(list(D4 = D4_result, D5 = D5_result, target = target_null_futility_rate,
             t_star_D5 = t_star_D5),
        "D4_D5_calibration_results.rds")
cat("\nSaved to D4_D5_calibration_results.rds\n")