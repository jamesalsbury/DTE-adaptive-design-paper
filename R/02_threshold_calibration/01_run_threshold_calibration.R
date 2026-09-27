#!/usr/bin/env Rscript
#
# Threshold calibration (Section 5.3.2) -- real manuscript parameters.
#
# Two stages:
#   Step 0: a quick, no-MCMC check of D2's (efficacy-only) real power under
#           S2 and S3, confirming the pre-specified power_floor has a real
#           margin below D2's ceiling before running the full grid search.
#   Step 1: the Eq. (11) constrained threshold calibration itself -- for
#           each candidate kappa, evaluate power/Type I error/ESS across
#           null/S2/S3, then select kappa* subject to the power floor.
#
# PRE-SPECIFIED before any results are seen (see revision-plan discussion):
#   power_floor    = 0.65 (applied identically to S2 and S3, chosen with
#                    reference to D2's own published S2 baseline of 69.5% --
#                    NOT retuned after seeing output from this script)
#   alt_scenarios  = S2 (delayed), S3 (immediate)
#   lambda_grid    = 0.05 to 0.50, step 0.05 (matches original Table 3 range)
#
# This script contains NO function definitions of its own -- only scenario
# inputs and calls to the package's own (patched, exported) functions.
#
# Usage: Rscript 01_run_threshold_calibration.R <n_sims> <seed>
#   n_sims = number of replicates simulated PER SCENARIO for the main grid
#            search (Step 1). Step 0's D2 check always uses a fixed 5,000
#            replicates, independent of this argument, since it is cheap
#            (no MCMC) and only needs to establish a ballpark margin, not
#            the same precision as the main calibration.
#   seed   = integer seed, used for reproducibility and output file suffix
#
# Submitted at paper scale by 01_run_threshold_calibration.sbatch
# (50 tasks x 300 = 15,000 replicates per scenario); combine the
# eq11_calibration_batch_<seed>.rds outputs with 02_aggregate_and_build_table3.R.
# For a quick local test before submitting at full scale, reduce NSIMS to ~5 and array size to 1-2 tasks.

library(DTEAssurance)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript 01_run_threshold_calibration.R <n_sims> <seed>")
}
n_sims <- as.numeric(args[1])
seed   <- as.numeric(args[2])

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))

# --- Scenario-independent inputs (match the manuscript case study exactly) ---

n_c <- 600
n_t <- 600
total_events <- 840
futility_IF <- 0.5   # the design's own chosen interim timing, confirmed via
# the real PP_Timing HPC run (5,000 replicates/candidate
# IF) and the corresponding Figure 1 histograms -- judgment-
# based per Section 4.5.2 Step 1, not re-optimised here

control_model <- list(dist = "Weibull",
                      parameter_mode = "Distribution",
                      t1 = 8,
                      t2 = 12,
                      t1_Beta_a = 1499.487,
                      t1_Beta_b = 1059.113,
                      diff_Beta_a = 1639.044,
                      diff_Beta_b = 8098.961)

effect_model <- list(delay_SHELF = SHELF::fitdist(c(2, 3, 4), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 12),
                     delay_dist = "gamma",
                     HR_SHELF = SHELF::fitdist(c(0.7, 0.8, 0.85), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1),
                     HR_dist = "gamma",
                     P_S = 0.9,
                     P_DTE = 0.7)

recruitment_model <- list(method = "power",
                          period = 24,
                          power = 1)

analysis_model <- list(method = "LRT",
                       alternative_hypothesis = "one.sided",
                       alpha = 0.025)

# --- Build the true design's future boundaries (efficacy look + final) ---

GSD_model_efficacy <- list(events = total_events,
                           alpha_spending = c(0.0125, 0.025),
                           alpha_IF = c(0.75, 1),
                           futility_type = "none")

rpact_out <- DTEAssurance:::make_rpact_design_from_GSD_model(GSD_model_efficacy)
design    <- rpact_out$design

future_boundaries <- list(
  list(events = ceiling(0.75 * total_events),
       crit   = design$criticalValues[which(abs(design$informationRates - 0.75) < 1e-8)]),
  list(events = total_events,
       crit   = design$criticalValues[which(abs(design$informationRates - 1) < 1e-8)])
)

cat(sprintf("Future boundaries: IF=0.75 -> Z>%.4f ; final -> Z>%.4f\n",
            future_boundaries[[1]]$crit, future_boundaries[[2]]$crit))
cat("Sanity check: these should match the manuscript's 2.241 / 2.047.\n\n")

# --- Three data-generating scenarios (Null / S2 delayed / S3 immediate) ---
# Defined here (rather than later) so both Step 0 and Step 1 can reuse them.

scenarios <- list(
  null = list(lambda_c = 0.07452199, gamma_c = 1.210833,
              delay_time = 0, post_delay_HR = 1),
  S2   = list(lambda_c = 0.07452199, gamma_c = 1.210833,
              delay_time = 3, post_delay_HR = 0.8),
  S3   = list(lambda_c = 0.07452199, gamma_c = 1.210833,
              delay_time = 0, post_delay_HR = 0.8)
)

# --- Pre-specified Eq. (11) settings (defined here so Step 0 can reference
#     power_floor when checking its margin against D2's ceiling) ---

power_floor   <- 0.65
alt_scenarios <- c("S2", "S3")

# =============================================================================
# STEP 0: D2 (efficacy-only) baseline check
#
# Confirms power_floor has a real margin below the ceiling no futility rule
# can exceed, before running the full grid search. No MCMC involved (D2 has
# no futility rule), so this is cheap -- runs in a couple of minutes even at
# n=5,000, independent of n_sims used for the main grid search below.
# =============================================================================

cat("=== Step 0: D2 baseline check (power_floor justification) ===\n\n")

n_D2_check_sims <- 5000

check_D2_power <- function(scenario_truth, label) {
  successes <- logical(n_D2_check_sims)
  
  for (i in seq_len(n_D2_check_sims)) {
    trial_data <- sim_dte(n_c, n_t, lambda_c = scenario_truth$lambda_c,
                          delay_time = scenario_truth$delay_time,
                          post_delay_HR = scenario_truth$post_delay_HR,
                          dist = "Weibull", gamma_c = scenario_truth$gamma_c)
    trial_data <- add_recruitment_time(trial_data, rec_method = recruitment_model$method,
                                       rec_period = recruitment_model$period,
                                       rec_power = recruitment_model$power)
    
    out <- DTEAssurance:::apply_GSD_to_trial(n_c = n_c, n_t = n_t, trial_data = trial_data,
                                             design = design, total_events = total_events,
                                             GSD_model = GSD_model_efficacy,
                                             analysis_model = analysis_model)
    
    successes[i] <- out$decision %in% c("Stop for efficacy", "Successful at final")
  }
  
  phat <- mean(successes)
  se <- sqrt(phat * (1 - phat) / n_D2_check_sims)
  ci <- phat + c(-1.96, 1.96) * se
  
  cat(sprintf("D2 power under %s: %.4f (95%% CI: %.4f - %.4f), n=%d\n",
              label, phat, ci[1], ci[2], n_D2_check_sims))
  
  list(power = phat, se = se, ci = ci)
}

D2_check_S2 <- check_D2_power(scenarios$S2, "S2 (delayed)")
D2_check_S3 <- check_D2_power(scenarios$S3, "S3 (immediate)")

cat(sprintf("\nPre-specified power_floor = %.2f\n", power_floor))
cat(sprintf("Margin below D2/S2 ceiling: %.4f\n", D2_check_S2$power - power_floor))
cat(sprintf("Margin below D2/S3 ceiling: %.4f\n\n", D2_check_S3$power - power_floor))

if (D2_check_S2$power <= power_floor || D2_check_S3$power <= power_floor) {
  warning("Step 0: power_floor is at or above D2's own baseline power under ",
          "one or both alternative scenarios -- the constrained selection ",
          "in Step 1 below is guaranteed to be infeasible. Check power_floor ",
          "before proceeding.")
}

# =============================================================================
# STEP 1: Eq. (11) constrained threshold calibration
# =============================================================================

cat("=== Step 1: Eq. (11) grid search ===\n\n")

update_priors_sims <- 1000
PP_sims <- 2000

grid_runs <- list()
for (scen_name in names(scenarios)) {
  cat(sprintf("Running scenario: %s (n_sims = %d)...\n", scen_name, n_sims))
  t0 <- Sys.time()
  
  grid_runs[[scen_name]] <- run_calibration_grid(
    n_c = n_c, n_t = n_t,
    control_model = control_model, effect_model = effect_model,
    recruitment_model = recruitment_model,
    data_generating_model = scenarios[[scen_name]],
    futility_IF = futility_IF, total_events = total_events,
    future_boundaries = future_boundaries, analysis_model = analysis_model,
    update_priors_sims = update_priors_sims, PP_sims = PP_sims,
    n_sims = n_sims, n_cores = n_cores, seed = seed
  )
  
  cat(sprintf("  %s took %.1f min\n\n", scen_name,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

# --- Summarize each scenario across the pre-specified lambda grid ---

lambda_grid <- seq(0.05, 0.50, by = 0.05)

summaries <- lapply(grid_runs, function(run) DTEAssurance:::summarize_grid_by_lambda(run$raw, lambda_grid))
names(summaries) <- names(grid_runs)

cat("=== Summaries by scenario ===\n")
for (scen_name in names(summaries)) {
  cat(sprintf("\n-- %s --\n", scen_name))
  print(summaries[[scen_name]])
}

# --- Eq. (11) selection: PRE-SPECIFIED floor and scenario set (see header) ---

selection <- DTEAssurance:::select_lambda_star(
  summary_by_scenario = summaries,
  power_floor = power_floor,
  null_scenario = "null",
  alt_scenarios = alt_scenarios
)

cat(sprintf("\n=== Eq. (11) selection (power_floor = %.2f, alt_scenarios = %s) ===\n",
            power_floor, paste(alt_scenarios, collapse = ", ")))
print(selection)

# --- Save everything: Step 0 check, raw replicates, summaries, and the
#     selection result, so every number in the paper can be traced back to
#     this exact run. ---

output <- list(
  D2_baseline_check = list(S2 = D2_check_S2, S3 = D2_check_S3,
                           n_sims = n_D2_check_sims, power_floor = power_floor),
  grid_runs = grid_runs,       # raw per-replicate data (also contains settings/provenance)
  summaries = summaries,       # per-lambda summary per scenario (this IS Table 3, for free)
  lambda_grid = lambda_grid,
  power_floor = power_floor,
  alt_scenarios = alt_scenarios,
  selection = selection,
  seed = seed,
  n_sims = n_sims
)

outname <- sprintf("eq11_calibration_batch_%04d.rds", seed)
saveRDS(output, outname)
