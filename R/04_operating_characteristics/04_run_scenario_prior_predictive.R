#!/usr/bin/env Rscript
#
# Table 4/5, Scenario S4 (prior-predictive / assurance-style, using
# FIXED-PROPORTION STRATIFIED sampling across the three latent states --
# Z=1 null, Z=2 immediate, Z=3 delayed -- rather than natural P_S/P_DTE-
# driven sampling. This makes Table 5's stratified reweighting (Section
# 4.6.1) a genuine re-application of the formula to real stratified
# output. Each replicate is tagged with its true state for that purpose.
#
# simulate_trial_stratified() validated in isolation before use here:
# state=1 gives EXACT delay_time=0, HR=1 (0/200 deviations observed);
# state=2/3 correctly sample HR/delay from the elicited priors.
#
# Usage: Rscript 04_run_scenario_prior_predictive.R <n_sims> <seed>
#
# Submitted at paper scale by 04_run_scenario_prior_predictive.sbatch (50 tasks x 2,000 =
# 100,000 replicates). Writes table4_S4_batch_<seed>.rds to the working
# directory; combine with 05_aggregate_table4.R.
# For a quick local test before submitting at full scale, reduce NSIMS to ~5 and array size to 1-2 tasks.

library(DTEAssurance)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript 04_run_scenario_prior_predictive.R <n_sims> <seed>")
n_sims <- as.numeric(args[1])
seed   <- as.numeric(args[2])
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))

# --- Fixed manuscript parameters ---
n_c <- 600; n_t <- 600; total_events <- 840; futility_IF <- 0.5

control_model <- list(dist = "Weibull", parameter_mode = "Distribution",
                      t1 = 8, t2 = 12,
                      t1_Beta_a = 1499.487, t1_Beta_b = 1059.113,
                      diff_Beta_a = 1639.044, diff_Beta_b = 8098.961)

effect_model <- list(delay_SHELF = SHELF::fitdist(c(2, 3, 4), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 12),
                     delay_dist = "gamma",
                     HR_SHELF = SHELF::fitdist(c(0.7, 0.8, 0.85), probs = c(0.25, 0.5, 0.75), lower = 0, upper = 1),
                     HR_dist = "gamma", P_S = 0.9, P_DTE = 0.7)

recruitment_model <- list(method = "power", period = 24, power = 1)

analysis_model_LRT <- list(method = "LRT", alpha = 0.025, alternative_hypothesis = "one.sided")
analysis_model_MW  <- list(method = "MW", alpha = 0.025, alternative_hypothesis = "one.sided",
                           t_star = 2, s_star = NULL)

# --- Finalized calibrated design parameters ---
kappa_star <- 0.20
D4_boundary_Z <- 0.9386
D5_boundary_Z <- 0.9424

# --- Design objects ---
GSD_model_D2 <- list(events = total_events, alpha_spending = c(0.0125, 0.025),
                     alpha_IF = c(0.75, 1), futility_type = "none")
GSD_model_D3 <- list(events = total_events, alpha_spending = c(0.0125, 0.025),
                     alpha_IF = c(0.75, 1), futility_type = "BPP",
                     futility_IF = futility_IF, BPP_threshold = kappa_star)
GSD_model_D4 <- list(events = total_events, alpha_spending = c(0.0125, 0.025),
                     alpha_IF = c(0.75, 1), futility_type = "MatchedZ",
                     futility_IF = futility_IF, futility_boundary_Z = D4_boundary_Z)
GSD_model_D5 <- list(events = total_events, alpha_spending = c(0.0125, 0.025),
                     alpha_IF = c(0.75, 1), futility_type = "MatchedZ",
                     futility_IF = futility_IF, futility_boundary_Z = D5_boundary_Z)
GSD_model_D1 <- list(events = total_events, alpha_spending = c(0.025), alpha_IF = c(1), futility_type = "none")

design_D1    <- DTEAssurance:::make_rpact_design_from_GSD_model(GSD_model_D1)$design
design_D2345 <- DTEAssurance:::make_rpact_design_from_GSD_model(GSD_model_D2)$design

cat(sprintf("Efficacy boundaries: interim(0.75)=%.4f, final=%.4f\n",
            design_D2345$criticalValues[1], design_D2345$criticalValues[2]))
cat(sprintf("D1 final boundary: %.4f\n", design_D1$criticalValues[1]))
cat(sprintf("kappa*=%.2f, D4_Z=%.4f, D5_Z=%.4f\n\n", kappa_star, D4_boundary_Z, D5_boundary_Z))

# --- Stratified S4 trial simulator (validated in isolation: state=1 gives
#     EXACT delay_time=0/HR=1; state=2/3 correctly sample from the priors) ---

simulate_trial_stratified <- function(state, n_c, n_t, control_model,
                                      effect_model, recruitment_model) {
  sampledS1 <- stats::rbeta(1, control_model$t1_Beta_a, control_model$t1_Beta_b)
  sampledDelta <- stats::rbeta(1, control_model$diff_Beta_a, control_model$diff_Beta_b)
  sampledS2 <- sampledS1 - sampledDelta
  solution <- nleqslv::nleqslv(c(10, 1), function(params) {
    lambda <- params[1]; k <- params[2]
    c(exp(-(control_model$t1 / lambda)^k) - sampledS1,
      exp(-(control_model$t2 / lambda)^k) - sampledS2)
  })
  lambda_c_i <- 1 / solution$x[1]
  gamma_c_i <- solution$x[2]
  
  if (state == 1) {
    delay_time <- 0; post_delay_HR <- 1
  } else if (state == 2) {
    delay_time <- 0
    post_delay_HR <- SHELF::sampleFit(effect_model$HR_SHELF, n = 1)[, effect_model$HR_dist]
  } else {
    delay_time <- SHELF::sampleFit(effect_model$delay_SHELF, n = 1)[, effect_model$delay_dist]
    post_delay_HR <- SHELF::sampleFit(effect_model$HR_SHELF, n = 1)[, effect_model$HR_dist]
  }
  
  data <- sim_dte(n_c, n_t, lambda_c_i, delay_time, post_delay_HR,
                  dist = control_model$dist, gamma_c = gamma_c_i)
  data <- add_recruitment_time(data,
                               rec_method = recruitment_model$method,
                               rec_period = recruitment_model$period,
                               rec_power = recruitment_model$power,
                               rec_rate = recruitment_model$rate,
                               rec_duration = recruitment_model$duration)
  data
}

# Fixed-proportion allocator: equal thirds, deterministic cycling.
# True P_S/P_DTE proportions are reintroduced later, in aggregation, by
# reweighting each state's results per Eq. (9) -- NOT by how many
# replicates happen to land in each state here.
state_for_replicate <- function(i) {
  ((i - 1) %% 3) + 1
}

# --- Per-replicate ---
run_one_replicate <- function(i) {
  state <- state_for_replicate(i)
  trial_data <- simulate_trial_stratified(state, n_c, n_t, control_model,
                                          effect_model, recruitment_model)
  
  run_design <- function(design, GSD_model, analysis_model) {
    DTEAssurance:::apply_GSD_to_trial(n_c = n_c, n_t = n_t, trial_data = trial_data,
                                      design = design, total_events = total_events,
                                      GSD_model = GSD_model,
                                      control_model = control_model, effect_model = effect_model,
                                      recruitment_model = recruitment_model,
                                      analysis_model = analysis_model,
                                      update_priors_sims = 1000, n_BPP_sims = 2000)
  }
  
  out_D1 <- run_design(design_D1, GSD_model_D1, analysis_model_LRT)
  out_D2 <- run_design(design_D2345, GSD_model_D2, analysis_model_LRT)
  out_D3 <- run_design(design_D2345, GSD_model_D3, analysis_model_LRT)
  out_D4 <- run_design(design_D2345, GSD_model_D4, analysis_model_LRT)
  out_D5 <- run_design(design_D2345, GSD_model_D5, analysis_model_MW)
  
  extract <- function(out, label) {
    data.frame(design = label, decision = out$decision,
               success = as.numeric(out$decision %in% c("Stop for efficacy", "Successful at final")),
               early_fut = as.numeric(out$decision == "Stop for futility"),
               early_eff = as.numeric(out$decision == "Stop for efficacy"),
               sample_size = out$sample_size, duration = out$stop_time,
               converged = if (!is.null(out$converged)) out$converged else NA)
  }
  
  res <- rbind(extract(out_D1, "D1"), extract(out_D2, "D2"), extract(out_D3, "D3"),
               extract(out_D4, "D4"), extract(out_D5, "D5"))
  res$replicate <- i
  res$state <- state
  res
}

cat(sprintf("Running S4 (n_sims=%d)...\n", n_sims))
t0 <- Sys.time()
if (n_cores > 1) {
  reps <- parallel::mclapply(seq_len(n_sims), run_one_replicate, mc.cores = n_cores)
} else {
  reps <- lapply(seq_len(n_sims), run_one_replicate)
}
results <- do.call(rbind, reps)
cat(sprintf("S4 took %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat(sprintf("State distribution: %s\n", paste(table(results$state[results$design == "D1"]), collapse = " / ")))

settings <- list(scenario = "S4", n_c = n_c, n_t = n_t, total_events = total_events,
                 futility_IF = futility_IF, kappa_star = kappa_star,
                 D4_boundary_Z = D4_boundary_Z, D5_boundary_Z = D5_boundary_Z,
                 t_star_D5 = 2, n_sims = n_sims, n_cores = n_cores, seed = seed,
                 stratified = TRUE,
                 P_S = 0.9, P_DTE = 0.7,   # needed for downstream reweighting
                 package_version = tryCatch(as.character(utils::packageVersion("DTEAssurance")),
                                            error = function(e) NA_character_),
                 timestamp = as.character(Sys.time()))

saveRDS(list(results = results, settings = settings),
        sprintf("table4_S4_batch_%04d.rds", seed))
cat(sprintf("Saved to table4_S4_batch_%04d.rds\n", seed))
