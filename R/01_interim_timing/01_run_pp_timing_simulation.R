#!/usr/bin/env Rscript
#
# Interim timing calibration (Section 5.2; Table 2 / Figure 1) -- HPC driver.
#
# This script contains NO function definitions of its own. It only:
#   (1) defines the scenario-specific inputs (control/effect/recruitment
#       models, the true design's efficacy/final boundaries),
#   (2) calls DTEAssurance::calibrate_BPP_timing(), and
#   (3) saves the result (which already records its own settings/provenance).
#
# Usage: Rscript 01_run_pp_timing_simulation.R <n_sims> <seed>
#   n_sims = number of interim datasets simulated PER candidate IF
#   seed   = integer seed, used both for reproducibility and as the output
#            file suffix
#
# Submitted at paper scale by 01_run_pp_timing_simulation.sbatch
# (50 tasks x 100 = 5,000 datasets per candidate IF). Writes
# PP_timing_batch_<seed>.rds to the working directory.

library(DTEAssurance)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript 01_run_pp_timing_simulation.R <n_sims> <seed>")
}
n_sims <- as.numeric(args[1])
seed   <- as.numeric(args[2])

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))

# --- Scenario-specific inputs (match the manuscript case study exactly) ---

n_c <- 600
n_t <- 600

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
                       alpha = 0.025,
                       events = 840)

# --- Candidate interim timings for the sweep ---
# Trimmed to stop at IF = 0.7: candidates at or beyond the true design's own
# efficacy look (IF = 0.75) are not a well-posed question about this design
# (see discussion in the revision plan) and calibrate_BPP_timing() will
# error if a candidate at/after a supplied future boundary is included.

IA_model <- list(events = 840,
                 IF = seq(0.2, 0.7, by = 0.1))

# --- Build the true design's future boundaries (efficacy look + final) ---
# Uses the package's own design-construction helper so these boundaries are
# guaranteed consistent with the real design evaluated elsewhere (Table 4/5),
# rather than hardcoded numbers that could silently drift out of sync.

GSD_model_efficacy <- list(events = 840,
                           alpha_spending = c(0.0125, 0.025),
                           alpha_IF = c(0.75, 1),
                           futility_type = "none")

rpact_out <- DTEAssurance:::make_rpact_design_from_GSD_model(GSD_model_efficacy)
design    <- rpact_out$design

future_boundaries <- list(
  list(events = ceiling(0.75 * IA_model$events),
       crit   = design$criticalValues[which(abs(design$informationRates - 0.75) < 1e-8)]),
  list(events = IA_model$events,
       crit   = design$criticalValues[which(abs(design$informationRates - 1) < 1e-8)])
)

cat(sprintf("Future boundaries: IF=0.75 -> Z>%.4f ; final -> Z>%.4f\n",
            future_boundaries[[1]]$crit, future_boundaries[[2]]$crit))
cat("Sanity check: these should match the manuscript's 2.241 / 2.047.\n\n")

# --- Run the timing sweep via the package's own (patched) function ---

result <- DTEAssurance::calibrate_BPP_timing(
  n_c = n_c,
  n_t = n_t,
  control_model = control_model,
  effect_model = effect_model,
  recruitment_model = recruitment_model,
  IA_model = IA_model,
  analysis_model = analysis_model,
  future_boundaries = future_boundaries,
  update_priors_sims = 1000,
  PP_sims = 2000,
  n_sims = n_sims,
  n_cores = n_cores,
  seed = seed
)

outname <- sprintf("PP_timing_batch_%04d.rds", seed)
saveRDS(result, outname)
