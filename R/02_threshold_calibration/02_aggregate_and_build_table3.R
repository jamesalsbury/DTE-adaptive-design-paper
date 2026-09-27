#!/usr/bin/env Rscript
#
# Aggregation script for the Eq. (11) calibration SLURM array output.
# Combines all eq11_calibration_batch_*.rds files (one per array task) into
# a single set of raw per-replicate data frames, one per scenario (null,
# S2, S3), ready to pass to summarize_grid_by_lambda() / select_lambda_star().
#
# As with 01_interim_timing/02_aggregate_and_plot_figure1.R, this checks that every batch used
# consistent settings before combining, and aborts rather than silently
# mixing batches run under different settings.
#
# Usage: Rscript 02_aggregate_and_build_table3.R <directory containing batch files> <output_file>
# e.g.:  Rscript 02_aggregate_and_build_table3.R . eq11_combined.rds

args <- commandArgs(trailingOnly = TRUE)
batch_dir <- if (length(args) >= 1) args[1] else "."
outname   <- if (length(args) >= 2) args[2] else "eq11_combined.rds"

batch_files <- list.files(batch_dir, pattern = "^eq11_calibration_batch_\\d+\\.rds$", full.names = TRUE)

if (length(batch_files) == 0) {
  stop("No files matching 'eq11_calibration_batch_*.rds' found in '", batch_dir, "'.")
}

cat(sprintf("Found %d batch files.\n", length(batch_files)))

batches <- lapply(batch_files, readRDS)

scenario_names <- names(batches[[1]]$grid_runs)

# --- Consistency checks across batches ---

check_field_consistent <- function(field_name, extractor) {
  ref <- extractor(batches[[1]])
  mismatched <- which(!vapply(batches, function(b) identical(extractor(b), ref), logical(1)))
  if (length(mismatched) > 0) {
    warning(sprintf(
      "Inconsistent '%s' across batches: files %s differ from the first batch (%s).",
      field_name,
      paste(basename(batch_files[mismatched]), collapse = ", "),
      basename(batch_files[1])
    ))
    return(FALSE)
  }
  TRUE
}

ok_lambda_grid <- check_field_consistent("lambda_grid", function(b) b$lambda_grid)
ok_power_floor <- check_field_consistent("power_floor", function(b) b$power_floor)
ok_alt_scenarios <- check_field_consistent("alt_scenarios", function(b) b$alt_scenarios)

ok_per_scenario <- TRUE
for (scen in scenario_names) {
  ok_per_scenario <- ok_per_scenario &
    check_field_consistent(paste0(scen, " future_boundaries"),
                           function(b) b$grid_runs[[scen]]$settings$future_boundaries) &
    check_field_consistent(paste0(scen, " futility_IF"),
                           function(b) b$grid_runs[[scen]]$settings$futility_IF) &
    check_field_consistent(paste0(scen, " data_generating_model"),
                           function(b) b$grid_runs[[scen]]$settings$data_generating_model) &
    check_field_consistent(paste0(scen, " update_priors_sims"),
                           function(b) b$grid_runs[[scen]]$settings$update_priors_sims) &
    check_field_consistent(paste0(scen, " PP_sims"),
                           function(b) b$grid_runs[[scen]]$settings$PP_sims)
}

pkg_versions <- unique(vapply(batches, function(b) {
  v <- b$grid_runs[[scenario_names[1]]]$settings$package_version
  if (is.null(v)) NA_character_ else v
}, character(1)))
if (length(pkg_versions) > 1) {
  warning("Batches were produced under multiple DTEAssurance package versions: ",
          paste(pkg_versions, collapse = ", "),
          ". Results may not be directly comparable across batches.")
}

all_consistent <- ok_lambda_grid && ok_power_floor && ok_alt_scenarios && ok_per_scenario
if (!all_consistent) {
  stop("Aborting aggregation: one or more settings were inconsistent across ",
       "batches (see warnings above). Resolve before combining.")
}

cat("All batches consistent on lambda_grid, power_floor, alt_scenarios, and per-scenario settings.\n")
if (length(pkg_versions) == 1) {
  cat(sprintf("All batches produced under DTEAssurance version %s.\n", pkg_versions))
}

# --- Combine: for each scenario, rbind raw data across all batches ---

combined_raw <- list()
for (scen in scenario_names) {
  combined_raw[[scen]] <- do.call(rbind, lapply(batches, function(b) b$grid_runs[[scen]]$raw))
  cat(sprintf("  %s: %d combined replicates (mean BPP_val = %.3f)\n",
              scen, nrow(combined_raw[[scen]]), mean(combined_raw[[scen]]$BPP_val)))
}

combined <- list(
  raw_by_scenario = combined_raw,
  lambda_grid = batches[[1]]$lambda_grid,
  power_floor = batches[[1]]$power_floor,
  alt_scenarios = batches[[1]]$alt_scenarios,
  n_batches_combined = length(batch_files),
  source_files = basename(batch_files)
)

saveRDS(combined, outname)
cat(sprintf("\nSaved combined result to '%s'.\n", outname))
cat("\nNext steps:\n")
cat("  summaries <- lapply(combined$raw_by_scenario, summarize_grid_by_lambda, lambda_grid = combined$lambda_grid)\n")
cat("  names(summaries) <- names(combined$raw_by_scenario)\n")
cat("  selection <- select_lambda_star(summaries, power_floor = combined$power_floor,\n")
cat("                                 null_scenario = 'null', alt_scenarios = combined$alt_scenarios)\n")
