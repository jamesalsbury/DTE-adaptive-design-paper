#!/usr/bin/env Rscript
#
# Aggregates all table4_S{scenario}_batch_*.rds files for one scenario into
# a combined dataset and a Table-4-formatted summary (one row per design).
# Works for S1/S2/S3 (fixed scenarios) and S4 (stratified) alike -- S4's
# extra `state` column is carried through if present, but not required.
#
# Run once per scenario (S1-S4). For S4 the Table 4 row comes from
# 06_reweight_for_table4_and_table5.R, which reads this script's output.
#
# Usage: Rscript 05_aggregate_table4.R <scenario e.g. S1> <directory> <output_file>
# e.g.:  Rscript 05_aggregate_table4.R S1 . table4_S1_combined.rds

args <- commandArgs(trailingOnly = TRUE)
scenario  <- if (length(args) >= 1) args[1] else stop("Specify scenario, e.g. S1")
batch_dir <- if (length(args) >= 2) args[2] else "."
outname   <- if (length(args) >= 3) args[3] else sprintf("table4_%s_combined.rds", scenario)

pattern <- sprintf("^table4_%s_batch_\\d+\\.rds$", scenario)
batch_files <- list.files(batch_dir, pattern = pattern, full.names = TRUE)

if (length(batch_files) == 0) {
  stop("No files matching '", pattern, "' found in '", batch_dir, "'.")
}

cat(sprintf("Found %d batch files for scenario %s.\n", length(batch_files), scenario))

batches <- lapply(batch_files, readRDS)

# --- Consistency checks across batches ---

ref_settings <- batches[[1]]$settings

check_field_consistent <- function(field_name, extractor) {
  ref <- extractor(batches[[1]])
  mismatched <- which(!vapply(batches, function(b) identical(extractor(b), ref), logical(1)))
  if (length(mismatched) > 0) {
    warning(sprintf("Inconsistent '%s' across batches: files %s differ from the first batch.",
                    field_name, paste(basename(batch_files[mismatched]), collapse = ", ")))
    return(FALSE)
  }
  TRUE
}

ok <- check_field_consistent("total_events",   function(b) b$settings$total_events) &
  check_field_consistent("futility_IF",    function(b) b$settings$futility_IF) &
  check_field_consistent("kappa_star",     function(b) b$settings$kappa_star) &
  check_field_consistent("D4_boundary_Z",  function(b) b$settings$D4_boundary_Z) &
  check_field_consistent("D5_boundary_Z",  function(b) b$settings$D5_boundary_Z) &
  check_field_consistent("t_star_D5",      function(b) b$settings$t_star_D5)

pkg_versions <- unique(vapply(batches, function(b) b$settings$package_version, character(1)))
if (length(pkg_versions) > 1) {
  warning("Batches produced under multiple DTEAssurance versions: ",
          paste(pkg_versions, collapse = ", "))
}

seeds <- vapply(batches, function(b) b$settings$seed, numeric(1))
if (length(unique(seeds)) != length(seeds)) {
  warning("Duplicate seeds detected across batch files -- check for accidental double-counting.")
}

if (!ok) {
  stop("Aborting: inconsistent design settings across batches (see warnings above). ",
       "Resolve before aggregating.")
}

cat("All batches consistent on total_events, futility_IF, kappa_star, D4/D5 boundaries, t_star.\n")
if (length(pkg_versions) == 1) cat(sprintf("Package version: %s\n", pkg_versions))
cat(sprintf("Seeds combined: %s\n\n", paste(sort(seeds), collapse = ",")))

# --- Combine, with globally unique replicate IDs (seed * 1e6 + local replicate) ---

combined <- do.call(rbind, lapply(batches, function(b) {
  df <- b$results
  df$global_replicate <- b$settings$seed * 1e6 + df$replicate
  df
}))

n_total_replicates <- length(unique(combined$global_replicate))
cat(sprintf("Total combined replicates: %d\n\n", n_total_replicates))

# --- Table 4-style summary: one row per design ---

summary_tab <- do.call(rbind, lapply(sort(unique(combined$design)), function(d) {
  sub <- combined[combined$design == d, ]
  data.frame(
    Design = d,
    `P(Reject H0)` = mean(sub$success),
    `P(Early Fut.)` = if (d == "D1") NA else mean(sub$early_fut),
    `P(Early Eff.)` = if (d == "D1") NA else mean(sub$early_eff),
    ESS = mean(sub$sample_size),
    Duration = mean(sub$duration),
    n = nrow(sub),
    check.names = FALSE
  )
}))

cat(sprintf("=== Table 4, Scenario %s (n=%d replicates per design) ===\n", scenario, n_total_replicates))
print(summary_tab, row.names = FALSE, digits = 4)

# --- D3 convergence rate, if recorded ---
d3_conv <- combined$converged[combined$design == "D3"]
if (any(!is.na(d3_conv))) {
  cat(sprintf("\nD3 MCMC convergence rate: %.2f%% (%d/%d fits)\n",
              100 * mean(d3_conv, na.rm = TRUE), sum(d3_conv, na.rm = TRUE),
              sum(!is.na(d3_conv))))
}

saveRDS(list(combined = combined, summary = summary_tab, settings = ref_settings,
             n_batches = length(batch_files), source_files = basename(batch_files)),
        outname)
cat(sprintf("\nSaved to %s\n", outname))
