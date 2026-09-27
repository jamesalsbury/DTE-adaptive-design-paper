#!/usr/bin/env Rscript
#
# Local analysis of PP_Timing results.
#
# Run this from the directory containing either:
#   (a) the 50 raw PP_timing_batch_*.rds files, OR
#   (b) an already-combined PP_timing_combined.rds
# It will auto-detect which you have and aggregate if needed.
#
# Produces:
#   - a console summary (mean, SD, quartiles) per candidate IF
#   - a PNG with one histogram panel per candidate IF (Figure 1 draft)
#
# Usage: Rscript analyze_PP_timing_local.R [directory]  (default: current dir)

args <- commandArgs(trailingOnly = TRUE)
work_dir <- if (length(args) >= 1) args[1] else "."

combined_file <- file.path(work_dir, "PP_timing_combined.rds")
batch_files <- list.files(work_dir, pattern = "^PP_timing_batch_\\d+\\.rds$", full.names = TRUE)

if (file.exists(combined_file)) {
  
  cat("Found existing combined file:", combined_file, "\n")
  result <- readRDS(combined_file)
  
} else if (length(batch_files) > 0) {
  
  cat(sprintf("Found %d raw batch files -- aggregating now.\n", length(batch_files)))
  
  batches <- lapply(batch_files, readRDS)
  
  # --- Consistency checks (same logic as aggregate_PP_timing.R) ---
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
  
  ok <- check_field_consistent("future_boundaries", function(b) b$settings$future_boundaries) &
    check_field_consistent("IA_model", function(b) b$settings$IA_model) &
    check_field_consistent("update_priors_sims", function(b) b$settings$update_priors_sims) &
    check_field_consistent("PP_sims", function(b) b$settings$PP_sims)
  
  if (!ok) {
    stop("Aborting: inconsistent settings across batches (see warnings above). ",
         "Resolve before analysing.")
  }
  cat("All batches consistent on future_boundaries, IA_model, update_priors_sims, PP_sims.\n")
  
  pkg_versions <- unique(vapply(batches, function(b) b$settings$package_version, character(1)))
  if (length(pkg_versions) > 1) {
    warning("Batches produced under multiple package versions: ", paste(pkg_versions, collapse = ", "))
  } else {
    cat(sprintf("All batches produced under DTEAssurance version %s.\n", pkg_versions))
  }
  
  n_candidates <- length(ref_settings$IA_model$IF)
  combined_outcome_list <- vector("list", n_candidates)
  for (i in seq_len(n_candidates)) {
    combined_outcome_list[[i]] <- list(
      BPP_values = unlist(lapply(batches, function(b) b$outcome_list[[i]]$BPP_values)),
      cens_time  = unlist(lapply(batches, function(b) b$outcome_list[[i]]$cens_time)),
      IF         = ref_settings$IA_model$IF[i]
    )
  }
  
  result <- list(outcome_list = combined_outcome_list, settings = ref_settings,
                 n_batches_combined = length(batch_files))
  
  saveRDS(result, combined_file)
  cat("Saved combined result to", combined_file, "\n")
  
} else {
  stop("No PP_timing_batch_*.rds or PP_timing_combined.rds files found in '", work_dir, "'.")
}

# --- Summary table ---

cat("\n=== Per-candidate-IF summary ===\n")
summary_df <- do.call(rbind, lapply(result$outcome_list, function(x) {
  q <- quantile(x$BPP_values, probs = c(0.10, 0.25, 0.50, 0.75, 0.90), na.rm = TRUE)
  data.frame(
    IF = x$IF,
    n = length(x$BPP_values),
    mean = mean(x$BPP_values, na.rm = TRUE),
    sd = sd(x$BPP_values, na.rm = TRUE),
    q10 = q[1], q25 = q[2], median = q[3], q75 = q[4], q90 = q[5],
    prop_extreme = mean(x$BPP_values < 0.1 | x$BPP_values > 0.9, na.rm = TRUE),
    mean_cens_time = mean(x$cens_time, na.rm = TRUE)
  )
}))
print(summary_df, row.names = FALSE, digits = 3)

# --- Histograms: one panel per candidate IF ---

n_candidates <- length(result$outcome_list)
png_file <- file.path(work_dir, "PP_timing_histograms.png")
png(png_file, width = 1000, height = 700, res = 120)
par(mfrow = c(2, ceiling(n_candidates / 2)), mar = c(4, 4, 3, 1))

for (x in result$outcome_list) {
  hist(x$BPP_values, breaks = 20, xlim = c(0, 1),
       main = sprintf("IF = %.1f (n=%d)", x$IF, length(x$BPP_values)),
       xlab = "PP", col = "grey80", border = "white")
}

dev.off()
cat("\nSaved histogram panel to", png_file, "\n")
cat("\nLook for: increasingly bimodal/polarised shape as IF increases,\n")
cat("consistent with 'approximately uniform early, polarising toward 0/1 later'.\n")