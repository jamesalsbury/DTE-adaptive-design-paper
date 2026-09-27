#!/usr/bin/env Rscript
#
# Reweights S4's fixed-proportion-stratified output by the TRUE elicited
# (or alternative, for Table 5) P_S/P_DTE, per Equation (9):
#   A_overall = (1-P_S)*A_null + P_S*P_DTE*A_delay + P_S*(1-P_DTE)*A_PH
#
# Table 4's S4 row uses the manuscript's actual elicited values
# (P_S=0.9, P_DTE=0.7). Table 5 repeats this at alternative (P_S,P_DTE)
# pairs, reusing the SAME simulated data (no re-simulation needed) --
# this is only valid because S4 was deliberately stratified into equal
# thirds by state, so each state's summary is estimated independently
# and precisely, then recombined with different weights.
#
# Usage: Rscript 06_reweight_for_table4_and_table5.R <path to table4_S4_combined.rds>
# (table4_S4_combined.rds is produced by 05_aggregate_table4.R S4)

args <- commandArgs(trailingOnly = TRUE)
combined_file <- if (length(args) >= 1) args[1] else "table4_S4_combined.rds"

result <- readRDS(combined_file)
combined <- result$combined

if (!"state" %in% names(combined)) {
  stop("No 'state' column found -- is this really the S4 combined output?")
}

cat(sprintf("Loaded %d rows, %d unique replicates.\n",
            nrow(combined), length(unique(combined$global_replicate))))
cat(sprintf("State distribution (should be ~equal thirds): %s\n\n",
            paste(table(combined$state[combined$design == "D1"]), collapse = " / ")))

# --- Per-state, per-design summary (the building blocks for reweighting) ---

per_state_summary <- function(data, design_label) {
  sub <- data[data$design == design_label, ]
  do.call(rbind, lapply(1:3, function(s) {
    ss <- sub[sub$state == s, ]
    data.frame(
      state = s,
      n = nrow(ss),
      P_reject = mean(ss$success),
      P_early_fut = if (design_label == "D1") NA else mean(ss$early_fut),
      P_early_eff = if (design_label == "D1") NA else mean(ss$early_eff),
      ESS = mean(ss$sample_size),
      Duration = mean(ss$duration)
    )
  }))
}

# --- Reweight per Eq. (9), given P_S, P_DTE ---
# state 1 = null (weight 1-P_S), state 2 = immediate (weight P_S*(1-P_DTE)),
# state 3 = delayed (weight P_S*P_DTE)

reweight <- function(per_state, P_S, P_DTE) {
  w <- c(1 - P_S, P_S * (1 - P_DTE), P_S * P_DTE)
  stopifnot(abs(sum(w) - 1) < 1e-8)
  
  cols <- c("P_reject", "P_early_fut", "P_early_eff", "ESS", "Duration")
  out <- sapply(cols, function(col) {
    vals <- per_state[[col]]
    if (all(is.na(vals))) return(NA_real_)
    sum(w * vals)
  })
  as.data.frame(t(out))
}

designs <- c("D1", "D2", "D3", "D4", "D5")

per_state_all <- lapply(designs, function(d) per_state_summary(combined, d))
names(per_state_all) <- designs

cat("=== Per-state summaries (building blocks -- sanity check) ===\n\n")
for (d in designs) {
  cat(sprintf("-- %s --\n", d))
  print(per_state_all[[d]], row.names = FALSE, digits = 4)
  cat("\n")
}

# --- Table 4's S4 row: reweight at the TRUE elicited P_S=0.9, P_DTE=0.7 ---

cat("=== Table 4, Scenario S4 (reweighted at true P_S=0.9, P_DTE=0.7) ===\n\n")

table4_S4 <- do.call(rbind, lapply(designs, function(d) {
  rw <- reweight(per_state_all[[d]], P_S = 0.9, P_DTE = 0.7)
  data.frame(Design = d, rw, check.names = FALSE)
}))
print(table4_S4, row.names = FALSE, digits = 4)

# --- Table 5: reweight at alternative (P_S, P_DTE) pairs ---

cat("\n\n=== Table 5: alternative (P_S, P_DTE) specifications ===\n\n")

alt_scenarios <- list(
  c(P_S = 0.8, P_DTE = 0.7),
  c(P_S = 0.9, P_DTE = 0.8),
  c(P_S = 0.8, P_DTE = 0.8)
)

table5 <- do.call(rbind, lapply(alt_scenarios, function(sc) {
  do.call(rbind, lapply(designs, function(d) {
    rw <- reweight(per_state_all[[d]], P_S = sc["P_S"], P_DTE = sc["P_DTE"])
    data.frame(P_S = sc["P_S"], P_DTE = sc["P_DTE"], Design = d, rw,
               check.names = FALSE, row.names = NULL)
  }))
}))
print(table5, row.names = FALSE, digits = 4)

saveRDS(list(per_state_all = per_state_all, table4_S4 = table4_S4, table5 = table5),
        "table4_5_S4_reweighted.rds")
cat("\nSaved to table4_5_S4_reweighted.rds\n")
