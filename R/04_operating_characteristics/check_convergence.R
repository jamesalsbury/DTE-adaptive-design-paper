#!/usr/bin/env Rscript
#
# Check whether D3's non-converged interim MCMC fits are concentrated
# among particular decision outcomes, or spread roughly evenly.
#
# Run this in the same session as aggregate_table4_5.R (uses the
# `combined` object it produces), or load a saved *_combined.rds file.

# If not already in your session:
# result <- readRDS("table4_S2_combined.rds")
# combined <- result$combined

d3 <- combined[combined$design == "D3", ]

cat("=== D3 decisions, split by convergence status ===\n\n")

tab <- table(d3$decision, d3$converged, useNA = "ifany")
print(tab)

cat("\n=== Proportion converged within each decision type ===\n\n")
prop_tab <- prop.table(tab, margin = 1)
print(round(prop_tab, 4))

cat("\n=== Overall non-convergence rate by decision type ===\n\n")
by_decision <- aggregate(converged ~ decision, data = d3, FUN = function(x) {
  c(n = length(x), pct_converged = round(100 * mean(x), 2))
})
print(by_decision)

# Chi-square test: is non-convergence independent of decision type?
cat("\n=== Chi-square test: convergence independent of decision? ===\n\n")
chisq_result <- chisq.test(table(d3$decision, d3$converged))
print(chisq_result)

cat("\nInterpretation:\n")
cat("- If 'pct_converged' is similar (~99%) across all decision types,\n")
cat("  non-convergence is roughly evenly spread -- not concentrated in\n")
cat("  borderline/futility decisions specifically.\n")
cat("- If one decision type (e.g. 'Stop for futility') has a notably\n")
cat("  lower convergence rate than others, that's worth a closer look --\n")
cat("  though even then, the chi-square p-value tells you whether this\n")
cat("  difference is likely real or could be chance given 100,000 reps.\n")