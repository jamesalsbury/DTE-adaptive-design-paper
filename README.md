# DTE-adaptive-design-paper

Companion analysis repository for:

> Salsbury, J.A., Oakley, J.E., Julious, S.A., Hampson, L.V. "Adaptive
> clinical trial design with delayed treatment effects using elicited
> prior distributions."

This repository contains the reconstructed control-arm data, and every
script used to produce the tables and figures in the case study
(Section 5). All analyses use the `DTEAssurance` R package
(https://github.com/jamesalsbury/DTEAssurance).

## Requirements

```r
install.packages("devtools")
devtools::install_github("jamesalsbury/DTEAssurance@<tag>")  # pin to the version used for this paper
```

Also requires: `rjags` (and JAGS itself), `SHELF`, `nphRCT`, `nph`,
`survival`, `nleqslv`, `parallel`.

## Repository structure

```
data/                     Reconstructed IPD (ZODIAC, REVEL, INTEREST)
control_prior/             Pooling, MCMC fit, Beta-Beta approximation + validation
pp_timing/                 Table 2 / Figure 1 (interim timing calibration)
threshold_calibration/     Table 3, kappa* selection (Algorithm 3)
matched_comparators/       D4/D5 futility boundary calibration
table4_5/                  Table 4/5 (main operating-characteristics study)
posterior_states/          Appendix posterior state probability illustration
slurm/                     HPC submission scripts for all of the above
```

## Script-to-output mapping

| Output | Script(s) | Notes |
|---|---|---|
| Control-arm prior hyperparameters (Section 5.1) | `control_prior/pool_and_fit.R`, `control_prior/fit_control_prior_betas.R` | IPD reconstructed by hand from published KM curves + at-risk tables (Guyot et al. 2012); validity check included |
| Figure 1 / Table 2 (interim timing) | `pp_timing/PP_Timing.R`, `pp_timing/aggregate_PP_timing.R` | 50-task SLURM array, 5,000 replicates/candidate IF |
| Table 3 / kappa* selection (Section 5.3.2) | `threshold_calibration/Eq11_Calibration.R`, `threshold_calibration/aggregate_eq11.R` | 50-task SLURM array, 15,000 replicates/scenario |
| D4/D5 futility boundaries | `matched_comparators/D4_D5_calibration.R` | No HPC needed; n=20,000/scenario, local |
| Table 4 (S1-S3) | `table4_5/Table4_S1.R`, `Table4_S2.R`, `Table4_S3.R` + `aggregate_table4_5.R` | 50-task SLURM arrays each, 2,000 replicates/task |
| Table 4 (S4) + Table 5 | `table4_5/Table4_S4.R` (stratified) + `table4_5/reweight_S4.R` | S4 uses fixed-proportion stratified sampling; reweighted post-hoc for both the true elicited (P_S, P_DTE) and the sensitivity-analysis values |
| Appendix posterior-state tables | `posterior_states/representative_posterior_states.R`, `representative_posterior_states_multiseed.R` | No HPC needed |
| Convergence-by-decision check | `table4_5/check_convergence_by_decision.R` | Run once per scenario on each `table4_S{n}_combined.rds` |

## Reproducing a table from scratch

Each SLURM-based script in `slurm/` follows the same pattern:
1. A smoke test (`*_smoke.sbatch`, tiny `n_sims`) - confirm the pipeline runs correctly before committing real compute
2. The full array (`*_array.sbatch`) - real replicate counts as used in the paper
3. The corresponding `aggregate_*.R` script - combines batch outputs, checks settings consistency across all batches before combining

All real-scale runs record full provenance (settings, seed, package
version, timestamp) alongside the results, so any saved `.rds` output
can be traced back to exactly what produced it.
