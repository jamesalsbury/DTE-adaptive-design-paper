# DTE-adaptive-design-paper

Companion analysis repository for:

> Salsbury, J.A., Oakley, J.E., Julious, S.A., Hampson, L.V. "Adaptive
> clinical trial design with delayed treatment effects using elicited
> prior distributions."

This repository contains the reconstructed control-arm data and the
scripts that reproduce the tables and figures in the case study
(Section 5). All analyses use the `DTEAssurance` R package
(https://github.com/jamesalsbury/DTEAssurance); no package functions are
redefined here.

## Requirements

```r
install.packages("devtools")
devtools::install_github("jamesalsbury/DTEAssurance@<tag>")  # pin to the version used for this paper
```

Also requires: `rjags` (and JAGS itself), `SHELF`, `nphRCT`, `nph`,
`survival`, `nleqslv`, `parallel`.

## Repository structure

The numbered folders follow the order in which Section 5 builds its
results; each stage uses outputs of the stages before it.

```
data/
  zodiac_control_arm.csv        Reconstructed docetaxel control-arm IPD (Herbst et al.)
  revel_control_arm.csv         Reconstructed docetaxel control-arm IPD (Garon et al.)
  interest_control_arm.csv      Reconstructed docetaxel control-arm IPD (Kim et al.)

R/
  00_control_prior/             Pooled Weibull MCMC fit -> Beta-Beta control prior
  01_interim_timing/            Interim timing calibration (IF = 0.50)
  02_threshold_calibration/     PP threshold calibration (kappa* = 0.20)
  03_matched_comparators/       D4/D5 matched futility boundaries
  04_operating_characteristics/ Full D1-D5 evaluation under S1-S4
  05_posterior_state_illustration/  Posterior state probabilities (supplementary)
```

## Script-to-output mapping

| Output | Script(s) | Compute |
|---|---|---|
| Control-arm prior hyperparameters (Section 5.1) | `00_control_prior/01_pool_and_fit_weibull_posterior.R`, `00_control_prior/02_fit_beta_approximation.R` | Local |
| Table 2 / Figure 1 (interim timing) | `01_interim_timing/01_run_pp_timing_simulation.R` (+ `.sbatch`), `01_interim_timing/02_aggregate_and_plot_figure1.R` | HPC: 50 tasks x 100 = 5,000 datasets per candidate IF |
| Table 3 / kappa* selection (incl. D2 power-floor check) | `02_threshold_calibration/01_run_threshold_calibration.R` (+ `.sbatch`), `02_aggregate_and_build_table3.R` | HPC: 50 tasks x 300 = 15,000 per scenario |
| D4/D5 futility boundaries | `03_matched_comparators/01_calibrate_D4_D5_boundaries.R` | Local, 20,000 datasets per scenario |
| Table 4 (S1-S3) | `04_operating_characteristics/0{1,2,3}_run_scenario_*.R` (+ `.sbatch`), `05_aggregate_table4.R` | HPC: 50 tasks x 2,000 = 100,000 per scenario |
| Table 4 (S4) + Table 5 | `04_operating_characteristics/04_run_scenario_prior_predictive.R` (+ `.sbatch`), `05_aggregate_table4.R`, `06_reweight_for_table4_and_table5.R` | HPC: 50 tasks x 2,000 = 100,000 |
| D3 convergence-by-decision check | `04_operating_characteristics/07_check_convergence_by_decision.R` | Local |
| Posterior state tables (main text + Appendix robustness check) | `05_posterior_state_illustration/01_representative_and_multiseed_check.R` | Local |

## Reproducing the results

Run every script from its own folder; outputs are written to the working
directory.

1. **Control prior.** In `R/00_control_prior/`, run `01_...R` (writes
   `weibull_posterior_samples.rds`), then `02_...R` (Beta-Beta
   hyperparameters and validity checks).
2. **HPC stages.** Each HPC stage has one `.sbatch` file next to the R
   script it runs, set to the paper-reported scale (50-task array). Submit
   it from that folder, e.g.
   `cd R/01_interim_timing && sbatch 01_run_pp_timing_simulation.sbatch`.
   For a quick test, reduce `NSIMS` to ~5 and the array to 1-2 tasks.
3. **Aggregate.** Combine the batch files with the stage's aggregation
   script, which checks that every batch was run with identical settings
   before combining:
   - `Rscript 02_aggregate_and_plot_figure1.R` (in `01_interim_timing/`)
   - `Rscript 02_aggregate_and_build_table3.R` (in `02_threshold_calibration/`)
   - `Rscript 05_aggregate_table4.R S1` (and `S2`, `S3`, `S4`) in
     `04_operating_characteristics/`, then
     `Rscript 06_reweight_for_table4_and_table5.R table4_S4_combined.rds`
4. **Local stages.** `03_matched_comparators/` and
   `05_posterior_state_illustration/` run directly with `Rscript`.

All HPC outputs record their settings, seed, package version and
timestamp, so any saved `.rds` file can be traced back to what produced it.
