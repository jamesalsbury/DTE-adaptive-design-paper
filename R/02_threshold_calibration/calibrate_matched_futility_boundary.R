#' (internal) Simulate one interim dataset and compute Z at the futility look
#'
#' @keywords internal
single_matched_futility_rep <- function(i, n_c, n_t, data_generating_model,
                                        recruitment_model, futility_IF, total_events,
                                        analysis_model, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed * 10000 + i)
  
  if (is.null(data_generating_model$gamma_c)) {
    trial_data <- sim_dte(n_c, n_t, data_generating_model$lambda_c,
                          delay_time = data_generating_model$delay_time,
                          post_delay_HR = data_generating_model$post_delay_HR,
                          dist = "Exponential")
  } else {
    trial_data <- sim_dte(n_c, n_t, data_generating_model$lambda_c,
                          delay_time = data_generating_model$delay_time,
                          post_delay_HR = data_generating_model$post_delay_HR,
                          dist = "Weibull", gamma_c = data_generating_model$gamma_c)
  }
  
  trial_data <- add_recruitment_time(trial_data,
                                     rec_method   = recruitment_model$method,
                                     rec_period   = recruitment_model$period,
                                     rec_power    = recruitment_model$power,
                                     rec_rate     = recruitment_model$rate,
                                     rec_duration = recruitment_model$duration)
  
  n_events_interim <- ceiling(futility_IF * total_events)
  censored <- cens_data(trial_data[order(trial_data$pseudo_time), ],
                        cens_method = "Events", cens_events = n_events_interim)
  
  Z <- survival_test(censored$data,
                     analysis_method = analysis_model$method,
                     alpha           = analysis_model$alpha,
                     alternative     = analysis_model$alternative_hypothesis,
                     rho             = analysis_model$rho,
                     gamma           = analysis_model$gamma,
                     t_star          = analysis_model$t_star,
                     s_star          = analysis_model$s_star)$Z
  
  data.frame(Z = Z)
}


#' Calibrate a matched-null-rate futility boundary for a frequentist design
#'
#' Finds the Z-statistic threshold at the futility look such that stopping
#' whenever \code{Z < threshold} reproduces a target null futility-stopping
#' rate (typically taken to match another design's, e.g. the PP-based
#' design D3's, null futility rate, for a like-for-like comparison), then
#' reports the futility-triggering rate that same threshold produces under
#' each of the other supplied scenarios. This is the general mechanism
#' underlying both the "matched log-rank futility design" (D4, using
#' \code{analysis_model$method = "LRT"}) and its modestly-weighted
#' log-rank analogue (D5, using \code{analysis_model$method = "MW"}) --
#' the two designs differ only in \code{analysis_model}, not in this
#' calibration procedure.
#'
#' @param n_c,n_t Control/treatment sample sizes.
#' @param recruitment_model See \code{\link{add_recruitment_time}}.
#' @param futility_IF Information fraction of the futility look.
#' @param total_events Total planned event count (100% information fraction).
#' @param analysis_model The analysis method whose Z statistic defines the
#'   boundary -- see \code{\link{survival_test}}. Set
#'   \code{method = "LRT"} for the D4 analogue, \code{method = "MW"}
#'   (with \code{t_star} set) for D5.
#' @param target_null_futility_rate The null futility-stopping rate to
#'   match (e.g. D3's null futility rate from the real simulation study --
#'   this must be supplied, not assumed, and should be updated if D3's
#'   rate changes on re-run).
#' @param scenarios A named list of data-generating scenarios (each a list
#'   with \code{lambda_c}, \code{gamma_c}, \code{delay_time},
#'   \code{post_delay_HR}) to report the calibrated boundary's behaviour
#'   under. Must include an entry named \code{"null"}, used to calibrate
#'   the boundary itself.
#' @param n_sims Number of interim datasets to simulate per scenario
#'   (default 2000, matching the precision used for D4's original
#'   calibration).
#' @param n_cores Number of cores for \code{parallel::mclapply} (default 1).
#' @param seed Optional integer seed.
#'
#' @return A list with:
#'   \describe{
#'     \item{boundary}{The calibrated Z threshold.}
#'     \item{scenario_futility_rates}{A named numeric vector, one entry
#'       per scenario, giving \code{mean(Z < boundary)} under that
#'       scenario -- for \code{"null"} this equals (up to Monte Carlo
#'       error) \code{target_null_futility_rate} by construction; for
#'       other scenarios it is the false-futility rate the calibrated
#'       design would exhibit.}
#'     \item{raw_Z_by_scenario}{The raw simulated Z values per scenario,
#'       for further inspection/plotting.}
#'     \item{settings}{Provenance: all arguments, plus the installed
#'       \code{DTEAssurance} package version.}
#'   }
#'
#' @export
calibrate_matched_futility_boundary <- function(n_c, n_t,
                                                recruitment_model,
                                                futility_IF, total_events,
                                                analysis_model,
                                                target_null_futility_rate,
                                                scenarios,
                                                n_sims = 2000,
                                                n_cores = 1,
                                                seed = NULL) {
  
  if (!"null" %in% names(scenarios)) {
    stop("calibrate_matched_futility_boundary: 'scenarios' must include an ",
         "entry named 'null', used to calibrate the boundary.")
  }
  
  raw_Z_by_scenario <- list()
  
  for (scen_name in names(scenarios)) {
    run_one <- function(i) {
      single_matched_futility_rep(
        i, n_c = n_c, n_t = n_t,
        data_generating_model = scenarios[[scen_name]],
        recruitment_model = recruitment_model,
        futility_IF = futility_IF, total_events = total_events,
        analysis_model = analysis_model, seed = seed
      )
    }
    
    if (n_cores > 1) {
      result <- parallel::mclapply(seq_len(n_sims), run_one, mc.cores = n_cores)
    } else {
      result <- lapply(seq_len(n_sims), run_one)
    }
    
    raw_Z_by_scenario[[scen_name]] <- do.call(rbind, result)$Z
  }
  
  # Calibrate: boundary is the target_null_futility_rate-quantile of the
  # null Z distribution, since P(Z < boundary | null) = target rate by
  # definition of the quantile function.
  boundary <- stats::quantile(raw_Z_by_scenario[["null"]],
                              probs = target_null_futility_rate, na.rm = TRUE)
  boundary <- as.numeric(boundary)
  
  scenario_futility_rates <- vapply(raw_Z_by_scenario, function(Z) {
    mean(Z < boundary, na.rm = TRUE)
  }, numeric(1))
  
  settings <- list(
    futility_IF = futility_IF,
    total_events = total_events,
    analysis_model = analysis_model,
    target_null_futility_rate = target_null_futility_rate,
    n_sims = n_sims,
    n_cores = n_cores,
    seed = seed,
    package_version = tryCatch(as.character(utils::packageVersion("DTEAssurance")),
                               error = function(e) NA_character_),
    timestamp = as.character(Sys.time())
  )
  
  list(boundary = boundary,
       scenario_futility_rates = scenario_futility_rates,
       raw_Z_by_scenario = raw_Z_by_scenario,
       settings = settings)
}